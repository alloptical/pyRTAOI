''''
 Dependences:
 follow instructions in these links to install all dependences
 caiman  : https://github.com/flatironinstitute/CaImAn
 opencv  (pip install .whl) : https://www.lfd.uci.edu/~gohlke/pythonlibs/#opencv (caiman installs it)
 nidaqwx (pip install)
 pyqtgraph (pip install)

SPYDER HAS NO AUTOSAVE FUNCTION!  (only autosaves when script is run)
CAREFUL WHEN USING 'REPLACE ALL' - IT WILL QUIT, WITHOUT SAVING!!


 TO DO:
0. data stream stopped in the middle of a movie - not polling data from buffer fast enough

    0.1 nested thread - read data and emit signal when the correct number of sample is collected - bad idea
        - saved in pyRTAOI-20180620_1; this guaranteed data streaming from PV but it takes more time to process one frame, 35 ms without nested thread; 60 ms with
    0.2 use queue to communicate between worker and streamer - done
        - saved in pyRTAOI-20180620_2; < 20 ms when running offline; need to test on rig
    0.3 use gpu for motion correction? - caiman motion correction is below 10 ms; good enough

    - down sample by 2 kept up with 30Hz image stream

1. log data for offline anaysis - photostimulation targets and frame idx; save online analysis result - done, need test
    offline analysis - convert cnmf results to .mat file - done

2. deal with zoom - check
3. mark point setting for changing photostim power - done
    control voltage = (outMax/displayMax)*PV value; where outMax and displayMax are defined in prairie view configuration.xml
    need to map power to voltage - done
    send volt to pockels via ni card -  consider outsource it to holoblink?

4. Plot STA dF/F - send triggers - check
5. threshold sometimes becomes Inf - check
6. auto initialisation - take reference movie and load - done
7. photostim protocol - done, need to test on rig
8. delete target by right clicking  - done by left click
9. some memory leak - every time the app is run, memory increases by ~3-4 % - mostly solved by explicit deleting of big variables

10. release gpu memory after one movie - done
11. daq task issues - 'the specified resource is reserved error' - this only happens when triggering online; need to restart spyder after one run

----
12. add option to save out motion-corrected movie by image-saver and by worker (when online photostim is not enabled)
13. save out detected events in ROIs for making trigger-triggered average frames


# other notes:
- numpy.append is slower than list.append, so avoid using it in loops but there's not much difference if just appending a single value
    --> changed np arrays to lists for target info

- pkl files are slower to load than npz


'''

import sys
import random
import os

# Qt
import GUI
from PyQt5.QtCore import Qt,QObject, pyqtSignal, QThread, QTimer, QRectF, QUrl,QPoint, QRect, QSize,QPointF,QSizeF, QSettings, QCoreApplication
from PyQt5.QtWidgets import (QComboBox, QCheckBox, QLineEdit, QSpinBox, QLabel,
                             QDoubleSpinBox, QFileDialog, QApplication, QWidget,
                             QDesktopWidget, QMainWindow, QMessageBox, QTableWidgetItem)
from PyQt5.QtGui import QFont, QColor, QIcon, QPalette, QDesktopServices, QImage, QPixmap, QPainter,QPen,QBrush


# utilities
import gc
from skimage.external import tifffile
import numpy as np
import pyqtgraph as pg
import time
import cv2
from copy import deepcopy
from scipy.sparse import issparse, spdiags, coo_matrix, csc_matrix
import pickle
import logging
import json
from itertools import compress

# prairie link
import pvlink
from pvlink import *
import ctypes
from ctypes import *

# blink
from bLink import bLink

# timer
if sys.platform == 'win32':
    time_ = time.clock
else:
    time_ = time.time

# caiman libraries
import caiman as cm
from caiman_func.initialisation import initialise
from caiman.motion_correction import motion_correct_iteration_fast
from caiman.base.rois import com
from caiman.utils.utils import load_object, save_object
from caiman.base.rois import extract_binary_masks_from_structural_channel as extract_mask
from caiman.paths import caiman_datadir

# plot libs
import matplotlib
matplotlib.use('QT5Agg', force=True) # Ensure using PyQt5 backend
import matplotlib.pyplot as plt
import pylab as pl

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.widgets import Slider

# ni max
# pydaqmx badly documented, use nidaqmx instead
import nidaqmx
import nidaqmx.task as daqtask
from nidaqmx import stream_writers

# buffer queue
from queue import Queue

# power control voltage
from loadPowerFile import get_power_params

# sta
from utils import sta

# configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler('errors.log')  # create a file handler
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')  # create a logging format
handler.setFormatter(formatter)
logger.addHandler(handler)  # add the handlers to the logger
logger.info('Started application')

# Interpret image data as row-major instead of col-major
pg.setConfigOptions(imageAxisOrder = 'row-major')

class CONSTANTS():
    # constants
    PHOTO_NONE = 0
    PHOTO_FIX_FRAMES   = 1
    PHOTO_ABOVE_THRESH = 2
    PHOTO_BELOW_THRESH = 3
    PHOTO_UPON_DETECT  = 4

#%%
class MyThread(QThread):
    def run(self):
        self.exec_()

#%% stream thread
class DataStream(QObject):
    processFrame_signal = pyqtSignal(object)
    stream_finished_signal = pyqtSignal()

    def __init__(self, **kwargs ):
        super(DataStream,self).__init__()

        # get parameters from Main thread
        self.pl = kwargs.get('prairie',[])
        self.framesRecv = 0
        self.FLAG_PV_CONNECTED = p['FLAG_PV_CONNECTED']
        global STOP_JOB_FLAG
        STOP_JOB_FLAG = False

        # movie name
        if self.FLAG_PV_CONNECTED:
            self.movie_name = os.path.basename(self.pl.get_movie_name())
        else:
            self.movie_name = os.path.splitext(p['moviePath'])[0]
        print('movie name:' + self.movie_name)
        print('stream object created')

    def stream(self):
        print('this is stream func')
        p['FLAG_END_LOADING'] = False

 #%% streaming data from PV
        if self.FLAG_PV_CONNECTED and not p['FLAG_OFFLINE']:
            # PYCUDA
            import pycuda.driver as cuda
            from pycuda.tools import make_default_context
            import pycuda.gpuarray as gpuarray
            from pycuda.compiler import SourceModule

            # initialise cuda
            cuda.init()
            self.context = make_default_context()
            print('starting CUDA module')
            mod = SourceModule("""

            __global__ void sample_mean(short* matrix, int pixelsPerLine,
                            int linesPerFrame, int samplesPerPixel, int flipEvenRows, long* result)
            {
                int idx = threadIdx.x + blockDim.x*blockIdx.x;
                int result_idx = 0;
                int col = 0;
                int num_sample = 0;
                int this_value = 0;

                if (idx<pixelsPerLine*linesPerFrame*samplesPerPixel - samplesPerPixel+1){
                    if ((idx - (idx / samplesPerPixel)*samplesPerPixel) == 0){
                        result_idx = idx / samplesPerPixel;
                        col = result_idx - (result_idx / pixelsPerLine)*pixelsPerLine;
                        if ((result_idx / pixelsPerLine) - ((result_idx / pixelsPerLine) / 2) * 2 != flipEvenRows){
                            result_idx = result_idx + pixelsPerLine - 2 * col - 1;
                        }

                        for (int i = 0; i < samplesPerPixel; i++){
                            if (matrix[idx + i]>8192){
                                this_value += matrix[idx + i] - 8192;
                                num_sample += 1;
                            }
                        }

                        if (num_sample>0){ result[result_idx] = this_value / num_sample; }
                    }
                }
            }

            """, cache_dir = False)

            # get cuda func
            sample_mean = mod.get_function('sample_mean')
            print('CUDA module initialised')

            # allocate memory for image buffer
            [samplesPerPixel, pixelsPerLine, linesPerFrame] = self.pl.get_frame_size()
            samplesPerFrame = samplesPerPixel*pixelsPerLine*linesPerFrame
            BUFFER_SIZE = p['pvBufferFrames']*samplesPerFrame
            thisID = os.getpid()
            buf = np.ndarray(BUFFER_SIZE, dtype = c_int16)
            buf_add, read_only_flag = buf.__array_interface__['data']
            _rrd = make_command('-rrd',str(thisID),str(c_int64(buf_add).value),str(BUFFER_SIZE))
            print('cpu buffer created')

            # assign destination buffer in gpu - this is frame size
            dest = np.ndarray(pixelsPerLine*linesPerFrame, dtype = 'int32')
            dest_gpu = gpuarray.to_gpu(dest.astype(np.int32))
            print('gpu initiated')

            # initialise local parameters before start scanning
            self.FLAG_SCANNING  = False
            loop = 0
            self.framesRecv = 0
            recvTime = 0
            incorrectCount = 0
            incorrectLimit = 5000

            # start imaging
            num_samples = self.pl.send_recv(_rrd) # flush
            if not self.pl.send_done(self.pl._ts): # t-series
                self.FLAG_SCANNING = True
                print('started t series')

            while self.FLAG_SCANNING and (not STOP_JOB_FLAG):
                loop += 1
    #                print('This is loop '+ str(loop))
                loop_start_time = time.time()
                num_samples = self.pl.send_recv(_rrd)

                if not num_samples:
                    self.FLAG_SCANNING  = False

                elif int(num_samples) == samplesPerFrame:
                    recvTime += time.time()-loop_start_time
#                    recv_time = time.time()
                    # -----  use cuda -------
                    buf_gpu = gpuarray.to_gpu(buf.astype(np.int16))
                    sample_mean( buf_gpu, pixelsPerLine, linesPerFrame, samplesPerPixel, p['flipEvenRows'], dest_gpu,
                        block=(1024,1,1), grid = (int(samplesPerFrame/1024),1,1) )
                    dest_gpu = dest_gpu.reshape(linesPerFrame,pixelsPerLine)
                    thisFrame = dest_gpu.get()
#                    print('convert time:' + str("%.4f"%(time.time()-recv_time)))

                    # put frame in global queue
#                    qcopy_time = time.time()
                    qbuffer.put(thisFrame.copy().astype(np.float32))
#                    print('frame to queue:' + str("%.4f"%(time.time()-qcopy_time)))
                    self.framesRecv += 1
                    print('frames recv ='+str(self.framesRecv))
                    incorrectCount = 0
                else:
                    incorrectCount +=1

                if incorrectCount>incorrectLimit:
                    print('too many incorrect frames, quit loop')
                    self.FLAG_SCANNING = False

            if (not self.FLAG_SCANNING) or STOP_JOB_FLAG:  # - need a flag outside of this thread to stop it. clicking stop works
                self.context.pop()
                del self.context
                from pycuda.tools import clear_context_caches
                clear_context_caches()

                # free and delete gpu arrays
                del buf_gpu
                del dest_gpu
                del sample_mean
                print('stopped scanning') # gets here
                p['FLAG_END_LOADING'] = True
                self.stream_finished_signal.emit()

#%% offline, read frames from tiff
        elif p['FLAG_OFFLINE']:
            # load movie
            if p['moviePath'] != 'U:/simulate_movie/20170823.tif':
                print('Loading video')
                Y_ = cm.load(p['moviePath'], subindices=slice(0,10000,None))
                print('Video loaded')
            else:
                print('No video provided!')
                return

            for frame_count, frame in enumerate(Y_):        # now process each file
                if np.isnan(np.sum(frame)):
                    raise Exception('Frame ' + str(frame_count) + ' contains nan')
                qbuffer.put(frame.copy().astype(np.float32))
                self.framesRecv += 1

            self.stream_finished_signal.emit()
            print('stream finished')
            p['FLAG_END_LOADING'] = True
        else:
            print('check pv connection - warning from streamer')

    def stop(self):
        print('stream stop')
        self.FLAG_SCANNING = False
        global STOP_JOB_FLAG
        STOP_JOB_FLAG = True

#%% save image thread
class imageSaver(QObject):
    # save image in the global queue from data steam thread to multi page tiff file
    # setup signals
    finished_signal          = pyqtSignal()

    def __init__(self, **kwargs ):
        super(imageSaver,self).__init__()
        self.mptiff_path = kwargs.get('save_movie_path',[])
        self.frameSaved = 0
        p['FLAG_END_LOADING'] = False

    def saveImage(self):
        # buffer for save multipage tiff
        mptiff_path = self.mptiff_path
        print('imageSaver, movie name='+str(mptiff_path))
        with tifffile.TiffWriter(mptiff_path, bigtiff=True) as tif:
            while ((not p['FLAG_END_LOADING']) or (not qbuffer.empty())):
                try:
                    frame_in = qbuffer.get(timeout = 5)
                except:
                    break
                tif.save(frame_in)
                self.frameSaved+=1
                print('number frame saved = '+str(self.frameSaved))

        # finishing
        print('Total number of frames saved '+str(self.frameSaved))
        self.finished_signal.emit()

#%% process thread
class Worker(QObject):
    # setup signals
    status_signal            = pyqtSignal(str,name ='statusSignal')
    roi_signal               = pyqtSignal(object,name = 'roiSignal')
    refreshPlot_signal       = pyqtSignal(object,name = 'refreshPlot')
    thresh_signal            = pyqtSignal(object,name = 'threshSignal')
    sta_amp_signal           = pyqtSignal(object,name = 'staAmpSignal')
    frame_signal             = pyqtSignal(object,name = 'frameSignal')
    refreshScatter_signal    = pyqtSignal(object)
    showROIIdx_signal        = pyqtSignal(name = 'showROIIdx')
    sendTTLTrigger_signal    = pyqtSignal()
    sendCoords_signal        = pyqtSignal(name = 'sendCoordsSignal')
    sendPhotoStimTrig_signal = pyqtSignal()
    getROImask_signal        = pyqtSignal(object, object, name = 'getROImaskSignal')
    transDataToMain_signal   = pyqtSignal(object,object,object,object,name = 'transDataToMain')
    updateTargetROIs_signal  = pyqtSignal()
    finished_signal          = pyqtSignal()


    def __init__(self, **kwargs ):
        super(Worker,self).__init__()

        # get parameters from Main thread
        self.c = kwargs.get('caiman',{})
        self.pl = kwargs.get('prairie',[])
        self.ROIlist = kwargs.get('ROIlist',dict())
        print('worker pl empty? ' + str(self.pl == []))

        # sta traces
        self.sta_traces = np.zeros([p['MaxNumROIs'],p['numberStims'],p['staPreFrame']+p['staPostFrame']])
        self.sta_traces.fill(np.nan)
        p['stimStopFrame'] = p['stimStartFrame']+p['numberStims']*p['interStimInterval']+p['staPreFrame']
        self.stim_frames = np.arange(p['stimStartFrame'],p['stimStopFrame']+1,p['interStimInterval'])

        # count processed frames
        self.framesProc = 0

        # FLAGS
        self.UseONACID = p['UseONACID']
        print('use onacid? ', self.UseONACID)
        self.FLAG_PV_CONNECTED = p['FLAG_PV_CONNECTED']
        self.FLAG_OFFLINE = p['FLAG_OFFLINE']
        print('pv connected? ' + str(self.FLAG_PV_CONNECTED))
        print('offline analysis ' + str(self.FLAG_OFFLINE))

        if self.c == {}:
            self.UseONACID = False
            print('No reference movie provided')
            
        if self.pl == []:
            print('No Prairie object provided')

        # initialise sta
        self.sta_start_frames = self.stim_frames -p['staPreFrame']
        self.photo_stim_frames = self.stim_frames + p['photoWaitFrames']
        self.flag_sta_recording = False
        self.sta_frame_idx = 0
        self.sta_trace_length = p['staPreFrame']+p['staPostFrame']
        self.num_stims = p['numberStims']

        # get current parameters
        self.BufferLength = p['BufferLength']
        self.RoiBuffer = p['RoiBuffer']
        self.BufferPointer = p['BufferPointer']
        self.ROIlist_threshold = np.zeros(p['MaxNumROIs'])
        self.dist_thresh = 21 # reject new rois too close to existing rois
        self.display_shift = False  # option to display motion correction shift

        self.tottime = []

        # lateral motion
        self.shifts = []

        # make Temp folder
        self.save_dir = p['temp_save_dir']

        global STOP_JOB_FLAG
        STOP_JOB_FLAG = False

        p['FLAG_END_LOADING'] = False


#%% process frame
    def work(self):
        # local counts
        sta_stim_idx = 0
        sta_frame_idx = 0
        LastPlot = 0
        framesProc = self.framesProc
        refreshFrame = p['refreshFrame']
        num_photostim = 0
        current_target_idx = []
        last_target_idx = []
        FLAG_TRIG_PHOTOSTIM = False
        FLAG_SEND_COORDS = False

        # local param
        try:
            FLAG_USE_ONACID = self.UseONACID
            FLAG_SAVE_TIFF = p['saveAsTiff']
            BufferPointer = self.BufferPointer
            sta_start_frames = self.sta_start_frames
            photo_stim_frames = self.photo_stim_frames # predefined fixed frames
            stim_frames = self.stim_frames # stim at fixed intervals
    
            # Local buffer for recording online protocol
            online_photo_frames = []
            online_photo_targets = []
        except Exception as e:
            print(e)

        # get power control prameters if provided
        try:
            power_polyfit_p = p['power_polyfit_p']
            photoPowerPerCell = p['photoPowerPerCell']
            NI_UNIT_POWER_ARRAY = p['NI_UNIT_POWER_ARRAY']
        except:
            pass

        if FLAG_USE_ONACID:
            # Use locally scoped variables to speed up
            self.status_signal.emit('Using OnACID')
    
            # Extract initialisation parameters
            mot_corr = True # hard coded now; using motion correction from onacid
            cnm2 = self.c['cnm2'] #deepcopy(cnm_init)
            img_norm = self.c['img_norm'].copy().astype(np.float32)
            img_min = self.c['img_min'].copy().astype(np.float32)
            ds_factor = self.c['ds_factor']
            dims = self.c['cnm_init'].dims
            keep_prev = self.c['keep_prev']
            
            # Define OnACID parameters
            max_shift = np.ceil(10./ds_factor).astype('int')  # max shift allowed
            accepted = [ix+cnm2.gnb for ix in cnm2.accepted]
            rejected_count = cnm2.N - len(accepted)
            expect_components = True
            self.com_count = len(accepted) #self.c['cnm2'].N
            com_count = self.com_count
            
            # Local record of all roi coords
            print('number of initial rois '+str(len(accepted)))
            ROIx = np.asarray([item.get("ROIx") for item in self.ROIlist[:com_count]])   # could also predefine ROIx/y as arrays of bigger size (e.g. 100) and fill them in instad of append
            ROIy = np.asarray([item.get("ROIy") for item in self.ROIlist[:com_count]])
            
            # Keep or discard previous online cells if rerunning OnACID
            keep_full_traces = True  # keep full trace history of the previous session (may be necessary for deconv of last minibatch, otherwise errors possible)
            if keep_prev:
                coms_init = self.c['coms']
                if keep_full_traces:
                    t_cnm = cnm2.t
#                else:
##                    mbs = cnm2.minibatch_shape
#                    init_shape = cnm2.initbatch
#                    t_cnm = init_shape # + mbs
##                    cnm2.C_on[:cnm2.N+cnm2.gnb,init_shape:t_cnm] = cnm2.C_on[:cnm2.N+cnm2.gnb,-mbs:]
#                    cnm2.C_on[:cnm2.N+cnm2.gnb,t_cnm:] = 0  # clears previous results -- should keep at least last minibatch if clearing
            else:
                coms_init = self.c['coms_init']
                t_cnm = cnm2.initbatch
            coms = coms_init.copy()
            
            
            # Extract opsin mask info constructed from c1v1 image
            try:
                opsin_mask = self.c['opsin_mask']
                opsin_thresh = self.c['opsin_thresh']
                opsin = self.c['cnm2'].opsin
                overlap = self.c['overlap']

            except Exception as e:
                opsin_mask = np.array([])
                opsin_positive = True  # default for when no opsin mask
                print(e)

        
        # prepare to save frames to tiff file
        if FLAG_SAVE_TIFF:
            try:
                MyTiffwriter =  tifffile.TiffWriter(p['currentMoviePath'], bigtiff=True)
                print('movie will be saved as '+p['currentMoviePath'])
            except Exception as e:
                print('tiffwriter error:' + e)
                FLAG_SAVE_TIFF = False
                
        
        if p['FLAG_OFFLINE']:
            max_wait = None  # wait till video loaded and buffer starts filling
        else:
            max_wait = 10 # timeout


        # keep processing frames in qbuffer
        while ((not p['FLAG_END_LOADING']) or (not qbuffer.empty())) and not STOP_JOB_FLAG:
            # get data from queue
            try:
                frame_in = qbuffer.get(timeout = max_wait)
            except Exception: # timeout exception
                print('Timeout Exception: empty qbuffer')
                break # this may be a problem for online too? stops reading the stream completely if timeout happens

            try:
                t0 = time.time()
                framesProc = framesProc+1
    
                if FLAG_USE_ONACID:
                    # move trace buffer pointer
                    if BufferPointer==self.BufferLength-1:
                        BufferPointer = 0
                    else:
                        BufferPointer +=1
    
                    # process current frame
                    if ds_factor > 1:
                        frame_in = cv2.resize(frame_in, img_norm.shape[::-1])   # downsampling
    
                    frame_in -= img_min                                       # make data non-negative
    
                    if mot_corr:                                            # motion correct
#                        mot_corr_start = time.time()
                        templ = cnm2.Ab.dot(cnm2.C_on[:cnm2.M, t_cnm - 1]).reshape(cnm2.dims, order='F') * img_norm
                        frame_cor, shift = motion_correct_iteration_fast(frame_in, templ, max_shift, max_shift)
                        self.shifts.append(shift)
    
#                        print('caiman motion correction time:' + str("%.4f"%(time.time()-mot_corr_start)))
                    else:
                        frame_cor = frame_in
    
                    frame_cor = frame_cor / img_norm                            # normalize data-frame
                    cnm2.fit_next(t_cnm, frame_cor.reshape(-1, order='F'))      # run OnACID on this frame
    
                    # detect new compunents
                    if expect_components:
                        update_comp_time = time.time()
                        if cnm2.N - (com_count+rejected_count) == 1:
                            new_coms = com(cnm2.Ab[:, -1], dims[0], dims[1])[0]
    
                            # Check for repeated components
                            close = abs(coms - new_coms) < self.dist_thresh/ds_factor
                            repeated = any(np.all(close,axis=1))
    
                            if repeated == False: # add component to ROI
                                coms = np.vstack((coms, new_coms))
                                y, x = new_coms   # reversed
                                ROIx = np.append(ROIx,x*ds_factor)  # ~0.11 ms // filling in empty array: ~0.07 ms
                                ROIy = np.append(ROIy,y*ds_factor)
    
                                com_count += 1
                                accepted.append(cnm2.N)
                                print('New cell detected (' + str(cnm2.N-rejected_count) + ')')
    
                                # Check cell for c1v1
    #                            tt = time_()
                                if opsin_mask.size:
                                    cell_A = np.array(cnm2.Ab[:,-1].todense())
                                    cell_mask = (np.reshape(cell_A, dims, order='F') > 0).astype('int')
                                    cell_pix = sum(sum(cell_mask == 1))
    
                                    inter = cv2.bitwise_and(opsin_mask, cell_mask)
                                    inter_pix = sum(sum(inter))
                                    cell_overlap = inter_pix/cell_pix
    
                                    overlap.append(cell_overlap)
                                    opsin_positive = cell_overlap > opsin_thresh
                                    opsin.append(opsin_positive)
    #                                cnm2.opsin.append(cell_overlap > opsin_thresh)
    #                            print(time_()-tt)
    
                                # add new ROI to photostim target, if required
                                try:
                                    if p['FLAG_BLINK_CONNECTED'] and p['FLAG_AUTO_ADD_TARGETS']:
                                        if opsin_positive:  # add target only if opsin present
                                            p['currentTargetX'].append(x*ds_factor)
                                            p['currentTargetY'].append(y*ds_factor)
                                            current_target_idx = range(com_count)
                                            num_stim_targets = len(current_target_idx)
                                            self.sendCoords_signal.emit()
                                            self.updateTargetROIs_signal.emit()
                                            FLAG_SEND_COORDS = False
    
                                    self.getROImask_signal.emit(x,y) # add roi coords to list in main
                                    print('add new component time:' + str("%.4f"%(time.time()-update_comp_time)))
                                except Exception as e:
                                    print(e)
    
                                if com_count == p['MaxNumROIs']:
                                    expect_components = False
                                    print('Not accepting more components')
                            else:
                                print('Repeated component found!')
                                rejected_count += 1
    
                    # add data to buffer
                    try:
                        self.RoiBuffer[:com_count, BufferPointer] = cnm2.C_on[accepted,t_cnm] # cnm2.noisyC is without deconvolution
                        self.ROIlist_threshold[:com_count] = np.nanmean(self.RoiBuffer[:com_count,:], axis=1) + 3*np.nanstd(self.RoiBuffer[:com_count,:], axis=1)
                    except Exception as e:
                        print(e)
                        logger.exception(e)
                        print(self.RoiBuffer[:com_count,:])
    
    
    
                # trigger photostim
                if p['photoProtoInx'] == CONSTANTS.PHOTO_FIX_FRAMES and framesProc == photo_stim_frames[num_photostim] and num_photostim < self.num_stims:
                    FLAG_TRIG_PHOTOSTIM = True

                elif p['photoProtoInx'] == CONSTANTS.PHOTO_ABOVE_THRESH: # TODO: no check for opsin here
                    photostim_flag = self.RoiBuffer[:com_count, BufferPointer]-self.ROIlist_threshold[:com_count]
                    last_target_idx = np.copy(current_target_idx)
                    current_target_idx = np.nonzero(photostim_flag>0)
                    num_stim_targets = len(current_target_idx)

                    if (num_stim_targets>0):
                        FLAG_TRIG_PHOTOSTIM = True
                        if np.array_equal(last_target_idx,current_target_idx):
                            FLAG_SEND_COORDS = True

                    p['currentTargetX'] = ROIx[current_target_idx]
                    p['currentTargetY'] = ROIy[current_target_idx]

                elif p['photoProtoInx'] == CONSTANTS.PHOTO_BELOW_THRESH:
                    photostim_flag = self.ROIlist_threshold[:com_count] - self.RoiBuffer[:com_count, BufferPointer]
                    last_target_idx = np.copy(current_target_idx)
                    current_target_idx = np.nonzero(photostim_flag>0)
                    num_stim_targets = len(current_target_idx)

                    if (num_stim_targets>0):
                        FLAG_TRIG_PHOTOSTIM = True
                        if np.array_equal(last_target_idx,current_target_idx):
                            FLAG_SEND_COORDS = True

                    p['currentTargetX'] = ROIx[current_target_idx]
                    p['currentTargetY'] = ROIy[current_target_idx]

                if FLAG_TRIG_PHOTOSTIM:
                    # update phase mask
                    if FLAG_SEND_COORDS:
                        self.sendCoords_signal.emit()

                    # trigger spiral
                    p['NI_2D_ARRAY'][1,:] = NI_UNIT_POWER_ARRAY *np.polyval(power_polyfit_p,photoPowerPerCell*num_stim_targets)
                    self.sendPhotoStimTrig_signal.emit()
                    FLAG_TRIG_PHOTOSTIM = False

                    # take notes
                    online_photo_frames.append(framesProc)
                    online_photo_targets.append(current_target_idx)

                    # update display
                    if FLAG_SEND_COORDS:
                        self.updateTargetROIs_signal.emit()
                        FLAG_SEND_COORDS = False
                    num_photostim +=1

                # Sta recording and sensory stim - not necessary to do online, move it elsewhere
                if sta_stim_idx < self.num_stims:
                    if p['FLAG_STIM_TRIG_ENABLED'] and framesProc == stim_frames[sta_stim_idx]: # send TTL
                        self.sendTTLTrigger_signal.emit()

                    if not self.flag_sta_recording:
                        if framesProc == sta_start_frames[sta_stim_idx]:
                            self.flag_sta_recording = True
                            sta_frame_idx = 0
                    else:
                        self.sta_traces[:com_count,sta_stim_idx,sta_frame_idx] = cnm2.C_on[accepted,t_cnm]
#                        sta_frame_idx +=1
                        if sta_frame_idx == self.sta_trace_length:
                            self.flag_sta_recording = False
                            sta_stim_idx += 1
    
                    # Update display
                    if framesProc > refreshFrame-1: #frame_count>self.BufferLength-1:
                        if LastPlot == refreshFrame:
                            if p['plotOn']:
                                plot_time = time.time()
                                self.refreshPlot_signal.emit(self.RoiBuffer[:com_count,:])
                                print('update plot time = ' +str(time.time()-plot_time))
                            LastPlot = 0
    
                        elif LastPlot == refreshFrame-1:
                            if p['displayOn']:
                                # display current frame
                                if p['denoiseOn']:
#                                    denoise_time = time.time()
                                    A, b = cnm2.Ab[:, cnm2.gnb:], cnm2.Ab[:, :cnm2.gnb].toarray()
                                    C_t, f_t = cnm2.C_on[cnm2.gnb:cnm2.M, t_cnm], cnm2.C_on[:cnm2.gnb, t_cnm]
                                    comps_frame = A.dot(C_t).reshape(cnm2.dims, order = 'F')*img_norm/np.max(img_norm)
                                    bgkrnd_frame = b.dot(f_t).reshape(cnm2.dims, order = 'F')*img_norm/np.max(img_norm)
                                    frame = comps_frame + bgkrnd_frame   # denoised frame = component activity + background
                                    denoised_frame =  np.repeat(frame[:,:,None],3,axis=-1)
                                    denoised_frame = np.minimum((denoised_frame*255.),255).astype('u1')
    
                                    if ds_factor == 1:
                                        self.roi_signal.emit(denoised_frame)
                                    else:
                                        res_denoised_frame = cv2.resize(denoised_frame, (512, 512), interpolation=cv2.INTER_CUBIC)
                                        self.roi_signal.emit(res_denoised_frame)
#                                    print('generate denoise frame time:' + str("%.4f"%(time.time()-denoise_time)))
    
                                else:
                                    if ds_factor == 1:
                                        self.roi_signal.emit(frame_cor)
                                    else:
                                        res_frame_cor = cv2.resize(frame_cor, (512, 512), interpolation=cv2.INTER_CUBIC)
                                        self.roi_signal.emit(res_frame_cor)
    
                              # refresh roi scatter
                                shift_ = [round(shift[0]), round(shift[1])]
                                if shift_ != [0,0]:
                                    self.display_shift_time = time.time()
                                    if self.display_shift == True: #denoiseOn == False:
                                        self.refreshScatter_signal.emit(shift_)
                                        print('display shift time:' + str("%.4f"%(time.time()-self.display_shift_time)))
    
                                self.thresh_signal.emit(self.ROIlist_threshold[:com_count])
                                LastPlot += 1
                        else:
                            LastPlot += 1
                    t_cnm +=1

                # else not using onACID
                else:
                    temp_time = time.time()
                    self.roi_signal.emit(frame_in)
                    print('display time '+str("%.4f"%(time.time()-temp_time)))
            except Exception as e:
                print(e)

            # record timing
            self.tottime.append(time.time() - t0)                       # store time for each frame
#            print('process frame time: ' + str("%.4f"%(time.time()- t0)))
            print('frames proc = ' + str(framesProc))
            self.frame_signal.emit(framesProc)
            qbuffer.task_done()
        # end of while loop


#%% post-loop finishing
        print('finishing work')
        print('number of photostims = ' + str(num_photostim))
        self.status_signal.emit('Mean processing time is ' + str(np.nanmean(self.tottime))[:6] + ' sec.')

        if FLAG_SAVE_TIFF:
            MyTiffwriter.close()
        
        if self.UseONACID:
            self.BufferPointer = BufferPointer
            self.BufferPointer = 0
            # show indices on viewbox
            self.showROIIdx_signal.emit()

            accepted = [idx-1 for idx in accepted] # count from 0
            cnm2.accepted = accepted # for easier access in MainWindow
            self.c['cnm2'].t = t_cnm
            self.c['coms'] = coms
            if opsin_mask.size:
                cnm2.opsin = opsin

            # save results to Temp folder
            self.movie_name = p['currentMoviePath']

            save_dict = dict()
            save_dict['cnm2'] = cnm2  # opsin info a part of cnm struct for now
            save_dict['init_com_count'] = self.com_count # com_count from init file (in case any cells removed from init file)
            save_dict['accepted'] = accepted  # accepted currently stored inside cnm2 as well
            save_dict['t_cnm'] = t_cnm
            save_dict['tottime'] = self.tottime
            save_dict['coms'] = coms
            save_dict['online_photo_frames'] = online_photo_frames
            save_dict['online_photo_targets'] = online_photo_targets
            save_dict['ds_factor'] = ds_factor

            try:
                save_dict['opsin_thresh'] = opsin_thresh
                save_dict['opsin_mask'] = opsin_mask
                save_dict['overlap'] = overlap
                # opsin currently as a property of cnm2
            except: pass

            try:
                print('Saving onacid output')
                save_separately = True # temp flag to save in a new folder inside movie folder
                timestr = time.strftime("%Y%m%d-%H%M%S")  # can add to title (day-time)

                if save_separately:  # TODO: what if pv connected? -- test
                    filename = os.path.basename(self.movie_name)[:-4] + '_OnlineProc.pkl'  # + '_DS_' + str(ds_factor) 
                    movie_folder = os.path.dirname(self.movie_name)
                    save_folder = os.path.join(movie_folder, 'pyrtaoi_results')  # save init result in a separate folder
                    if not os.path.exists(save_folder):
                        os.makedirs(save_folder)
                    save_path = os.path.join(save_folder, filename)
                else:
                    save_path = self.movie_name + '_OnlineProc.pkl'   # save results in the same folder # + '_DS_' + str(ds_factor) 
                print('OnACID result saved as:' + save_path)
                save_object(save_dict, save_path)
            except Exception as e:
                print(e)

            # save sta traces to Temp folder - havent testing in live stream mode
            sta_trial_avg = np.nanmean(self.sta_traces,1)
            if(p['staAvgStopFrame']>p['staAvgStartFrame']):
                sta_trial_avg_amp = np.nanmean(sta_trial_avg[:,p['staAvgStartFrame']+p['staPreFrame']:p['staAvgStopFrame']+p['staPreFrame']],1)
                self.sta_amp_signal.emit(sta_trial_avg_amp[:self.com_count])
            else:
                self.status_signal.emit('Check sta average range')
            save_name = str(self.save_dir) + 'sta_traces.npy'
            print('sta traces saved as: '+save_name)
            np.save(save_name,self.sta_traces)


            # transfer data to main and show traces in plot tab
#            self.c['test'] = 1
            self.transDataToMain_signal.emit(cnm2, accepted, t_cnm, sta_trial_avg)

            # delete big variables
            del self.c
            del ROIx
            del ROIy
#        del self.pl

        # finishing
        self.finished_signal.emit()
        print('worker finished')

#    def photostimTargets(self, coords):

    def stop(self):
        self.status_signal.emit('Stop clicked')
        print('worker stop')
        global STOP_JOB_FLAG
        STOP_JOB_FLAG = True


#%%
class MainWindow(QMainWindow, GUI.Ui_MainWindow,CONSTANTS):
    '''
    The GUI window

    '''

    def __init__(self):

        QMainWindow.__init__(self)
        self.setupUi(self)

        # set up install paths
        self.install_dir = os.getcwd()
        self.presets_dir = os.path.join(self.install_dir, 'Presets')
        if not os.path.exists(self.presets_dir):
            os.makedirs(self.presets_dir)


        # prairie default settings - maybe move to a xml and load from there
        self.PV_IP = '128.40.156.161'  # bruker 2: TCP_IP = '128.40.202.220'
        self.PV_PORT = 1236
        self.flipEvenRows = np.int32(1)
        self.pvBufferFrames = 1 # worked for 25/   -- didn't work larger than 1 20180619
        self.FLAG_PV_CONNECTED = False
        self.pl = []
        self.refZoom = 1.14


        # Blink settings
        self.BLINK_IP = '128.40.156.163'
        self.BLINK_PORT = 8888
        self.bl = []

        # setup image window
        self.ImageWindow.ci.layout.setContentsMargins(0, 0, 0, 0)
        self.ImageWindow.ci.layout.setSpacing(0)
        self.graphicsView = self.ImageWindow.addViewBox(enableMouse=False, invertY = True)
        self.ImageWindow.addItem(self.graphicsView)

        self.scatterbrush = QBrush()
        self.scatterbrush.setStyle(0)

        self.imageItem = pg.ImageItem(border='w', invertY=False)
        self.graphicsView.addItem(self.imageItem)
        self.graphicsView.setAspectLocked(True)  # lock the aspect ratio so pixels are always square
        self.graphicsView.setRange(QRectF(0, 0, 512, 512), padding=0)  # Set initial view bounds

        self.BlankImage = np.zeros(shape = (512,512))
        self.imageItem.setImage(self.BlankImage)
        self.movie_path = '' #'U:/simulate_movie/20170823.tif'
        self.ref_movie_path = '' #'U:/simulate_movie/20170823_ref.tif'
        self.movie = np.atleast_3d(self.BlankImage)
        self.movie = np.swapaxes(self.movie,0,2)   # SWAPPED axes of the movie before. not anymore. should it?
        self.movie_folder = ''

        self.opsin_img_path = 'U:/simulate_movie/20170823_Ch01.tif'
        self.opsin_img = np.array([])
        self.opsin_mask = np.array([])
        self.A_opsin = np.array([])
        self.opsinMaskOn = False
        self.A_loaded = np.array([])
        self.mask_path = ''

        # ROI selection
#        self.drawOn = False
        self.RoiCenter = QPoint(244,244)
        self.RoiRadius = 10
        self.thisROI = pg.CircleROI(self.RoiCenter,self.RoiRadius*2)
        self.removeModeOn = False
        self.removeIdx = []

        # initialise ROI list
        self.InitNumROIs = 3
        self.MaxNumROIs = 30  # desired number of cells
        self.thisROIIdx = 0
        self.ALL_ROI_SELECTED = False
        self.ROIlist = dict()
        self.resetROIlist()

        # mean intensity of ROI
        self.BufferLength = 200
        self.RoiBuffer = np.zeros([self.MaxNumROIs,self.BufferLength])
        self.BufferPointer = 0

        # ROI contours
        self.ROIpen = QPen()
        self.ROIpen.setWidthF = 0.5
        self.ROIcontour_item = pg.ScatterPlotItem() # brush for cells on current frame
        self.ROIcontour_item.setBrush(self.scatterbrush)  # this gives open circles
        self.graphicsView.addItem(self.ROIcontour_item)

        self.TargetPen = pg.mkPen(color = (255,112,75), width = 3, style = Qt.DotLine)
        self.RemovePen = pg.mkPen(color = (255,255,0), width = 3, style = Qt.DotLine)
        self.Targetcontour_item = pg.ScatterPlotItem()
        self.Targetcontour_item.setBrush(self.scatterbrush)
        self.graphicsView.addItem(self.Targetcontour_item)

        # dictionary of ROIs
        self.ROIlist = [dict() for i in range(self.MaxNumROIs)]
        self.Thresh_tableWidget.setColumnCount(4)
        self.Thresh_tableWidget.setHorizontalHeaderLabels(['X','Y','Threshold','STA'])
        self.updateTable()

        # display settings
        p['plotOn'] = True
        p['displayOn'] = True
        p['denoiseOn'] = True
        self.ds_factor = 1

        # display traces
        self.plotItem = self.plotArea_graphicsLayoutWidget.addPlot(antialias=True)
        self.plotArea_graphicsLayoutWidget.setAntialiasing(True)
        self.plotItem.disableAutoRange()
        self.plotItem.setXRange(0,self.BufferLength-1)

        # plot tab
        self.fig = Figure()
        self.figCanvas = FigureCanvas(self.fig)
        self.figCanvas.setFocusPolicy(Qt.ClickFocus) # to enable key press detection
        self.figCanvas.setFocus() # to enable key press detection
        self.plot_gridLayout.addWidget(self.figCanvas)

        # caiman values
        self.c = {}
        self.cnm2 = {}
#        self.proc_cnm = {}

        # buffer for sta traces
        self.sta_traces = np.array([])
        self.sta_amp = None

        # photostim target list
        self.numTargets = 0
        self.bl = []
        self.TargetIdx = [] #np.array([],dtype ='uint16') # indices in ROIlist
        self.TargetX = [] #np.array([]) # all target coords
        self.TargetY = [] # np.array([])
        p['ExtraTargetX'] = [] # np.array([]) #selected targets (not in the ROIlist) -- TODO
        p['ExtraTargetY'] = []  # np.array([])
        p['currentTargetX'] = [] # np.array([]) # keep track of what is currently on SLM
        p['currentTargetY'] = [] # np.array([])
        p['FLAG_AUTO_ADD_TARGETS'] = False
        p['FLAG_BLINK_CONNECTED'] = False
        p['targetScaleFactor']  = 1 # currentZoom/refZoom

        # offline
        self.IsOffline = False

        # initialise FLAGs
        p['FLAG_STIM_TRIG_ENABLED'] = False
        self.photoProtoInx = CONSTANTS.PHOTO_NONE

        # make dir to save data temporally
        p['temp_save_dir'] = os.path.join(os.path.dirname(__file__),'/Temp/')
        if not os.path.exists(p['temp_save_dir']):
            os.makedirs(p['temp_save_dir'])
        p['saveResultPath'] = p['temp_save_dir'] +str('temp.pkl')
        self.saveResultPath_lineEdit.setText(p['saveResultPath'])

        # power control
        self.POWER_CONTROL_READY = False


        # get gui elements
        self.getValues()

        # configuration in GUI
        self.trial_config = {}


        # photostim power control
        try:
            p['power_polyfit_p'] = get_power_params()
            p['power_volts'] = np.polyval( p['power_polyfit_p'], p['photoPowerPerCell']) # initialise with one cell
        except Exception as e:
            print(e)

        # thread and task objects
        self.workerObject = None
        self.streamObject = None
        self.imageSaverObject = None
        self.niStimWriter = None
        self.niPhotoStimWriter = None
        self.niPhotostimFullWriter = None


        # daq config
        NI_SAMPLE_RATE = 20000
        NI_TTL_NUM_SAMPLES = int(0.01*NI_SAMPLE_RATE)
        NI_STIM_NUM_SAMPLES = int(p['photoDuration']*0.001*NI_SAMPLE_RATE)
        p['NI_SAMPLE_RATE'] = NI_SAMPLE_RATE
        p['NI_TTL_NUM_SAMPLES'] = NI_TTL_NUM_SAMPLES

        self.daq_array = np.ones((NI_TTL_NUM_SAMPLES-1), dtype = np.float)*5 # to write to single output channel
        self.daq_array = np.append(self.daq_array,0)
        self.niStimTask = daqtask.Task()
        self.niPhotoStimTask = daqtask.Task() #  TTL trigger only
        self.niPhotoStimFullTask = daqtask.Task() # TTL trigger plus power control
        p['NI_1D_ARRAY'] = self.daq_array

        try: # sensory stim
            self.niStimTask.ao_channels.add_ao_voltage_chan(p['stimDaqDevice'])
            self.niStimWriter = stream_writers.AnalogSingleChannelWriter(self.niStimTask.out_stream,True)
            self.niStimWriter.write_one_sample(0,10.0)
            print('sensory stim channel created')
        except Exception as e:
            print(str(e))
            logger.exception(e)
            self.updateStatusBar(str(e))
            time.sleep(2)

        try: # photo stim
            # single-channel writer - only send triggers
            self.niPhotoStimTask.ao_channels.add_ao_voltage_chan('Dev5/ao2','photostim_simple_trig') # hard coded for now
            self.niPhotoStimWriter= stream_writers.AnalogSingleChannelWriter(self.niPhotoStimTask.out_stream,True)
            self.niPhotoStimWriter.write_one_sample(0,10.0)
            print('single channel photostim trigger set')

            # multi-channel writer - send triggers and voltage to power control device  - need test
            self.niPhotoStimFullTask.ao_channels.add_ao_voltage_chan(p['photostimDaqDevice'],'photostim_trig')
            self.niPhotoStimFullTask.ao_channels.add_ao_voltage_chan(p['powerDaqDevice'],'photostim_power')
            self.niPhotostimFullWriter = stream_writers.AnalogMultiChannelWriter(self.niPhotoStimFullTask.out_stream,auto_start = True)
            self.niPhotoStimFullTask.timing.cfg_samp_clk_timing(rate = NI_SAMPLE_RATE,sample_mode= nidaqmx.constants.AcquisitionType.FINITE, samps_per_chan = max(NI_TTL_NUM_SAMPLES+1,NI_STIM_NUM_SAMPLES+1))


            init_output = np.zeros([2, max(NI_TTL_NUM_SAMPLES+1,NI_STIM_NUM_SAMPLES+1)])
            init_output.ravel()[:len(self.daq_array)] = self.daq_array  # use ravel to data between array (fastest way)
            init_output[1,:-1] = 1
            p['NI_2D_ARRAY'] = np.copy(init_output)
            p['NI_UNIT_POWER_ARRAY'] = np.copy(init_output[1,:])
            p['NI_2D_ARRAY'][1,:] = p['NI_UNIT_POWER_ARRAY'] *np.polyval(p['power_polyfit_p'],p['photoPowerPerCell']) # single cell power
            self.niPhotostimFullWriter.write_many_sample(init_output)
            while(not self.niPhotoStimFullTask.is_task_done()):
                pass
            self.niPhotoStimFullTask.stop()
            self.niPhotostimFullWriter = stream_writers.AnalogMultiChannelWriter(self.niPhotoStimFullTask.out_stream,auto_start = True)
            print('multi channel photostim trigger set')
            self.POWER_CONTROL_READY = True

        except Exception as e:
            print('power control initialisation error')
            print(str(e))
            self.updateStatusBar(str(e))


        # make threads
        self.workerThread = MyThread()
        self.streamThread = MyThread()
        self.imageSaverThread = MyThread()

        # flag reading/streaming data
        p['FLAG_END_LOADING'] = False

        # signal/slot connections
        self.setConnects()

        # methods for elements within MainWindow
        self.ImageWindow.enterEvent = self.displayOpsinImg.__get__(self.ImageWindow)  # ImageWindow when opsinMaskOn
        self.ImageWindow.leaveEvent = self.displayOpsinMask.__get__(self.ImageWindow)
        
        self.figCanvas.keyPressEvent = self.arrow_key_image_control.__get__(self.figCanvas) # changing cell display in the plot tab with arrows


    def random_color(self):
        r = random.randrange(0, 255)
        g = random.randrange(0, 255)
        b = random.randrange(0, 255)
        return QColor(r, g, b)

    def initialiseROI(self):
        self.ALL_ROI_SELECTED = False
        self.thisROIIdx = 0
#        self.MaxNumROIs = self.MaxNumROIs_spinBox.value()
        NumROIs = p['MaxNumROIs']
        print('Number of ROIs to be tracked: ' + str(NumROIs))
        self.RoiBuffer = np.zeros([NumROIs,self.BufferLength])

        # dictionary of ROIs
        self.resetROIlist()

        # threshold table in GUI
        self.updateTable()

        # clear scatter plot
        self.ROIcontour_item.clear()

    def updateTable(self):
        self.Thresh_tableWidget.clear()  # empty the table
        self.Thresh_tableWidget.setHorizontalHeaderLabels(['X','Y','Threshold','STA'])  # column names disappear with clearing

        NumROIs = self.MaxNumROIs
        self.Thresh_tableWidget.setRowCount(NumROIs)
        self.Thresh_labels = [QTableWidgetItem() for i in range(NumROIs)]
        self.STA_labels = [QTableWidgetItem() for i in range(NumROIs)]
        self.X_labels = [QTableWidgetItem() for i in range(NumROIs)]
        self.Y_labels = [QTableWidgetItem() for i in range(NumROIs)]

    def showMousePosition(self,event):
        x = event.x()
        y = event.y()
        self.xPos_label.setText(str(x))
        self.yPos_label.setText(str(y))

    def getMousePosition(self,event):
        pointSize = self.RoiRadius*2+5
        x = event.pos().x()
        y = event.pos().y()

        if self.removeModeOn:
            det_dist = 10

            for idx in range(self.thisROIIdx):
                ROI_x = self.ROIlist[idx]["ROIx"]
                ROI_y = self.ROIlist[idx]["ROIy"]
                detected = abs(x - ROI_x) <= det_dist and abs(y - ROI_y) <= det_dist

                if detected:
                    if idx in self.removeIdx:
                        # reset all remove list works
#                            self.Targetcontour_item.clear()
#                            self.updateTargets()
#                            self.removeIdx = []

                        # unselect individual remove items
                        selected = self.Targetcontour_item.data
                        selected_xy = [[value[0], value[1]] for value in selected][self.numTargets:]

                        ROI_xy = [ROI_x, ROI_y]
                        selected_idx = selected_xy.index(ROI_xy)

                        del selected_xy[selected_idx]
                        remove_x = [item[0] for item in selected_xy]
                        remove_y = [item[1] for item in selected_xy]

                        self.Targetcontour_item.clear()
                        self.updateTargets()
                        self.Targetcontour_item.addPoints(x = remove_x, y = remove_y, pen = self.RemovePen, size = pointSize)

                        self.removeIdx.remove(idx)
                    else:
                        self.removeIdx.append(idx)
                        self.Targetcontour_item.addPoints(x = [ROI_x], y = [ROI_y], pen = self.RemovePen, size = pointSize)

        else:
            print(str(x)+' '+str(y))
            det_dist = 20  # detection distance
            detected = 0

            # check caiman ROIs - if select all, can deselect some and keep the rest
            if self.selectAll_checkBox.isChecked():
                for idx in range(len(self.TargetIdx)):
                    detected = abs(x - self.TargetX[idx]) <= det_dist and abs(y - self.TargetY[idx]) <= det_dist
                    if detected:
                        del self.TargetIdx[idx]
                        self.updateTargets()
                        return

            # check extra targets
            for idx in range(len(p['ExtraTargetX'])):
                detected = abs(x - p['ExtraTargetX'][idx]) <= det_dist and abs(y - p['ExtraTargetY'][idx]) <= det_dist
                if detected:
                    del p['ExtraTargetX'][idx]
                    del p['ExtraTargetY'][idx]
                    self.updateTargets()
                    return

            if not detected:
                p['ExtraTargetX'].append(x)
                p['ExtraTargetY'].append(y)
                print('selected x = ' + str(x))
                print('selected y = ' + str(y))
                self.updateTargets()


    def displayOpsinImg(self,event):
        if self.opsinMaskOn:
            self.updateImage(cv2.resize(np.squeeze(self.opsin_img), (512, 512), interpolation=cv2.INTER_CUBIC))

    def displayOpsinMask(self,event):
        if self.opsinMaskOn:
            self.updateImage(cv2.resize(np.squeeze(self.opsin_mask).astype('u1'), (512, 512), interpolation=cv2.INTER_CUBIC))

    def tempTest(self):
        print('button clicked')
        self.plotSTAonMasks()

    def setConnects(self):
        # load and save
        self.loadMoviePath_pushButton.clicked.connect(self.loadMoviePath)
        self.loadRefMoviePath_pushButton.clicked.connect(self.loadRefMoviePath)
        self.loadOpsinImgPath_pushButton.clicked.connect(self.loadOpsinImg)
        self.createMask_pushButton.clicked.connect(self.createOpsinMask)
        self.showCellsOnMask_pushButton.clicked.connect(self.showCellsOnMask)
        self.loadMask_pushButton.clicked.connect(self.loadMask)
        self.initialise_pushButton.clicked.connect(self.initialiseCaiman)

        self.saveConfig_pushButton.clicked.connect(self.saveConfig)
        self.loadConfig_pushButton.clicked.connect(self.loadConfig)
        self.saveResult_pushButton.clicked.connect(self.saveResults)
        self.browseResultPath_pushButton.clicked.connect(self.browseResultPath)

        self.removeMode_pushButton.clicked.connect(self.removeModeController)

        # start worker
        self.run_pushButton.clicked.connect(self.clickRun)

        # display
        self.plotOn_checkBox.clicked.connect(self.switch_plotOn)
        self.displayOn_checkBox.clicked.connect(self.switch_displayOn)
        self.denoiseOn_checkBox.clicked.connect(self.switch_denoiseOn)
        self.plotSTA_pushButton.clicked.connect(self.plotSTA)
        self.plotSTAonMasks_pushButton.clicked.connect(lambda: self.plotSTAonMasks(self.sta_amp))
        self.showComponents_pushButton.clicked.connect(lambda: self.plotSTAonMasks(None))
        self.showFOV_pushButton.clicked.connect(self.showFOV)

        # triggers
        self.enableStimTrigger_checkBox.clicked.connect(self.enableStimTrigger)
        self.testTTLTrigger_pushButton.clicked.connect(self.testTTLTrigger)
        self.testPhotoStimTrigger_pushButton.clicked.connect(self.testPhotoStimTrigger)

        # running mode
        self.IsOffline_radioButton.toggled.connect(self.switch_IsOffline)
        self.UseOpsinMask_radioButton.toggled.connect(self.switch_useOpsinMask)
        self.UseOnacidMask_radioButton.toggled.connect(self.switch_useOnacidMask)
        self.photoProto_comboBox.currentIndexChanged.connect(self.photoProtoChanged)

        # connections
        self.connectPV_pushButton.clicked.connect(self.connectPV)
        self.connectBlink_pushButton.clicked.connect(self.connectBlink)
        self.SendCoords_pushButton.clicked.connect(self.sendCoords)

        # targets
        self.addNewROIsToTarget_checkBox.clicked.connect(self.autoAddClicked)
        self.selectAll_checkBox.clicked.connect(self.selectAllROIsClicked)
        self.clearTargets_pushButton.clicked.connect(self.clearTargets)
        self.updateTargets_pushButton.clicked.connect(self.updateTargets)

        # mouse events
        self.graphicsView.scene().sigMouseClicked.connect(self.getMousePosition)
        self.graphicsView.scene().sigMouseMoved.connect(self.showMousePosition)

        # reference movie
        self.takeRefMovie_pushButton.clicked.connect(self.takeRefMovie)

        # others
        self.test_pushButton.clicked.connect(self.tempTest)
        self.reset_pushButton.clicked.connect(self.resetAll)

        # auto add connects to update p and trial config plot whenever anything changes
        widgets = (QComboBox, QCheckBox, QLineEdit, QSpinBox, QDoubleSpinBox)
        for obj in self.findChildren(widgets):
            if isinstance(obj, QComboBox):
                obj.currentIndexChanged.connect(self.getValues)
            if isinstance(obj, QCheckBox):
                obj.stateChanged.connect(self.getValues)
            if isinstance(obj, QLineEdit):
                obj.textChanged.connect(self.getValues)
            if isinstance(obj, QSpinBox):
                obj.valueChanged.connect(self.getValues)
            if isinstance(obj, QDoubleSpinBox):
                obj.valueChanged.connect(self.getValues)

    def getValues(self):
        # extract gui values store in self.p
        widgets = (QComboBox, QCheckBox, QLineEdit, QSpinBox, QDoubleSpinBox)

        for obj in self.findChildren(widgets):
            fullname = str(obj.objectName())
            trimmed_name = fullname.split('_')[0]
            if isinstance(obj, QComboBox):
                p[trimmed_name] = str(obj.currentText())
            if isinstance(obj, QCheckBox):
                p[trimmed_name] = bool(obj.isChecked())
            if isinstance(obj, QLineEdit):
                if 'spinbox' not in fullname:
                    p[trimmed_name] = str(obj.text())
            if isinstance(obj, QSpinBox):
                p[trimmed_name] = int(obj.value())
            if isinstance(obj, QDoubleSpinBox):
                p[trimmed_name] = float(obj.value())
            if isinstance(obj, QLabel):
                p[trimmed_name] = float(obj.text())

        # parameters not in widgets
        p['BufferLength'] = self.BufferLength
        p['RoiBuffer'] = self.RoiBuffer
        p['BufferPointer'] = self.BufferPointer
        p['FLAG_PV_CONNECTED'] = self.FLAG_PV_CONNECTED
        p['FLAG_OFFLINE'] = self.IsOffline
        p['photoProtoInx'] = self.photoProtoInx

#        self.ds_factor = p['dsFactor']
        self.MaxNumROIs = p['MaxNumROIs']
        self.flipRowChanged()
        self.currentZoomChanged()

        if self.POWER_CONTROL_READY:
            self.updatePowerControl()


    def updatePowerControl(self):

        NI_TTL_NUM_SAMPLES = p['NI_TTL_NUM_SAMPLES']
        NI_STIM_NUM_SAMPLES = int(p['photoDuration']*0.001*p['NI_SAMPLE_RATE'])
        init_output = np.zeros([2, max(NI_TTL_NUM_SAMPLES+1,NI_STIM_NUM_SAMPLES+1)])
        init_output.ravel()[:len(p['NI_1D_ARRAY'])] = p['NI_1D_ARRAY']
        init_output[1,:-1] = 1
        p['NI_2D_ARRAY'] = np.copy(init_output)
        p['NI_UNIT_POWER_ARRAY'] = np.copy(init_output[1,:])
        p['NI_2D_ARRAY'][1,:] = p['NI_UNIT_POWER_ARRAY'] *np.polyval(p['power_polyfit_p'],p['photoPowerPerCell'])*10

        # reconfigure sample number on task
        self.niPhotoStimFullTask.timing.cfg_samp_clk_timing(rate = p['NI_SAMPLE_RATE'],sample_mode= nidaqmx.constants.AcquisitionType.FINITE, samps_per_chan = max(NI_TTL_NUM_SAMPLES+1,NI_STIM_NUM_SAMPLES+1))

    def testTTLTrigger(self):
        try:
            self.sendTTLTrigger()
            self.updateStatusBar('All settings updated, test trigger sent')
        except Exception as e:
            print(e)

    def saveResults(self):  # TODO: add coms?
        try:
            save_dict = dict()
            save_dict['cnm2'] = self.c['cnm2'] #self.cnm2
            save_dict['com_count_init'] = self.K_init
            save_dict['accepted'] = self.accepted
            save_dict['t_cnm'] = self.t_cnm
            save_object(save_dict, p['saveResultPath'])
        except Exception as e:
            self.updateStatusBar('Save results error: '+str(e))

    def browseResultPath(self):
        saveResultsPath = QFileDialog.getSaveFileName (self, 'Select save path', p['saveResultPath'], '.pkl')
        p['saveResultPath'] = saveResultsPath[0]+saveResultsPath[1]
        self.saveResultPath_lineEdit.setText(p['saveResultPath'])

    def getPathFromPV(self):
        if p['FLAG_PV_CONNECTED']:
            p['saveResultPath'] = self.pl.get_movie_name()
        self.saveResultPath_lineEdit.setText(p['saveResultPath'])

    def testPhotoStimTrigger(self):
        try:
            self.sendPhotoStimTrigger()
            self.updateStatusBar('All settings updated, test trigger sent')
        except Exception as e:
            print(e)

    def sendTTLTrigger(self):
        numsent = self.niStimWriter.write_many_sample(self.daq_array,10.0)
        print('stim trigger sent')
        return numsent

    def sendPhotoStimTrigger(self):
        try:
            if not p['SendPowerVolt']:
               self.niPhotoStimWriter.write_many_sample( p['NI_1D_ARRAY'],10.0)
            else:
                self.niPhotostimFullWriter.write_many_sample(p['NI_2D_ARRAY'],10.0)
                while(not self.niPhotoStimFullTask.is_task_done()):
                    pass
                self.niPhotoStimFullTask.stop()
                self.niPhotostimFullWriter = stream_writers.AnalogMultiChannelWriter(self.niPhotoStimFullTask.out_stream,auto_start = True)
            print('photostim trigger sent')
        except Exception as e:
            print('send photostim trigger error')
            print(e)

    def clearTargets(self):
        msg = QMessageBox()
        msg.setText("Clear all targets?")
        msg.setWindowTitle('pyRTAOI Message')
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        retval = msg.exec_()
        if(retval==QMessageBox.Ok):
            self.Targetcontour_item.clear()
            self.TargetIdx = [] # indices in ROIlist
            self.TargetX = [] # all target coords
            self.TargetY = []
            p['ExtraTargetX'] = [] #selected targets (not in the ROIlist) -- TODO
            p['ExtraTargetY'] = [] 
            p['currentTargetX'] = [] # keep track of what is currently on SLM
            p['currentTargetY'] = []

            self.numTargets = 0
            self.selectAll_checkBox.setChecked(False)

    def deleteTextItems(self):
        for obj in self.graphicsView.addedItems[:]:
            if isinstance(obj,pg.TextItem):
                self.graphicsView.removeItem(obj)

    def resetROIlist(self):
        self.ROIlist = [dict() for i in range(self.MaxNumROIs)]
        for ROIIdx in range(0,self.MaxNumROIs):
            self.ROIlist[ROIIdx]["ROIx"]= list()
            self.ROIlist[ROIIdx]["ROIy"] = list()
            self.ROIlist[ROIIdx]["ROIcolor"] = self.random_color()
            self.ROIlist[ROIIdx]["threshold"] = list()
            self.ROIlist[ROIIdx]["STA"] = list()
        # print(self.ROIlist[ROIIdx]["ROIcolor"])

    def resetAll(self):
        msg = QMessageBox()
        msg.setText("RESET ALL??")
        msg.setWindowTitle('pyRTAOI Message')
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        retval = msg.exec_()
        if(retval==QMessageBox.Ok):
            # update settings
            self.getValues()

            # stim targets
            self.clearTargets()

            # detected ROIs
            self.ROIcontour_item.clear()
            self.thisROIIdx = 0
            self.resetROIlist()
            self.deleteTextItems()
            
            self.IsOffline_radioButton.setChecked(False)
            self.IsOffline_radioButton.setEnabled(False)

            # buffer for sta traces
            self.sta_traces = np.array([])

            # clear table
            self.updateTable()
            
            # clear plots and images
            self.plotItem.clear()
            self.resetMask()
            self.resetFigure()

            self.updateFrameInfo(0)
            self.xPos_label.setText(' ')
            self.yPos_label.setText(' ')
            self.updateStatusBar('')

        msg.setText("RESET CAIMAN??")
        msg.setWindowTitle('pyRTAOI Message')
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
        retval = msg.exec_()
        if(retval==QMessageBox.Ok):
            self.c = {}
            self.cnm2 = {}
#            self.UseONACID_checkBox.setEnabled(False)
            self.UseONACID_checkBox.setChecked(False)
            self.imageItem.setImage(self.BlankImage)


    def connectBlink(self):
        self.updateStatusBar('connecting Blink')
        self.bl = bLink(self.BLINK_IP,self.BLINK_PORT)
        print('bl created')
        if self.bl.CONNECTED:
            p['FLAG_BLINK_CONNECTED'] = True
        print('bl connected = '+str(self.bl.CONNECTED))


    def sendCoords(self):
        print('sending coords')
        if p['FLAG_BLINK_CONNECTED']:
            # scale targets coordinates (to refZoom with which the SLM transform matrix was computed);
            # currentTargetXY should be coordinates in 512x512 frames
            currentTargetX = [int((item-255.0)*p['targetScaleFactor']+255.0) for item in p['currentTargetX']]
            currentTargetY = [int((item-255.0)*p['targetScaleFactor']+255.0) for item in p['currentTargetY']]
            if(self.bl.CONNECTED):
                if not self.bl.send_coords(currentTargetX, currentTargetY):
                    print('Phase mask updated')
                else:
                    print('Update phase mask ERROR!')
                    p['FLAG_BLINK_CONNECTED'] = False

            else:
                self.updateStatusBar('Send coords faild, check blink connection')
                p['FLAG_BLINK_CONNECTED'] = False


    def selectAllROIs(self):
        if self.selectAll_checkBox.isChecked():

            # add all ROIlist to the targetlist
            self.TargetIdx = list(range(0,self.thisROIIdx))
            self.TargetIdx = list(sorted(set(self.TargetIdx)))

            opsin_only = True  # select only opsin expressing ROIs
            if opsin_only and self.A_opsin.size:
                self.TargetIdx = list(compress(self.TargetIdx, self.c['cnm2'].opsin))

            self.TargetX = [self.ROIlist[i]["ROIx"] for i in self.TargetIdx]
            self.TargetY = [self.ROIlist[i]["ROIy"] for i in self.TargetIdx]

            # check for existing extra targets
            det_dist = 10 # remove extra targets only if ROI very close
            idx_remove = []

            for roi_idx in range(len(self.TargetX)):
                for extra_idx in range(len(p['ExtraTargetX'])):
                    x = p['ExtraTargetX'][extra_idx]
                    y = p['ExtraTargetY'][extra_idx]

                    detected = abs(x - self.TargetX[roi_idx]) <= det_dist and abs(y - self.TargetY[roi_idx]) <= det_dist
                    if detected:
                        print('A ROI overlaps with an extra target. Removing the extra duplicate.')
                        idx_remove.append(extra_idx)

            for idx in sorted(idx_remove, reverse=True):
                del p['ExtraTargetX'][idx]
                del p['ExtraTargetY'][idx]

            self.updateTargets()

        else:
            self.TargetIdx = []
            self.TargetX = []
            self.TargetY = []
            self.updateTargets()

    def updateTargets(self):
        # update cell targets based on TargetIdx
        self.TargetX = [self.ROIlist[i]["ROIx"] for i in self.TargetIdx]
        self.TargetY = [self.ROIlist[i]["ROIy"] for i in self.TargetIdx]

        # save current targets to p
        p['currentTargetX'] = self.TargetX + p['ExtraTargetX']
        p['currentTargetY'] = self.TargetY + p['ExtraTargetY']
        self.numTargets = len(self.TargetX) + len(p['ExtraTargetX'])
        if not self.removeModeOn:
            print('Number of current targets: ', self.numTargets)

        self.updateTargetROIs()

    def updateTargetROIs(self):
        # redraw targeted rois in imaging window - to do: draw what is saved in p
        self.Targetcontour_item.clear()
        self.Targetcontour_item.addPoints(x = p['currentTargetX'], y = p['currentTargetY'], pen = self.TargetPen, size = self.RoiRadius*2+5)

    def removeModeController(self):
        self.removeModeOn = not self.removeModeOn

        if self.removeModeOn:
            self.startRemoveMode()
        else:
            self.exitRemoveMode()

    def startRemoveMode(self):
        print('Remove mode on!')

        disable = [self.Prairie_groupBox,self.caiman_groupBox,
                   self.Blink_groupBox, self.opsinMask_groupBox,
                   self.StimOptions_groupBox, self.OfflineAnalysis_groupBox,
                   self.DisplayOptions_groupBox, self.photostim_groupBox,
                   self.groupBox, self.config_groupBox, self.run_groupBox]

        for item in disable:
            item.setEnabled(False)

        self.removeMode_pushButton.setEnabled(True)
        self.removeMode_pushButton.setText('Remove selected')

    def exitRemoveMode(self):
        if self.removeIdx:
            msg = QMessageBox()
            msg.setText("Do you want to remove selected cells?")
            msg.setWindowTitle('pyRTAOI Message')
            msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            retval = msg.exec_()
            if(retval==QMessageBox.Ok):
                self.updateInitialisation()
            else:
                self.updateTargetROIs()
                self.removeIdx = []

        self.removeMode_pushButton.setText('Start remove')

        enable = [self.Prairie_groupBox, self.caiman_groupBox,
           self.Blink_groupBox, self.opsinMask_groupBox,
           self.StimOptions_groupBox, self.OfflineAnalysis_groupBox,
           self.DisplayOptions_groupBox, self.photostim_groupBox,
           self.groupBox, self.config_groupBox, self.run_groupBox]

        for item in enable:
            item.setEnabled(True)

        print('Remove mode off')


    def connectPV(self):
        self.updateStatusBar('connecting PV')
        self.pl = pvlink(self.PV_IP,self.PV_PORT)
        print('pl created')
        ERROR = self.pl.send_done('-sam','Resonant Galvo')
        if(not ERROR):
            self.pl.init_scan()
            self.updateStatusBar('PV connected')
            self.FLAG_PV_CONNECTED = True
        else:
            self.updateStatusBar('connecting to PV failed')

    def flipRowChanged(self):
        if p['flipRow'] is 'Even':
            self.flipEvenRows = np.int32(1)
            print('flip even rows')
        elif p['flipRow'] is 'Odd':
            self.flipEvenRows = np.int32(0)
            print('flip odd rows')
        p['flipEvenRows'] = self.flipEvenRows

    def currentZoomChanged(self):
        p['currentZoom'] = float(p['currentZoom'])
        p['targetScaleFactor'] = self.refZoom/p['currentZoom']

    def photoProtoChanged(self):
        self.photoProtoInx = self.photoProto_comboBox.currentIndex()
        self.updateStatusBar('Photostim protocol changed to' + str(self.photoProtoInx))
        p['photoProtoInx'] = self.photoProtoInx

    def loadMoviePath(self):
#        movie_folder = os.path.dirname(os.path.dirname(self.ref_movie_path))
        movie_path = str(QFileDialog.getOpenFileName(self, 'Load movie', self.movie_folder, 'MPTIFF (*.tif);;All files (*.*)')[0])  # changed '' (default) to movie_folder
        self.movie_path = movie_path
        p['moviePath'] = movie_path
        if movie_path:
            self.IsOffline = True
            self.IsOffline_radioButton.setEnabled(True)
            self.IsOffline_radioButton.setChecked(True)
            self.moviePath_lineEdit.setText(movie_path)
            movie_ext = os.path.splitext(p['moviePath'])[1]
            if movie_ext == '.tif':
                print('Correct video found')

    def showFOV(self):
        try:
            self.imageItem.setImage(cv2.resize(self.c['img_norm'],(512,512),interpolation=cv2.INTER_CUBIC)) # display FOV
        except Exception as e:
            print(e)


    def plotSTAonMasks(self,sta_amp):
        # show cells detected by caiman
        # make a image with sta levels
        self.opsinMaskOn = False

#        cnm = self.proc_cnm
        cnm = self.c['cnm2']

        A, b = cnm.Ab[:, cnm.gnb:], cnm.Ab[:, :cnm.gnb].toarray()

        if issparse(A):
            A = np.array(A.todense())
        else:
            A = np.array(A)

        d, nr = np.shape(A)

        # do not show rejected cells
        nr = self.c['cnm2'].N
        print('N',nr)
        accepted = self.c['cnm2'].accepted

        # use sta value, otherwise use one
        if sta_amp is None: # will show scaled amplitude of A
            sta_amp = np.ones((nr,))*255
        else: # normalise within component before multiply with sta
            for i in accepted: # range(nr):
                A[:,i] = A[:,i]/sum(A[:,i])

        # put value into mask
        cellMask = np.zeros((cnm.dims[1]*cnm.dims[0],))

        j = 0 # separate incrementer for sta_amp (all sta_amp traces are accepted)

        for i in accepted: # range(np.minimum(len(sta_amp),nr)):
            if not np.isnan(sta_amp[j]):
                cellMask+=A[:,i].flatten()*sta_amp[j]
                j += 1

        cellMask2D = np.reshape(cellMask,cnm.dims,order='F')
        cellMask2D = cellMask2D/max(cellMask)*255
        print(cellMask2D.shape)

        norm = plt.Normalize(0,1)
        im = plt.imshow(norm(cellMask2D),aspect = 'equal',cmap = 'Greys')
        plt.colorbar(im, orientation='horizontal')
        plt.show()

        cellMask2D = np.repeat(cellMask2D[:,:,None],3,axis=-1)
        self.imageItem.setImage(cv2.resize(cellMask2D, (512, 512), interpolation=cv2.INTER_CUBIC))


    def plotSTA(self):
        # load file in Temp folder
        try:
            dirname = dirname = os.path.dirname(__file__)
            print(dirname)
            filename = os.path.join(dirname,'/Temp/sta_traces.npy')
            sta_traces = np.load(filename)
        except:
            self.updateStatusBar('STA file not found')
            try:
                filename = str(QFileDialog.getOpenFileName(self, 'Select STA traces file', '', 'npy (*.npy);;All files (*.*)')[0])
            except:
                return
        sta_trial_avg = np.nanmean(sta_traces,1)
        sta_trial_avg_amp = np.nanmean(sta_trial_avg[:,p['staAvgStartFrame']+p['staPreFrame']:p['staAvgStopFrame']+p['staPreFrame']],1)


        num_rois = sta_traces.shape[0]
        num_trials = sta_traces.shape[1]
        plot_rows = np.ceil(np.sqrt(num_rois))
        plot_cols = np.ceil(30/plot_rows)


        fig,axs = plt.subplots(int(plot_rows),int(plot_cols),
                               facecolor='w', edgecolor='k', frameon=False)
        axs = axs.ravel()
        for i in range(num_rois):
            for j in range(num_trials):
                axs[i].plot(sta_traces[i,j,:],color = (.5,.5,.5))
            axs[i].plot(sta_trial_avg[i,:], color=(0,0,0))
            axs[i].set_title('roi'+str(i)+' amp = '+str("%.2f"%sta_trial_avg_amp[i]))
            axs[i].set_aspect('auto')
            axs[i].set_frame_on(False)
            axs[i].set_adjustable('box')

        self.sta_amp = sta_trial_avg_amp
        print(self.sta_amp.shape)

#        plt.show()
        #  add figure to layout
        self.staTracesFig = fig
        self.staTracesCanvas = FigureCanvas(self.staTracesFig)
        self.plot_gridLayout.addWidget(self.staTracesCanvas)
        self.updateStatusBar('STA traces plotted on plot tab')


    def resetFigure(self):
        plt.close('all')

        self.ax1.cla()
        self.ax2.cla()
        self.ax3.cla()
        self.axcomp.cla()

        self.figCanvas.draw()
#            self.figCanvas.update()
        self.figCanvas.flush_events()

        g = gc.collect()
        print('gc: ', g)


    def transDataToMain(self, cnm_struct, accepted, t, sta_amp, plot=True):
#        print('test')  # TODO: changing c inside worker automatically changes it inside mainwindow -- this is sent/shared somehow
#        print(self.c['test'])
#        self.cnm2 = cnm_struct
        self.c['cnm2'] = cnm_struct
        self.accepted = accepted
        self.t_cnm = t
        self.sta_amp = sta_amp

        if self.A_opsin.size:
            print('Checking opsin overlap')
            self.checkOpsin()

        if plot:
            self.plotOnacidTraces(t)


    def plotOnacidTraces(self, t):

        try: self.resetFigure()
        except: pass

        cnm_struct = self.c['cnm2']
        img = self.c['Cn_init']

        C, f = cnm_struct.C_on[cnm_struct.gnb:cnm_struct.M], cnm_struct.C_on[:cnm_struct.gnb]
        A, b = cnm_struct.Ab[:, cnm_struct.gnb:cnm_struct.M], cnm_struct.Ab[:, :cnm_struct.gnb]

        # Convert A to numpy array
        if issparse(A):
            A = np.array(A.todense())
        else:
            A = np.array(A)

        # Trim longer arrays
        f = f[:,:t]
        C = C[:,:t]

        noisyC = cnm_struct.noisyC[:, t - t // 1:t]
        YrA = noisyC[cnm_struct.gnb:cnm_struct.M,:t] - C


        if 'csc_matrix' not in str(type(A)):
            A = csc_matrix(A)
        if 'array' not in str(type(b)):
            b = b.toarray()

        nr, T = C.shape
        nb = f.shape[0]

        Y_r = YrA + C
        d1, d2 = cnm_struct.dims

        all_ = list(range(0, cnm_struct.N))
        try:
            accepted = cnm_struct.accepted
            rejected = [item for item in all_ if item not in accepted]
        except:
            accepted = all_
            rejected = []

        self.axcomp = self.fig.add_axes([0.05, 0.05, 0.9, 0.03])
        self.ax1 = self.fig.add_axes([0.05, 0.55, 0.4, 0.4])
        self.ax3 = self.fig.add_axes([0.55, 0.55, 0.4, 0.4])
        self.ax2 = self.fig.add_axes([0.05, 0.1, 0.9, 0.4])
        self.s_comp = Slider(self.axcomp, 'Component', 0, nr + nb - 1, valinit=0)

        vmax = np.percentile(img, 95)

        def update(val):
            i = np.int(np.round(self.s_comp.val))

            if i < nr:

                if i < len(accepted):
                    j = accepted[i]
                    rej = False
                else:
                    j = rejected[i-len(accepted)]
                    rej = True

                self.ax1.cla()
                imgtmp = np.reshape(A[:, j].toarray(), (d1, d2), order='F')
                self.ax1.imshow(imgtmp, interpolation='None', cmap=pl.cm.gray, vmax=np.max(imgtmp)*0.5)
                if not rej:
                    self.ax1.set_title('Spatial component ' + str(i+1))
                else:
                    self.ax1.set_title('Rejected spatial component ' + str(i-len(accepted)+1))
                self.ax1.axis('off')


                self.ax2.cla()
                self.ax2.plot(np.arange(T), Y_r[j], 'c', linewidth=2)
                self.ax2.plot(np.arange(T), C[j], 'r', linewidth=2)
                if not rej:
                    self.ax2.set_title('Temporal component ' + str(i+1))
                else:
                    self.ax2.set_title('Rejected temporal component ' + str(i-len(accepted)+1))

                self.ax2.legend(labels=['Filtered raw data', 'Inferred trace'])


                self.ax3.cla()
                self.ax3.imshow(img, interpolation='None', cmap=pl.cm.gray, vmax=vmax)
                imgtmp2 = imgtmp.copy()
                imgtmp2[imgtmp2 == 0] = np.nan
                self.ax3.imshow(imgtmp2, interpolation='None', alpha=0.5, cmap=pl.cm.hot)
                self.ax3.axis('off')

            else:
                self.ax1.cla()
                bkgrnd = np.reshape(b[:, i - nr], (d1, d2), order='F')
                self.ax1.imshow(bkgrnd, interpolation='None')
                self.ax1.set_title('Spatial background ' + str(i + 1 - nr))
                self.ax1.axis('off')

                self.ax2.cla()
                self.ax2.plot(np.arange(T), np.squeeze(np.array(f[i - nr, :])))
                self.ax2.set_title('Temporal background ' + str(i + 1 - nr))

        self.nr = nr
        self.nb = nb

        self.s_comp.on_changed(update)
        self.s_comp.set_val(0)
        self.figCanvas.draw()
        self.figCanvas.flush_events()
#        self.updateStatusBar('OnACID traces plotted in plot tab')


    def arrow_key_image_control(self, event):
        print(event.key())
        if event.key() == Qt.Key_Left:
            new_val = np.round(self.s_comp.val - 1)
            if new_val < 0:
                new_val = 0
            self.s_comp.set_val(new_val)

        elif event.key() == Qt.Key_Right:
            new_val = np.round(self.s_comp.val + 1)
            if new_val > self.nr + self.nb:
                new_val = self.nr + self.nb
            self.s_comp.set_val(new_val)
        else:
            pass


    def loadRefMoviePath(self):
        if self.UseOpsinMask_radioButton.isChecked() or self.UseOnacidMask_radioButton.isChecked():
            ref_movie_path = str(QFileDialog.getOpenFileName(self, 'Load reference', self.movie_folder, 'MPTIFF (*.tif);;All files (*.*)')[0])
        else:
            ref_movie_path = str(QFileDialog.getOpenFileName(self, 'Load reference', self.movie_folder, 'PKL(*.pkl);;MPTIFF (*.tif);;All files (*.*)')[0])

        if ref_movie_path and ref_movie_path != 'U:/simulate_movie/20170823.tif':
            self.ref_movie_path = ref_movie_path
            self.movie_folder = os.path.dirname(os.path.dirname(self.ref_movie_path))

            self.refMoviePath_lineEdit.setText(self.ref_movie_path)
            print('Reference movie selected')

    def resetMask(self):
        self.UseOnacidMask_radioButton.setChecked(False)
#        self.UseOpsinMask_radioButton.setChecked(False)
        self.A_loaded = np.array([])
        self.updateImage(np.zeros([512,512]))

    def loadMask(self):
        self.updateStatusBar('')
        
        mask_path = str(QFileDialog.getOpenFileName(self, 'Load reference', '', 'PKL (*.pkl);;NPZ (*.npz)')[0])
        self.mask_path = mask_path
        self.movie_folder = os.path.dirname(os.path.dirname(self.mask_path))
        
        if self.mask_path:
            loaded_mask_ext = os.path.splitext(self.mask_path)[1]

        if loaded_mask_ext == '.npz' or loaded_mask_ext == '.pkl':
            self.updateStatusBar('Loading the mask...')

        if 'init' in self.mask_path:
            init_file = True
        else:
            init_file = False

        if loaded_mask_ext == '.npz': # faster than pkl for loading mask
            # TODO: problem when loading cnmf npz files - issparse = False but matrix sparse
            data_loaded = np.load(self.mask_path)
            try:
                A_loaded = data_loaded['cnm_A']
            except: 
                try:
                    A_loaded = data_loaded['A']
                except:
                    print('No mask found in the npz file')
                    self.resetMask()
                    return
            dims = data_loaded['dims']

        elif loaded_mask_ext == '.pkl':
            data_loaded = load_object(self.mask_path)
#            if 'cnm_init' in list(data_loaded.keys()):
#                init_file = True
#            else: 
#                init_file = False
            try:
                if init_file:
                    cnm = data_loaded['cnm_init']
                    A_loaded = cnm.A
                else:
                    cnm = data_loaded['cnm2']
                    A_loaded = cnm.Ab[:, cnm.gnb:]
            except:
                print('No mask found in the pkl file')
                self.resetMask()
                return
            dims = cnm.dims


        if issparse(A_loaded):
            A_loaded = np.array(A_loaded.todense())
            
        self.A_loaded = np.array(A_loaded)


        if self.A_loaded.size:
            ds_factor = self.dsFactor_doubleSpinBox.value()
            expected_dim = int(512/ds_factor)
            equal_size = list(dims) == [expected_dim for i in range(2)]
            if not equal_size:
                self.resetMask()
                self.updateStatusBar('Mask dimension mismatch. Change downsampling to ' + str(512/dims[0])[:4] +
                                     ' or select a different mask')
                return
            
            try:
                self.opsinMaskOn = False
                self.dims = dims
                loaded_mask = np.reshape(np.array(self.A_loaded.mean(axis=1)), dims, order='F')
                self.updateImage(cv2.resize(np.squeeze(loaded_mask), (512, 512), interpolation=cv2.INTER_CUBIC))
                self.updateStatusBar('Mask with ' + str(self.A_loaded.shape[-1]) + ' cells was loaded')
                self.UseOnacidMask_radioButton.setEnabled(True)
                self.UseOnacidMask_radioButton.setChecked(True)
            except Exception as e:
                print(e)
        else:
            self.updateStatusBar('No mask was not found in the file selected')
            self.UseOnacidMask_radioButton.setChecked(False)


    def initialiseCaiman(self):
        self.deleteTextItems()

        if self.ref_movie_path:
            movie_ext = os.path.splitext(p['refMoviePath'])[1]
            self.updateStatusBar('Initialising...')
        else:
            self.updateStatusBar('No reference movie selected!')
            return

        
        # init parameters
        if movie_ext == '.tif':    
            ds_factor = self.dsFactor_doubleSpinBox.value()
            K = self.InitNumROIs_spinBox.value()
        
            cell_radius = self.cellRadius_spinBox_2.value()
            gSig = (cell_radius, cell_radius)
            min_SNR = self.minSNR_doubleSpinBox.value()
            rval_thr = self.rval_thresh_doubleSpinBox.value()  # for component quality
            merge_thresh = self.merge_thresh_doubleSpinBox.value()
            thresh_overlap = self.overlap_thresh_doubleSpinBox.value()


        # init method to be used
        opsin_seeded = self.UseOpsinMask_radioButton.isChecked()
        mask_seeded = self.UseOnacidMask_radioButton.isChecked()
        cnmf_init = not (opsin_seeded or mask_seeded)

        # process image or load from pkl file directly
        if cnmf_init:
            if movie_ext == '.tif':
                K = self.InitNumROIs_spinBox.value()
                print('Starting CNMF initialisation with tiff file')
                lframe, init_values = initialise(self.ref_movie_path, init_method='cnmf', 
                                                 K=K, gSig=gSig, initbatch=500, 
                                                 ds_factor=ds_factor, rval_thr=rval_thr,
                                                 thresh_overlap=thresh_overlap,
                                                 save_init=True, mot_corr=True,
                                                 merge_thresh=merge_thresh,
                                                 decay_time=0.2, min_SNR=min_SNR)

            elif movie_ext == '.pkl':
                print('Loading initialisation file')
                init_values = load_object(self.ref_movie_path)
                ds_factor = init_values['ds_factor']
                self.dsFactor_doubleSpinBox.setValue(ds_factor)


        elif opsin_seeded:
            if self.A_opsin.size:
                if movie_ext == '.tif':
                    expected_dim = int(512/ds_factor)
                    mask_dim = int(np.sqrt(self.A_opsin.shape[0]))
                    equal_size = mask_dim == expected_dim
                    if not equal_size:
                        print('Change of donwsampling detected; resizing opsin mask')
                        self.createOpsinMask()

                    print('Starting seeded initialisation using the opsin mask')
                    lframe, init_values = initialise(self.ref_movie_path, init_method='seeded',
                                                     Ain=self.A_opsin, gSig=gSig,
                                                     ds_factor=ds_factor, initbatch=500,
                                                     rval_thr=rval_thr, merge_thresh=merge_thresh,
                                                     thresh_overlap=thresh_overlap,
                                                     save_init=True, mot_corr=True,
                                                     min_SNR=min_SNR, CNN_filter=True)

                else:
                    self.updateStatusBar('A non-tif file was provided - select a tif movie to initialise with a mask')
                    return
            else:
                self.updateStatusBar('No opsin mask created!')


        elif mask_seeded:
            if movie_ext == '.tif':
                expected_dim = int(512/ds_factor)
                mask_dim = int(np.sqrt(self.A_loaded.shape[0]))
                equal_size = mask_dim == expected_dim
                if not equal_size:
                    self.updateStatusBar('Mask dimension mismatch. Change downsampling back to ' + str(512/mask_dim)[:4] +
                                         ' or select a different mask')
                    return

                print('Starting seeded initialisation using the mask provided')
                lframe, init_values = initialise(self.ref_movie_path, init_method='seeded', Ain=self.A_loaded,
                                                 ds_factor=ds_factor, initbatch=500, rval_thr=0.85,
                                                 thresh_overlap=0.1, save_init=True, mot_corr=True,
                                                 merge_thresh=0.85)

            else:
                self.updateStatusBar('A non-tif file was provided - select a tif movie to initialise with a mask')
                return

        # Prepare object for OnACID
        idx_components = init_values['idx_components']
        path_to_cnn_residual = os.path.join(caiman_datadir(), 'model', 'cnn_model_online.h5')

        cnm2 = deepcopy(init_values['cnm_init'])
        cnm2._prepare_object(np.asarray(init_values['Yr']), init_values['T1'],
                             init_values['expected_comps'], idx_components=idx_components,
                             min_num_trial=2, N_samples_exceptionality=int(init_values['N_samples']),
                             path_to_model=path_to_cnn_residual)

        cnm2.opsin = None
        cnm2.accepted = list(range(0,cnm2.N))
        cnm2.t = cnm2.initbatch

        # Extract number of cells detected
        K = cnm2.N
        self.InitNumROIs_spinBox.setValue(K)
        print('Number of components initialised: ' + str(K))

        if self.MaxNumROIs_spinBox.value() < K+5:
            self.MaxNumROIs_spinBox.setValue(K+10)
            
        if movie_ext == '.pkl':
            cnm = init_values['cnm_init']
            cell_radius = cnm.gSig[0]*ds_factor
            rval_thr = cnm.rval_thr
            merge_thresh = cnm.merge_thresh
            thresh_overlap = cnm.thresh_overlap
            try:
                min_SNR = init_values['min_SNR']
            except:
                min_SNR = 2.5 # default    
            
            self.cellRadius_spinBox_2.setValue(cell_radius)
            self.rval_thresh_doubleSpinBox.setValue(rval_thr)
            self.merge_thresh_doubleSpinBox.setValue(merge_thresh)
            self.overlap_thresh_doubleSpinBox.setValue(thresh_overlap)
            self.minSNR_doubleSpinBox.setValue(min_SNR)
                
        
        self.cell_radius = cell_radius
#        self.gSig = (cell_radius, cell_radius)
        self.min_SNR = min_SNR
        self.rval_thr = rval_thr
        self.merge_thresh = merge_thresh
        self.thresh_overlap = thresh_overlap
        self.ds_factor = ds_factor
        self.K_init = K

        self.c = init_values
        self.c['cnm2'] = cnm2
        self.c['removed_idx'] = []

        if idx_components is not None:
            coms = self.c['coms_init']
            self.c['coms_init'] = coms[idx_components]


        self.dims = init_values['cnm_init'].dims
        self.InitNumROIs = K
        self.opsinMaskOn = False
        self.imageItem.setImage(cv2.resize(self.c['img_norm'],(512,512),interpolation=cv2.INTER_CUBIC)) # display FOV

        if self.A_opsin.size:
            print('Checking opsin overlap')
            self.checkOpsin()

        self.initialiseROI()
        for i in range(K):
            y, x = init_values['coms_init'][i]  # reversed
            self.getROImask(thisROIx = x, thisROIy = y)

        self.UseONACID_checkBox.setEnabled(True)
        self.UseONACID_checkBox.setChecked(True)
        if self.selectAll_checkBox.isChecked():
            self.selectAllROIs()

        self.updateStatusBar('Initialision completed')
        show_traces = 1
        if show_traces:
            self.showROIIdx()
            cnm_init = init_values['cnm_init']
            self.plotOnacidTraces(t=cnm_init.initbatch)


    def updateInitialisation(self):
        keep_idx = list(set(range(self.thisROIIdx)) - set(self.removeIdx))
        
        # reset previous settings
        self.ROIcontour_item.clear()
        self.thisROIIdx = 0
        self.resetROIlist()
        self.deleteTextItems()
        self.updateTable()
        self.plotItem.clear()
        self.resetFigure()

        # remove chosen cells
        self.c['cnm2'].remove_components(self.removeIdx)

        K = self.c['cnm2'].N
        print('new K:', K)
        self.K_init = K

        self.c['cnm2'].accepted = self.c['cnm2'].accepted[:K] # all cells accepted post initialisation

        coms = self.c['coms_init']
        self.c['coms_init'] = coms[keep_idx]

        self.InitNumROIs_spinBox.setValue(K)
        self.InitNumROIs = K

        if self.MaxNumROIs_spinBox.value() < K+5:
            self.MaxNumROIs_spinBox.setValue(K+10)

        self.initialiseROI()
        for i in range(K):
            y, x = self.c['coms_init'][i]  # reversed
            self.getROImask(thisROIx = x, thisROIy = y)

        show_traces = 1
        if show_traces:
            self.showROIIdx()
            cnm_init = self.c['cnm_init']
            self.plotOnacidTraces(t=cnm_init.initbatch)

        # after new roi list initialised
#            print('target idx list', self.TargetIdx)
        for idx in sorted(self.removeIdx, reverse=True):
            if idx in self.TargetIdx:
#                    print('idx', idx)
                self.TargetIdx.remove(idx)
            self.TargetIdx = [target-1 if target>idx else target for target in self.TargetIdx]

        if self.A_opsin.size:
            self.checkOpsin()

        self.updateTargets()
        self.updateStatusBar('Removed cells: ' + str([value+1 for value in self.removeIdx])[1:-1])

        self.removeIdx = []


    def loadOpsinImg(self):
        opsin_img_path = str(QFileDialog.getOpenFileName(self, 'Load a C1V1 image', self.movie_folder, 'MPTIFF (*.tif)')[0])

        self.opsin_img_path = opsin_img_path

        if self.opsin_img_path:
            self.opsinImgPath_lineEdit.setText(self.opsin_img_path)
            opsin_img = cm.load(opsin_img_path, subindices = slice(0,1,None))
            self.updateImage(opsin_img)
            self.opsin_img = opsin_img


    def createOpsinMask(self):
        if self.opsin_img.size:
            ds_factor = self.dsFactor_doubleSpinBox.value()

            opsin_img = self.opsin_img[np.newaxis,:,:].resize(1./ds_factor,1./ds_factor) # moved here for init recreating

            min_area = self.minArea_spinBox.value()
            radius = self.cellRadius_spinBox.value() # must be radius % 2 == 1 and > 1

            if radius % 2 == 0:
                radius += 1
            self.cellRadius_spinBox.setValue(radius)

            dims = opsin_img.shape[-2:]

            A_opsin, mr_opsin = extract_mask(opsin_img, gSig=radius, min_area_size=min_area, min_hole_size=15)
            self.A_opsin = A_opsin.astype('int')
            self.opsin_mask = np.reshape(np.array(self.A_opsin.max(axis=1)), dims, order='F').astype('int')  # max because it's a binary mask
            self.updateImage(cv2.resize(np.squeeze(self.opsin_mask).astype('u1'), (512, 512), interpolation=cv2.INTER_CUBIC))

            self.opsinMaskOn = True
            self.UseOpsinMask_radioButton.setEnabled(True)
            self.updateStatusBar('Opsin mask created')

            # screen current components for c1v1
            if self.c != {}:
                self.checkOpsin()

        else:
            self.updateStatusBar('No opsin image loaded!')
            

    def checkOpsin(self, Ain=None):
        if Ain is None:
            cnm_struct = self.c['cnm2']
            A = cnm_struct.Ab[:, cnm_struct.gnb:cnm_struct.M]
#            print('a shape', A.shape)
            accepted = self.c['cnm2'].accepted

        else:
            A = Ain
            accepted = range(list(0, A.shape[-1]))

        dims = self.dims # tuple([np.sqrt(A.shape[0])]*2)
        overlap = []
        opsin = []
        self.opsin_thresh = self.overlapRatio_doubleSpinBox.value()

        # Convert A to numpy array
        if issparse(A):
            A = np.array(A.todense())
        else:
            A = np.array(A)

        onacid_mask = (deepcopy(A)>0).astype('int')  # binarise the onacid output mask

        for cell in accepted: # range(onacid_mask.shape[-1]):
            cell_mask = (np.reshape(onacid_mask[:,cell], dims, order='F'))
            cell_pix = sum(sum(cell_mask == 1))

            inter = cv2.bitwise_and(self.opsin_mask, cell_mask)
            inter_pix = sum(sum(inter))
            cell_overlap = inter_pix/cell_pix
            overlap.append(cell_overlap)

            if cell_overlap <= self.opsin_thresh:
                onacid_mask[:,cell][onacid_mask[:,cell] == 1] = -3
            else:
                onacid_mask[:,cell][onacid_mask[:,cell] == 1] = 3

            opsin.append(cell_overlap > self.opsin_thresh)

        if Ain is None:
            # store info on opsin
            self.c['cnm2'].opsin = opsin  # True/False based on threshold
            self.c['overlap'] = overlap   # value
            self.c['opsin_mask'] = self.opsin_mask
            self.c['opsin_thresh'] = self.opsin_thresh
            self.onacid_mask = onacid_mask

        if self.selectAll_checkBox.isChecked():
            self.selectAllROIs()  # update ROIs selected

        return onacid_mask


    def showCellsOnMask(self):
        try:
            if self.A_opsin.size and self.c != {}:
                self.checkOpsin()  # for overlap update

                # visualise all comps
                summed_A = np.hstack((self.A_opsin, self.onacid_mask))
                summed_mask = np.reshape(np.array(summed_A.sum(axis=1)), self.dims, order='F')
                pl.figure();pl.imshow(summed_mask)
                pl.colorbar()


            elif self.A_opsin.size and self.A_loaded.size:
                loaded_mask = self.checkOpsin(self.A_loaded)

                summed_A = np.hstack((self.A_opsin, loaded_mask))
                summed_mask = np.reshape(np.array(summed_A.sum(axis=1)), self.dims, order='F')
                pl.figure();pl.imshow(summed_mask)
                pl.colorbar()

            else:
                self.updateStatusBar('No cells to be shown!')
        except Exception as e:
            print(e)


    def getROImask(self, thisROIx, thisROIy):
        if self.ALL_ROI_SELECTED == False:
            if self.ds_factor > 1:
                thisROIx = round(self.ds_factor*thisROIx)
                thisROIy = round(self.ds_factor*thisROIy)
            else:
                thisROIx = round(thisROIx)
                thisROIy = round(thisROIy)

            qcolor = self.ROIlist[self.thisROIIdx]["ROIcolor"]

            try:
                self.ROIcontour_item.addPoints(x = [thisROIx], y = [thisROIy], pen = pg.mkPen(qcolor, width=2), size = self.RoiRadius*2)

                tempbrush = QBrush(qcolor)
                self.X_labels[self.thisROIIdx].setForeground(tempbrush)
                self.Y_labels[self.thisROIIdx].setForeground(tempbrush)
                self.Thresh_labels[self.thisROIIdx].setForeground(tempbrush)

                self.X_labels[self.thisROIIdx].setText(str('{:.0f}'.format(thisROIx)))
                self.Y_labels[self.thisROIIdx].setText(str('{:.0f}'.format(thisROIy)))

                self.Thresh_tableWidget.setItem(self.thisROIIdx,0,self.X_labels[self.thisROIIdx])
                self.Thresh_tableWidget.setItem(self.thisROIIdx,1,self.Y_labels[self.thisROIIdx])
                self.Thresh_tableWidget.setItem(self.thisROIIdx,2,self.Thresh_labels[self.thisROIIdx]) # threshold table
            except Exception as e:
                print(e)

            # save to list
            self.ROIlist[self.thisROIIdx]["ROIx"] = thisROIx
            self.ROIlist[self.thisROIIdx]["ROIy"] = thisROIy
            self.thisROIIdx += 1
            if self.thisROIIdx == self.MaxNumROIs:
                self.ALL_ROI_SELECTED = True

        else:
            print('All ROI selected')


    def addScatterROI(self, shift):
        self.ROIcontour_item.clear()
        self.ROIcontour_item.addPoints(x = [ROI["ROIx"]+shift[0] for ROI in self.ROIlist[:self.thisROIIdx]],
                                           y = [ROI["ROIy"]+shift[1] for ROI in self.ROIlist[:self.thisROIIdx]],
                                           pen = [pg.mkPen(ROI["ROIcolor"], width=2) for ROI in self.ROIlist[:self.thisROIIdx]],
                                           size = self.RoiRadius*2)

    def showROIIdx(self):
        print('show num rois = '+str(self.thisROIIdx))
        font = QFont("Arial", 12, QFont.Bold)
        for i in range(self.thisROIIdx):
            thisText = pg.TextItem(text=str(i+1), color=self.ROIlist[i]["ROIcolor"].getRgb()[:3],
                                   html=None, anchor=(0,0), border=None, fill=None, angle=0, rotateAxis=None)
            thisText.setPos(self.ROIlist[i]["ROIx"],self.ROIlist[i]["ROIy"] )
            thisText.setFont(font)
            self.graphicsView.addItem(thisText)

    def updateROI(self,img):
        # import pdb
        # pdb.set_trace()
        self.updateImage(img)

    def selectAllROIsClicked(self):
        if (self.selectAll_checkBox.isChecked):
            self.selectAllROIs()

    def autoAddClicked(self):
        p['FLAG_AUTO_ADD_TARGETS'] = self.addNewROIsToTarget_checkBox.isChecked()
        print('FLAG_AUTO_ADD_TARGETS = '+str(p['FLAG_AUTO_ADD_TARGETS']))

    def enterEvent(self,event):
        self.graphicsView.setCursor(Qt.CrossCursor)

    def takeRefMovie(self):
        # take reference movie by starting t-series in prairie and save as multipage tiffs
        global STOP_JOB_FLAG
        STOP_JOB_FLAG = False

        # clear items in qbuffer
        if not qbuffer.empty():
            qbuffer.clear()

        # disable onacid
        self.UseONACID_checkBox.setCheckState(Qt.Unchecked)
        self.getValues()

        if self.FLAG_PV_CONNECTED:
            print('this is take ref movie function')

            # setup thread objects
            save_path = self.pl.get_movie_name()+'_rtaoi.tif'
            print(save_path)

            # set this file as the ref movie
            self.ref_movie_path = save_path
            self.refMoviePath_lineEdit.setText(self.ref_movie_path)


            kwargs = {"save_movie_path": save_path}
            self.imageSaverObject = imageSaver(**kwargs)
            self.updateStatusBar('imageSaver object created')
            self.imageSaverObject.moveToThread(self.imageSaverThread)
            self.imageSaverThread.started.connect(self.imageSaverObject.saveImage)
            self.imageSaverObject.finished_signal.connect(self.imageSaverThread.exit)
            self.imageSaverThread.start()

            # start stream thread
            kwargs = {"prairie": self.pl}
            self.streamObject = DataStream(**kwargs)
            self.streamObject.moveToThread(self.streamThread)
            self.streamThread.started.connect(self.streamObject.stream)
            self.streamObject.stream_finished_signal.connect(self.streamThread.exit)
            self.streamThread.start()
            self.updateStatusBar('stream started in takeRefMovie function')
        else:
            print('check pv connection')
            self.updateStatusBar('Check PV connection')


    def clickRun(self):
        global STOP_JOB_FLAG
        STOP_JOB_FLAG = False
        
        if not qbuffer.empty():
#                qbuffer.clear()  # error: 'Queue' object has no attribute 'clear'
            with qbuffer.mutex:
                qbuffer.queue.clear()
                print('qbuffer cleared!')
        
        # connect to pv if not offline
        if (not self.IsOffline) and (not self.FLAG_PV_CONNECTED):
            self.connectPV

        # get current movie name if PV is connected
        if self.FLAG_PV_CONNECTED:
            p['currentMoviePath'] = self.pl.get_movie_name()+ '_DS_' + str(self.ds_factor) +'_rtaoi.tif'
            print('current movie:' + p['currentMoviePath'])
        else:
            p['currentMoviePath'] = os.path.splitext(p['moviePath'])[0]+ '_DS_' + str(self.ds_factor) +'_rtaoi.tif'


        if self.IsOffline or self.FLAG_PV_CONNECTED: # p['UseONACID'] or self.FLAG_PV_CONNECTED
            try: self.resetFigure()
            except: pass
        
            self.c['keep_prev'] = False # default
    
            if self.c['cnm2'].t > self.c['cnm2'].initbatch: # self.cnm2:
                msg = QMessageBox()
                msg.setText("Run again?")
                msg.setWindowTitle('pyRTAOI Message')
                msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
                retval = msg.exec_()
                if(retval==QMessageBox.Ok):
                    pass
                    msg.setText("Keep the cells from the previous session?")
                    msg.setWindowTitle('pyRTAOI Message')
                    msg.setStandardButtons(QMessageBox.Yes | QMessageBox.No)
                    retval = msg.exec_()
                    if(retval==QMessageBox.No):
                        curr_count = self.c['cnm2'].N
                        self.c['cnm2'].t = self.c['cnm2'].initbatch
                        self.removeIdx = np.arange(self.K_init,curr_count)
                        self.updateInitialisation()
                    else:
                        self.c['keep_prev'] = True
                else:
                    return
    
            self.deleteTextItems()

            if self.opsinMaskOn:
                self.opsinMaskOn = False
                self.updateImage(cv2.resize(self.c['img_norm'],(512,512),interpolation=cv2.INTER_CUBIC))

            # ensure correct param values following init are displayed
            self.dsFactor_doubleSpinBox.setValue(self.ds_factor)
            self.cellRadius_spinBox_2.setValue(self.cell_radius)
            self.minSNR_doubleSpinBox.setValue(self.min_SNR)
            self.rval_thresh_doubleSpinBox.setValue(self.rval_thr)
            self.merge_thresh_doubleSpinBox.setValue(self.merge_thresh)
            self.overlap_thresh_doubleSpinBox.setValue(self.thresh_overlap)

            self.getValues()

            # connect stop_thread to stop button when analysis running
            self.stop_pushButton.clicked.connect(self.stop_thread)

            def resetworkerObject():
                self.workerObject = None

            def resetstreamObject():
                self.streamObject = None
                

            kwargs = {"caiman": self.c,"prairie": self.pl,"ROIlist":self.ROIlist}
            self.workerObject = Worker(**kwargs)
            self.updateStatusBar('Worker created')
            self.workerObject.moveToThread(self.workerThread)

            # update display signals
            self.workerObject.status_signal.connect(self.updateStatusBar)
            self.workerObject.roi_signal.connect(self.updateROI)
            self.workerObject.refreshPlot_signal.connect(self.refreshPlot)
            self.workerObject.frame_signal.connect(self.updateFrameInfo)
            self.workerObject.thresh_signal.connect(self.updateThresh)
            self.workerObject.sta_amp_signal.connect(self.updateSTA)
            self.workerObject.refreshScatter_signal.connect(self.addScatterROI)
            self.workerObject.showROIIdx_signal.connect(self.showROIIdx)
            self.workerObject.getROImask_signal.connect(self.getROImask)
            self.workerObject.updateTargetROIs_signal.connect(self.updateTargetROIs)

            # start and finish
            self.workerObject.transDataToMain_signal.connect(self.transDataToMain)
            self.workerThread.started.connect(self.workerObject.work)
            self.workerObject.finished_signal.connect(self.workerThread.exit)

            # reset worker thread after finishing and disconnect stop button
            self.stopButtonDisconnect = lambda: self.stop_pushButton.clicked.disconnect(self.stop_thread)
            self.workerObject.finished_signal.connect(self.stopButtonDisconnect)
            self.workerObject.finished_signal.connect(self.workerObject.deleteLater)  # delete qt (C++) instance later
            self.workerObject.finished_signal.connect(resetworkerObject) # remove python reference to workerObject

            # triggers
            self.workerObject.sendTTLTrigger_signal.connect(self.sendTTLTrigger)
            self.workerObject.sendCoords_signal.connect(self.sendCoords)
            self.workerObject.sendPhotoStimTrig_signal.connect(self.sendPhotoStimTrigger)
            self.workerThread.start()
            self.updateStatusBar('Worker started')

            # start stream thread
            kwargs = {"prairie": self.pl}
#            if not self.streamObject:  # could use this if self.pl not updated/separate update function existed
            self.streamObject = DataStream(**kwargs)
            self.streamObject.moveToThread(self.streamThread)
            self.streamThread.started.connect(self.streamObject.stream)
            self.streamObject.stream_finished_signal.connect(self.streamThread.exit)
            self.streamObject.stream_finished_signal.connect(self.streamObject.deleteLater)  # seems it doesn't get deleted but stop disconnect sorts it
#            self.streamObject.stream_finished_signal.connect(resetstreamObject)  # looks like stream object not deleted though
            self.streamThread.start()
            self.updateStatusBar('stream started')

        else:
            self.updateStatusBar('No initialisation provided and PV not connected')
            return


    def refreshPlot(self, arr):
        self.plotItem.clear()
        for i in range(self.thisROIIdx):
            self.plotItem.plot(arr[i,:], antialias=True, pen=pg.mkPen(self.ROIlist[i]["ROIcolor"], width=1))
        try:
            self.plotItem.setYRange(np.round(arr.min())-1,np.round(arr.max())+1)
        except:
            self.plotItem.setYRange(0,30)


    def updateImage(self,img):
        self.imageItem.setImage(img)

    def updateStatusBar(self, msg):
        self.statusbar.showMessage(msg)

    def updateThresh(self,thresh_arr):
        try:
            for i in range(self.thisROIIdx):
                self.Thresh_labels[i].setText(str('{0:.1f}'.format(thresh_arr[i])))
        except Exception as e:
            print(e)

    def updateSTA(self, sta_arr):
        for i in range(0,len(sta_arr)):
#            self.STA_labels[i].setForeground(QBrush(self.ROIlist[i]["ROIcolor"]))
            self.STA_labels[i].setText(str('{0:.1f}'.format(sta_arr[i])))
            self.Thresh_tableWidget.setItem(i,3,self.STA_labels[i])
        print('updated sta value')

    def updateFrameInfo(self,FrameIdx):
        self.CurrentFrameIdx_label.setText(str(FrameIdx))

    def switch_plotOn(self):
        p['plotOn'] = self.plotOn_checkBox.isChecked()
        print(p['plotOn'])

    def switch_displayOn(self):
        p['displayOn'] = self.displayOn_checkBox.isChecked()
        print(p['displayOn'])

    def switch_denoiseOn(self):
        p['denoiseOn'] = self.denoiseOn_checkBox.isChecked()

    def switch_IsOffline(self):
        print('offline button')
        self.IsOffline = self.IsOffline_radioButton.isChecked()
        self.updateStatusBar('Flag offline = '+str(self.IsOffline))


    def switch_useOnacidMask(self):  # toggle the two buttons
        self.usesOnacidMask = self.UseOnacidMask_radioButton.isChecked()
        self.usesOpsinMask = self.UseOpsinMask_radioButton.isChecked()
        if self.usesOpsinMask and self.usesOnacidMask:
            self.UseOpsinMask_radioButton.setChecked(False)

    def switch_useOpsinMask(self):
        self.usesOpsinMask = self.UseOpsinMask_radioButton.isChecked()
        self.usesOnacidMask = self.UseOnacidMask_radioButton.isChecked()
        if self.usesOpsinMask and self.usesOnacidMask:
            self.UseOnacidMask_radioButton.setChecked(False)


    def enableStimTrigger(self):
        p['FLAG_STIM_TRIG_ENABLED'] = self.enableStimTrigger_checkBox.isChecked

    def setupPowerControl(self):
        self.getValues()
        NI_TTL_NUM_SAMPLES = p['NI_TTL_NUM_SAMPLES']
        NI_STIM_NUM_SAMPLES = p['NI_STIM_NUM_SAMPLES']
        p['SendPowerVolt'] = self.SendPowerVolt_checkBox.isChecked
        init_output = np.zeros([2, max(NI_TTL_NUM_SAMPLES+1,NI_STIM_NUM_SAMPLES+1)])
        init_output.ravel()[:len(self.daq_array)] = self.daq_array  # use ravel to data between array (fastest way)
        init_output[1,:-1] = 1
        p['NI_2D_ARRAY'] = init_output
        p['NI_UNIT_POWER_ARRAY'] = init_output[1,:]
        self.niPhotostimFullWriter.write_many_sample(init_output)

    def testPowerControl(self):
        self.getValues()
        p['NI_2D_ARRAY'][1,:] = p['NI_UNIT_POWER_ARRAY'] *np.polyval(p['power_polyfit_p'],p['photoPowerPerCell'])
        self.niPhotostimFullWriter.write_many_sample(p['NI_2D_ARRAY'])
        print('TTL with power volt sent')

    def stop_thread(self):  # TODO: aborting t-series correct? when aborted, no return from prairie so disabled waiting

        p['FLAG_END_LOADING'] = True
        
#        self.stop_pushButton.clicked.disconnect(self.stop_thread)
#        print('disconnected')

        # this worked a few times after which the timeout error occurred again
#        if p['FLAG_PV_CONNECTED']:
#            if not self.pl.send_done(self.pl._abort):
#                print('Prairie aborted!')

        try:
            if self.streamObject:
                self.streamObject.stop()
                self.streamThread.wait(1)
                self.streamThread.quit()
                self.streamThread.wait()
        except: pass
    
        try:
            if self.workerObject:
                self.workerObject.stop()
                self.workerThread.wait(1)
                self.workerThread.quit()
                self.workerThread.wait() # ensures worker finished before next run can be started
        except: pass


        self.stop_pushButton.clicked.disconnect(self.stop_thread)
        print('stop button disconnected')


    def closeEvent(self,event):
        # override closeEvent method in Qt
        print('form closing, PV connection = ' + str(self.FLAG_PV_CONNECTED))

        if self.FLAG_PV_CONNECTED:
            del self.pl
            print('pv disconnected')
        if p['FLAG_BLINK_CONNECTED']:
            self.bl.abort()
            del self.bl
            print('blink disconnected')

        if self.workerObject:
            try:
                print('stopping worker')
                self.workerObject.stop()
                self.workerThread.wait(1)
                self.workerThread.quit()
                self.workerObject.deleteLater()
                print('worker exit')
            except Exception as e:
                print(e)

        if self.streamObject:
            try:
                print('stopping stream')
                self.streamObject.stop()
                self.streamThread.wait(1)
                self.streamThread.quit()
                self.streamObject.deleteLater()
                print('stream exit')
            except Exception as e:
                print(e)

        # reset AO to zero, clear ni tasks
        try:
            self.niStimWriter.write_one_sample(0,10.0)
            self.niStimTask.stop()
            self.niStimTask.close()
            del self.niStimTask
        except Exception as e:
            print(e)

        try:
#                        # zero outputs
#            zero_output =  np.zeros(p['NI_2D_ARRAY'].shape)
#            self.niPhotostimFullWriter.write_many_sample(zero_output,10.0)
#            while(not self.niPhotoStimFullTask.is_task_done()):
#                pass
            self.niPhotoStimFullTask.stop()
            self.niPhotoStimFullTask.close()
            del self.niPhotoStimFullTask
        except Exception as e:
            print('photostim full task clean up error')
            print(e)

        try:
            self.niPhotoStimWriter.write_one_sample(0,10.0)
            self.niPhotoStimTask.stop()
            self.niPhotoStimTask.close()
            del self.niPhotoStimTask
        except Exception as e:
            print(e)

        # delete big structs to free memory
        del self.c
#        del self.cnm2
#        del self.proc_cnm
#            p = {}
        
        plt.close('all')

        # delete all locals
        for name in dir():
            if not name.startswith('_'):
                del locals()[name]


        g = gc.collect()
        print('gc exit: ', g)
        event.accept()


#%% save and load configuration
    def saveConfig(self):
        widgets = (QComboBox, QCheckBox, QLineEdit, QSpinBox, QDoubleSpinBox)
        for obj in self.findChildren(widgets):
            name = str(obj.objectName())
            if isinstance(obj, QComboBox):
                self.trial_config[name] = str(obj.currentText())
            if isinstance(obj, QCheckBox):
                self.trial_config[name] = bool(obj.isChecked())
            if isinstance(obj, QLineEdit):
                if 'spinbox' not in name:
                    self.trial_config[name] = str(obj.text())
            if isinstance(obj, QSpinBox):
                self.trial_config[name] = int(obj.value())
            if isinstance(obj, QDoubleSpinBox):
                self.trial_config[name] = float(obj.value())
        filepath = str(QFileDialog.getSaveFileName(self, 'Save as preset...', 'Configs', 'Config file (*.cfg)')[0])
        json.dump(self.trial_config, open(filepath, 'w'), sort_keys=True, indent=4)
        filename = os.path.basename(filepath)
        filename = os.path.splitext(filename)[0]


    def loadConfig(self):
        filepath = str(QFileDialog.getOpenFileName(self, 'Load configuration...', 'Configs', 'Config file (*.cfg)')[0])
        if os.path.exists(filepath):
            self.trial_config = json.load(open(filepath, 'r'))
            widgets = (QComboBox, QCheckBox, QLineEdit, QSpinBox, QDoubleSpinBox)
            for obj in self.findChildren(widgets):
                name = str(obj.objectName())

                try:
                    if isinstance(obj, QComboBox):
                        value = self.trial_config[name]
                        index = obj.findText(value)  # get the corresponding index for specified string in combobox
                        obj.setCurrentIndex(index)  # preselect a combobox value by index
                    if isinstance(obj, QLineEdit):
                        value = self.trial_config[name]
                        if 'spinbox' not in name:
                            obj.setText(value)  # restore lineEditFile
                    if isinstance(obj, QCheckBox):
                        value = self.trial_config[name]
                        if value is not None:
                            obj.setChecked(value)  # restore checkbox
                    if isinstance(obj, QSpinBox):
                        value = self.trial_config[name]
                        obj.setValue(value)  # restore lineEditFile
                    if isinstance(obj, QDoubleSpinBox):
                        value = self.trial_config[name]
                        obj.setValue(value)  # restore lineEditFile
                except:
                    continue

            self.getValues()


#%%
def main(argv):
    # create Qt application
    app = QCoreApplication.instance()
    if app is None:  # stops kernel crashing
        app = QApplication(argv)

    # create main window
    window = MainWindow()

    # show it and bring to front
    window.show()
    window.raise_()

    # start the app
    sys.exit(app.exec_())



if __name__ == '__main__':
    # setup diractories
    config_directory = 'Configs'
    if not os.path.exists(config_directory):
        os.makedirs(config_directory)

    results_directory = 'Results'
    if not os.path.exists(results_directory):
        os.makedirs(results_directory)

    # initialise parameters
    global p
    p = {}
    p['saveAsTiff'] = False # temp here

   # global data buffer
    global qbuffer
    qbuffer = Queue(maxsize = 0)

    # global stop flag
    global STOP_JOB_FLAG
    STOP_JOB_FLAG = False

    # launch program
    try:
        main(sys.argv)
    except Exception as e:
        logger.exception(e)
        print(str(e))
