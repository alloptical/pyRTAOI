''''
 Dependences:
 follow instructions in these links to install all dependences
 caiman  : https://github.com/flatironinstitute/CaImAn 
 opencv  (pip install .whl) : https://www.lfd.uci.edu/~gohlke/pythonlibs/#opencv
 PyDAQmx (pip) : https://github.com/clade/PyDAQmx

SPYDER HAS NO AUTOSAVE FUNCTION! 
CAREFUL WHEN USING 'REPLACE ALL' - IT WILL QUIT, WITHOUT SAVING!!


 TO DO:
0. data stream stopped in the middle of a movie - not polling data from buffer fast enough

    0.1 nested thread - read data and emit signal when the correct number of sample is collected 
        - saved in pyRTAOI-20180620_1; this guaranteed data streaming from PV but it takes more time to process one frame, 35 ms without nested thread; 60 ms with
    0.2 use queue to communicate between worker and streamer
        - saved in pyRTAOI-20180620_2; < 20 ms when running offline; need to test on rig
    0.3 use gpu for motion correction? - caiman motion correction is below 10 ms 
	
	- down sample by 2 kept up with 30Hz image stream
    
1. log data for offline anaysis - photostimulation targets and frame idx; save online analysis result
    offline analysis - convert cnmf results to .mat file - done
 
2. deal with zoom - check
3. mark point setting for changing photostim power; if not fast enough use dummy targets in zero blocker, or feed voltage to pockels directly
    control voltage = (outMax/displayMax)*PV value; where outMax and displayMax are defined in prairie view configuration.xml
    need to map power to voltage - done
    send volt to pockels via ni card
    
    
 
4. Plot STA dF/F - send triggers - check
5. threshold sometimes becomes Inf - check
6. auto initialisation - take reference movie and load
7. photostim protocol 
8. delete target by right clicking  - done by left click
9. some memory leak - every use of the interface increases memory usage by ~3-4% which can't be reset

'''

import sys
import random
import os

# Qt
import GUI
from PyQt5.QtCore import Qt,QObject, pyqtSignal, QThread, QTimer, QRectF, QUrl,QPoint, QRect, QSize,QPointF,QSizeF, QSettings, QCoreApplication
from PyQt5.QtWidgets import (QComboBox, QCheckBox, QLineEdit, QSpinBox, QLabel,
                             QDoubleSpinBox, QFileDialog, QApplication,
                             QDesktopWidget, QMainWindow, QMessageBox, QTableWidgetItem)
from PyQt5.QtGui import QFont, QColor, QIcon, QPalette, QDesktopServices, QImage, QPixmap, QPainter,QPen,QBrush


# utilities
import gc
from skimage.external import tifffile
import numpy as np
import pyqtgraph as pg
import time
import pylab as pl
import cv2
from copy import deepcopy
from scipy.sparse import issparse, spdiags, coo_matrix, csc_matrix
import pickle
import logging
import json

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

# Ensure using PyQt5 backend
matplotlib.use('QT5Agg', force=True)
import matplotlib.pyplot as plt
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
        if self.FLAG_PV_CONNECTED: 
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
                
            while self.FLAG_SCANNING:
                loop += 1
    #                print('This is loop '+ str(loop))
                loop_start_time = time.time()
                num_samples = self.pl.send_recv(_rrd)
    
                if not num_samples:
                    self.FLAG_SCANNING  = False
                
                elif int(num_samples) == samplesPerFrame:
                       
                    print('recv time: '+str("%.4f"%(time.time()-loop_start_time)))
                    recvTime += time.time()-loop_start_time
                    recv_time = time.time()
                    
                    # -----  use cuda -------
                    buf_gpu = gpuarray.to_gpu(buf.astype(np.int16))
                    sample_mean( buf_gpu, pixelsPerLine, linesPerFrame, samplesPerPixel, p['flipEvenRows'], dest_gpu,
                        block=(1024,1,1), grid = (int(samplesPerFrame/1024),1,1) )
                    dest_gpu = dest_gpu.reshape(linesPerFrame,pixelsPerLine)
                    thisFrame = dest_gpu.get()
                    print('convert time:' + str("%.4f"%(time.time()-recv_time)))

                    # put frame in global queue
                    qcopy_time = time.time()
                    qbuffer.put(thisFrame.copy().astype(np.float32))
                    print('frame to queue:' + str("%.4f"%(time.time()-qcopy_time)))
                    self.framesRecv += 1
                    print('frames recv ='+str(self.framesRecv))
                    incorrectCount = 0
                else:
                    incorrectCount +=1
                    
                if incorrectCount>incorrectLimit: 
                    print('too many incorrect frames, quit loop')
                    self.FLAG_SCANNING   = False
    
            if not self.FLAG_SCANNING:   # - need a flag outside of this thread to stop it
                self.context.pop()
                del self.context
                from pycuda.tools import clear_context_caches
                clear_context_caches()
                del buf_gpu
                del dest_gpu
                del sample_mean
                self.stream_finished_signal.emit()
                print('stopped scanning')
                p['FLAG_END_LOADING'] = True
#%% offline, read frames from tiff                
        else:
            # load movie
            print('Loading video')
            Y_ = cm.load(p['moviePath'], subindices=slice(0,2500,None))
            print('Video loaded')

            for frame_count, frame in enumerate(Y_):        # now process each file                
                if np.isnan(np.sum(frame)):
                    raise Exception('Frame ' + str(frame_count) + ' contains nan')                  
                qbuffer.put(frame.copy().astype(np.float32))
                self.framesRecv += 1
                
            p['FLAG_END_LOADING'] = True   
                
    def stop(self):
        self.FLAG_SCANNING = False        
            
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
    transDataToMain_signal   = pyqtSignal(object,object,name = 'transDataToMain')
    updateTargetROIs_signal  = pyqtSignal()
    finished_signal          = pyqtSignal()
    sendTraces_signal = pyqtSignal(object,object,object,name = 'sendTraces')



    def __init__(self, **kwargs ):
        super(Worker,self).__init__()
        
        # get parameters from Main thread
        self.c = kwargs.get('caiman',{})
        self.pl = kwargs.get('prairie',[])
        self.ROIlist = kwargs.get('ROIlist',dict())
        print('worker pl empty? ' + str(self.pl == []))
        
        # sta traces
        self.sta_traces = np.empty([p['MaxNumROIs'],p['numberStims'],p['staPreFrame']+p['staPostFrame']])
        self.sta_traces.fill(np.nan)
        p['stimStopFrame'] = p['stimStartFrame']+p['numberStims']*p['interStimInterval']+p['staPreFrame']
        self.stim_frames = np.arange(p['stimStartFrame'],p['stimStopFrame']+1,p['interStimInterval'])
                
        # count processed frames
        self.framesProc = 0
        
        # FLAGS
        self.STOP_JOB_FLAG = None
        self.UseONACID = p['UseONACID']
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
        
        # caiman related 
        self.com_count =  self.c['K']
        
        # make Temp folder 
        self.save_dir = p['temp_save_dir']
        


#%% process frame
    def work(self):
        # local counts
        sta_stim_idx = 0
        sta_frame_idx = 0
        LastPlot = 0
        framesProc = self.framesProc
        com_count = self.com_count
        refreshFrame = p['refreshFrame']

        # local param
        FLAG_USE_ONACID = self.UseONACID
        BufferPointer = self.BufferPointer
        sta_start_frames = self.sta_start_frames
        photo_stim_frames = self.photo_stim_frames
        stim_frames = self.stim_frames # stim at fixed intervals
        
        # local record of all roi coords
        print('number of initial rois '+str(self.com_count))
        ROIx = np.asarray([item.get("ROIx") for item in self.ROIlist[:self.com_count]]) 
        ROIy = np.asarray([item.get("ROIy") for item in self.ROIlist[:self.com_count]])
        print(ROIx)
        print(ROIy)

        if FLAG_USE_ONACID:
        # use locally scoped variables to speed up
            self.status_signal.emit('Using OnACID')
    
            # Extract initialisation parameters  
            mot_corr = True # hard coded now; using motion correction from onacid
            img_norm = self.c['img_norm'].copy().astype(np.float32)
            img_min = self.c['img_min'].copy().astype(np.float32)
            ds_factor = self.c['ds_factor']
            dims =  self.c['cnm_init'].dims
    
            # Define OnACID parameters
            max_shift = np.ceil(10./ds_factor).astype('int')  # max shift allowed
            cnm_init = self.c['cnm_init']
            cnm2 = deepcopy( cnm_init)
            t_cnm = cnm2.initbatch
            Cn_init = self.c['Cn_init']
            coms_init = self.c['coms_init']
            self.Cn = Cn_init.copy()
            coms = coms_init.copy()
            com_count = self.c['K']
            rejected = 0
            accepted = list(range(1,self.c['K']+1))
            expect_components = True

        # keep processing frames in qbuffer
        while ((not p['FLAG_END_LOADING']) or (not qbuffer.empty())) and not self.STOP_JOB_FLAG:
            
            # get data from queue
            frame_in = qbuffer.get() 
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
                    mot_corr_start = time.time()
                    templ = cnm2.Ab.dot(cnm2.C_on[:cnm2.M, t_cnm - 1]).reshape(cnm2.dims, order='F') * img_norm
                    frame_cor, shift = motion_correct_iteration_fast(frame_in, templ, max_shift, max_shift)
                    self.shifts.append(shift)
                    print(shift)                    
                    print('caiman motion correction time:' + str("%.4f"%(time.time()-mot_corr_start)))
                else:
                    frame_cor = frame_in
    
                frame_cor = frame_cor / img_norm                            # normalize data-frame
                cnm2.fit_next(t_cnm, frame_cor.reshape(-1, order='F'))      # run OnACID on this frame
                
                # detect new compunents 
                if expect_components:
                    update_comp_time = time.time()    
                    if cnm2.N - (com_count+rejected) == 1:
                        self.new_coms = com(cnm2.Ab[:, -1], dims[0], dims[1])[0]
    
                        # Check for repeated components
                        close = abs(coms - self.new_coms) < self.dist_thresh/ds_factor
                        repeated = any(np.all(close,axis=1))
    
                        if repeated == False: # add component to ROI
                            coms = np.vstack((coms, self.new_coms))
                            y, x = self.new_coms   # reversed
                            ROIx = np.append(ROIx,x*ds_factor)
                            ROIy = np.append(ROIy,y*ds_factor)
    
                            com_count += 1
                            accepted.append(cnm2.N)
                            print('New cell detected (' + str(cnm2.N-rejected) + ')')
                            
                            # add new ROI to photostim target, if required
                            if p['FLAG_BLINK_CONNECTED'] and p['FLAG_AUTO_ADD_TARGETS']:
                                p['currentTargetX'] = np.append(p['currentTargetX'],x*ds_factor)
                                p['currentTargetY'] = np.append(p['currentTargetY'],y*ds_factor)
                                self.sendCoords_signal.emit()
                                self.updateTargetROIs_signal.emit()
                            
                            self.getROImask_signal.emit(x,y) # add roi coords to list in main
                            print('add new component time:' + str("%.4f"%(time.time()-update_comp_time)))
                            if com_count == p['MaxNumROIs']:
                                expect_components = False
                                print('Not accepting more components')
                        else:
                            print('Repeated component found')
                            rejected += 1
                          
                
                # add data to buffer
                try:
                    self.RoiBuffer[:com_count, BufferPointer] = cnm2.C_on[accepted,t_cnm] # cnm2.noisyC is without deconvolution
                    self.ROIlist_threshold[:com_count] = np.nanmean(self.RoiBuffer[:com_count,:], axis=1) + 3*np.nanstd(self.RoiBuffer[:com_count,:], axis=1)
                except Exception as e:
                    print(e)
                    logger.exception(e)
                    print(self.RoiBuffer[:com_count,:])
                
                # trigger photostim

                if p['photoProtoInx'] == CONSTANTS.PHOTO_ABOVE_THRESH:
                    photostim_idx = self.RoiBuffer[:com_count, BufferPointer]-self.ROIlist_threshold[:com_count]
                    print(photostim_idx)
                    p['currentTargetX'] = ROIx[photostim_idx>0]
                    p['currentTargetY'] = ROIy[photostim_idx>0] 
                    print(p['currentTargetX'])
                    print(p['currentTargetX'])
                    self.sendCoords_signal.emit()
                    self.sendPhotoStimTrig_signal.emit()
                    self.updateTargetROIs_signal.emit() 
                    
                elif p['photoProtoInx'] == CONSTANTS.PHOTO_BELOW_THRESH:
                    photostim_idx = self.ROIlist_threshold[:com_count] - self.RoiBuffer[:com_count, BufferPointer]
                    p['currentTargetX'] = ROIx[photostim_idx>0]
                    p['currentTargetY'] = ROIy[photostim_idx>0]                 
                    self.sendCoords_signal.emit()
                    self.sendPhotoStimTrig_signal.emit()
                    self.updateTargetROIs_signal.emit() 
                    
                
                
                    
                # trigger sta recording 
                if sta_stim_idx <self.num_stims:
                    if p['FLAG_STIM_TRIG_ENABLED'] and framesProc == stim_frames[sta_stim_idx]: # send TTL 
                        self.sendTTLTrigger_signal.emit()
                    if p['photoProtoInx'] == CONSTANTS.PHOTO_FIX_FRAMES and framesProc == photo_stim_frames[sta_stim_idx]:
                        self.sendPhotoStimTrig_signal.emit()
                    
                    if not self.flag_sta_recording:
                        if framesProc == sta_start_frames[sta_stim_idx]:
                            self.flag_sta_recording = True
                            sta_frame_idx = 0
                    else:
                        self.sta_traces[:com_count,sta_stim_idx,sta_frame_idx] = cnm2.C_on[accepted,t_cnm]
                        sta_frame_idx +=1
                        if sta_frame_idx == self.sta_trace_length:
                            self.flag_sta_recording = False
                            sta_stim_idx += 1
                
                if framesProc>refreshFrame-1: #frame_count>self.BufferLength-1:
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
                                denoise_time = time.time()                            
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
                                print('generate denoise frame time:' + str("%.4f"%(time.time()-denoise_time)))
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
                print('display time'+str("%.4f"%(time.time()-temp_time)))
    
            self.tottime.append(time.time() - t0)                             # store time for each frame
            print('process frame time:' + str("%.4f"%(time.time()- t0)))
            print('frames proc =' + str(framesProc))
            self.frame_signal.emit(framesProc)

            qbuffer.task_done()

#%% post-loop finishing   
        print('finishing work')
        self.com_count = com_count
        self.BufferPointer = BufferPointer
        self.cnm2 = cnm2
        self.status_signal.emit('Cumulative processing time is ' + str(np.nanmean(self.tottime)) + ' sec.')
        if self.UseONACID:
            self.BufferPointer = 0
            # show indices on viewbox
            self.showROIIdx_signal.emit()
            
            
            # save results to Temp folder 
            if self.FLAG_PV_CONNECTED:
                self.movie_name = os.path.basename(self.pl.get_movie_name())
            else:
                self.movie_name = os.path.splitext(p['moviePath'])[0]
            
            save_dict = dict()
            save_dict['cnm2'] = self.cnm2
            save_dict['accepted'] = accepted
            save_dict['t_cnm'] = t_cnm
            save_object(save_dict, self.movie_name + '_DS_' + str(ds_factor) + '_OnlineProc.pkl')
    

            # save sta traces to Temp folder - havent testing in live stream mode
            sta_trial_avg = np.nanmean(self.sta_traces,1)
            if(p['staAvgStopFrame']>p['staAvgStartFrame']):
                sta_tiral_avg_amp = np.nanmean(sta_trial_avg[:,p['staAvgStartFrame']+p['staPreFrame']:p['staAvgStopFrame']+p['staPreFrame']],1)
                self.sta_amp_signal.emit(sta_tiral_avg_amp[:self.com_count])
            else:
                self.status_signal.emit('Check sta average range')
            save_name = str(self.save_dir) + 'sta_traces.npy'
            print('sta traces saved as: '+save_name)
            np.save(save_name,self.sta_traces)
            
            
            # transfer data to main
            self.transDataToMain_signal.emit(sta_trial_avg, self.cnm2)
            
            # finishing          
            self.finished_signal.emit()

#    def photostimTargets(self, coords):
        
    def stop(self):
        self.status_signal.emit('Stop clicked')
        self.STOP_JOB_FLAG = True

        
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
        self.movie_path = 'U:/simulate_movie/20170823.tif'
        self.ref_movie_path = 'U:/simulate_movie/20170823_ref.tif'
        self.movie = np.atleast_3d(self.BlankImage)
        self.movie = np.swapaxes(self.movie,0,2)   # TODO: SWAPPED axes of the movie before. not anymore. should it?
        
        self.opsin_img_path = 'U:/simulate_movie/20170823_Ch01.tif'
        self.opsin_img = np.array([])
        self.opsin_mask = np.array([])
        self.A_opsin = np.array([])
        self.opsinMaskOn = False
        self.A_loaded = np.array([])
        self.mask_path = ''
        
        # ROI selection --not used
        self.drawOn = False
        self.RoiCenter = QPoint(244,244)  #TODO: UNNECESSARY? but maybe it is
        self.RoiRadius = 10
        self.thisROI = pg.CircleROI(self.RoiCenter,self.RoiRadius*2)  # TODO: UNNECESSARY?

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
        
        self.Targetpen = pg.mkPen(color = (255,112,75), width = 3, style = Qt.DotLine)
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
        self.staPlot_gridLayout.addWidget(self.figCanvas)

        # caiman values
        self.c = {}
        
        # buffer for sta traces
        self.sta_traces = np.array([])
        self.sta_amp = None
        
        # photostim target list
        self.numTargets = 0
        self.bl = []
        self.TargetIdx = np.array([],dtype ='uint16')# indices in ROIlist
        self.TargetX = np.array([]) # all target coords
        self.TargetY = np.array([])
        p['ExtraTargetX'] = np.array([]) #selected targets (not in the ROIlist) -- TODO
        p['ExtraTargetY'] = np.array([])
        p['currentTargetX'] = np.array([]) # keep track of what is currently on SLM
        p['currentTargetY'] = np.array([])
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
        
        # get gui elements
        self.getValues()
        
        # configuration in GUI
        self.trial_config = {}
        
        # daq config
        self.daq_array = np.ones((99), dtype = np.float)*5
        self.daq_array = np.append(self.daq_array,0)
        self.niStimTask = daqtask.Task()
        self.niPhotoStimTask = daqtask.Task()
        
        try: # sensory stim
            self.niStimTask.ao_channels.add_ao_voltage_chan(p['stimDaqDevice'])
            self.niStimWriter = stream_writers.AnalogSingleChannelWriter(self.niStimTask.out_stream,True)
            self.niStimWriter.write_one_sample(0,10.0)
        except Exception as e:
            print(str(e))
            logger.exception(e)
            self.updateStatusBar(str(e))
            time.sleep(2)
            
        try: # photo stim
            self.niPhotoStimTask.ao_channels.add_ao_voltage_chan(p['photostimDaqDevice'])
            self.niPhotoStimWriter= stream_writers.AnalogSingleChannelWriter(self.niPhotoStimTask.out_stream,True)
            self.niPhotoStimWriter.write_one_sample(0,10.0)
        except Exception as e:
            print(str(e))
            self.updateStatusBar(str(e))
            time.sleep(2)
            
        # make threads
        self.workerThread = MyThread()
        self.streamThread = MyThread()

        
        # flag reading/streaming data
        p['FLAG_END_LOADING'] = False
        
        # power control
        p['power_polyfit_p'] = get_power_params()
        
        # signal/slot connections
        self.setConnects()
        
        # methods for ImageWindow when opsinMaskOn
        self.ImageWindow.enterEvent = self.displayOpsinImg.__get__(self.ImageWindow)
        self.ImageWindow.leaveEvent = self.displayOpsinMask.__get__(self.ImageWindow)

        
    def random_color(self):
        r = random.randrange(0, 255)
        g = random.randrange(0, 255)
        b = random.randrange(0, 255)
        return QColor(r, g, b)       
            
    def initialiseROI(self):
        self.ALL_ROI_SELECTED = False
        self.thisROIIdx = 0
#        self.MaxNumROIs = self.MaxNumROIs_spinBox.value()  # TODO: move elsewhere?
        NumROIs =p['MaxNumROIs']
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
        
    def getMousePosition(self,event):
        pointSize = self.RoiRadius*2+5
        x = event.pos().x()
        y = event.pos().y()
        print(str(x)+' '+str(y))
        
        selected = self.Targetcontour_item.data
        selected_x = [value[0] for value in selected]
        selected_y = [value[1] for value in selected]
        
        det_dist = 20  # detection distance
        detected = 0
        
        # unselect a target if distance to other target <= det_dist
        for target in selected:
            detected += abs(x - target[0]) <= det_dist and abs(y - target[1]) <= det_dist
            if detected:
                index = np.argwhere(selected==target)[0][0]
                self.Targetcontour_item.clear()
                del selected_x[index]
                del selected_y[index]
                self.Targetcontour_item.addPoints(x = selected_x, y = selected_y, pen = self.Targetpen, size = pointSize)
                
                #   works too but less efficient:
#                    points = list(range(len(selected_x)))
#                    points.remove(index)
#                    for point in points:
#                        self.Targetcontour_item.addPoints(x = [selected_x[point]], y = [selected_y[point]], pen = self.Targetpen, size = pointSize)
                break

        if not detected:
            self.Targetcontour_item.addPoints(x = [x], y = [y], pen = self.Targetpen, size = pointSize)
            p['ExtraTargetX'] = np.append(p['ExtraTargetX'],x)
            p['ExtraTargetY'] = np.append(p['ExtraTargetY'],y)
            print('selected x = ' + str(x))
            print('selected y = ' + str(y))
    
    def showMousePosition(self,event):
        x = event.x()
        y = event.y()  
        self.xPos_label.setText(str(x))
        self.yPos_label.setText(str(y))
    
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
        
        # start worker
        self.run_pushButton.clicked.connect(self.clickRun)
        
        # display
        self.plotOn_checkBox.clicked.connect(self.switch_plotOn)
        self.displayOn_checkBox.clicked.connect(self.switch_displayOn)
        self.denoiseOn_checkBox.clicked.connect(self.switch_denoiseOn)
        self.plotSTA_pushButton.clicked.connect(self.plotSTA)
        self.plotSTAonMasks_pushButton.clicked.connect(lambda: self.plotSTAonMasks(self.sta_amp))
        self.showComponents_pushButton.clicked.connect(lambda: self.plotSTAonMasks(None))
        
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
        
        # prairie
        
        # others
        self.test_pushButton.clicked.connect(self.tempTest)
        self.reset_pushButton.clicked.connect(self.resetAll)  # can be implemented in the future

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
        

        self.ds_factor = p['dsFactor']
        self.MaxNumROIs = p['MaxNumROIs']
        self.flipRowChanged()
        self.currentZoomChanged()

        #p['InitNumROIs'] = K  #TODO: keep one copy only
    def testTTLTrigger(self):
        self.getValues()
        try:
            self.niStimTask.ao_channels.add_ao_voltage_chan(p['stimDaqDevice'])
            self.niStimWriter = stream_writers.AnalogSingleChannelWriter(self.niStimTask.out_stream,True)
            self.niStimWriter.write_one_sample(0,10.0)
        except Exception as e:
            print(e)
            self.updateStatusBar('All settings updated, Error sending trgger:'+str(e))
        try:
            self.sendTTLTrigger()
            self.updateStatusBar('All settings updated, test trigger sent')
        except Exception as e:
            print(e)
            
    def saveResults(self):
        try:
            save_dict = dict()
            save_dict['cnm2'] = self.cnm2
            save_dict['accepted'] = self.accepted
            save_dict['t_cnm'] = self.t_cnm
            save_object(save_dict, p['saveResultPath'])
        except Exception as e:
            self.updateStatusBar('save results error:'+str(e))
            
    def browseResultPath(self):
        saveResultsPath = QFileDialog.getSaveFileName (self, 'Select save path', p['saveResultPath'], '.pkl')
        p['saveResultPath'] = saveResultsPath[0]+saveResultsPath[1]
        self.saveResultPath_lineEdit.setText(p['saveResultPath'])
        
    def getPathFromPV(self):
        if p['FLAG_PV_CONNECTED']:
            p['saveResultPath'] = self.pl.get_movie_name()
        self.saveResultPath_lineEdit.setText(p['saveResultPath'])    
            
    def testPhotoStimTrigger(self):
        self.getValues()
        try:
            self.niPhotoStimTask.ao_channels.add_ao_voltage_chan(p['photostimDaqDevice'])
            self.niPhotoStimWriter = stream_writers.AnalogSingleChannelWriter(self.niPhotoStimTask.out_stream,True)
            self.niPhotoStimWriter.write_one_sample(0,10.0)
        except Exception as e:
            print(e)
            self.updateStatusBar('All settings updated, Error sending photostim:'+str(e))
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
            self.niPhotoStimWriter.write_many_sample(self.daq_array,10.0)
            print('photostim trigger sent')
        except Exception as e:
            print(e)
          
    def clearTargets(self):
        msg = QMessageBox()
        msg.setText("Clear all targets?") 
        msg.setWindowTitle('pyRTAOI Message')
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)       
        retval = msg.exec_()
        if(retval==QMessageBox.Ok):
            self.Targetcontour_item.clear()
            self.TargetIdx = np.array([],dtype ='uint16')# indices in ROIlist
            self.TargetX = np.array([]) # all target coords
            self.TargetY = np.array([])
            p['ExtraTargetX'] = np.array([]) #selected targets (not in the ROIlist) -- TODO
            p['ExtraTargetY'] = np.array([])
            p['CurrentTargetX'] = np.array([])
            p['CurrentTargetY'] = np.array([])
            self.numTargets = 0
            
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
            self.ROIlist[ROIIdx]["threshold"] = list()   # ADDED THRESHOLD
            self.ROIlist[ROIIdx]["STA"] = list()
        print( self.ROIlist[ROIIdx]["ROIcolor"] )
        
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
               
            # buffer for sta traces
            self.sta_traces = np.array([])
            
            # clear table   
            self.updateTable()
            
            
        msg.setText("RESET CAIMAN??") 
        msg.setWindowTitle('pyRTAOI Message')
        msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)  
        retval = msg.exec_()
        if(retval==QMessageBox.Ok):
            self.c = {}
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
#            print('sending coords, bl connected = '+ str(self.bl.CONNECTED))
            # scale targets coordinates (to refZoom with which the SLM transform matrix was computed)
            currentTargetX = (p['currentTargetX']-255.0)*p['targetScaleFactor']+255.0
            currentTargetY = (p['currentTargetY']-255.0)*p['targetScaleFactor']+255.0
        
            if(self.bl.CONNECTED):                
                            
                if not self.bl.send_coords(currentTargetX.astype(int).tolist(), currentTargetY.astype(int).tolist() ):

                    self.updateStatusBar('Phase mask updated')
                else:
                    self.updateStatusBar('Update phase mask ERROR!')
                    
            else:
                self.updateStatusBar('Send coords faild, check blink connection')
                p['FLAG_BLINK_CONNECTED'] = False
                
    def sendCoordsWorker(self): # - NOT USED
        # called by worker
        if p['FLAG_BLINK_CONNECTED']:
            if(self.bl.CONNECTED):
                if not self.bl.send_coords(p['currentTargetX'].tolist(), p['currentTargetY'].tolist() ):
                    self.updateStatusBar('Phase mask updated')
                else:
                    self.updateStatusBar('Update phase mask ERROR!')
                    
    def selectAllROIs(self):
        # add all ROIlist to the targetlist        
        self.TargetIdx =np.append(self.TargetIdx, range(0,self.thisROIIdx))
        self.TargetIdx = np.unique(self.TargetIdx)
        self.TargetX = np.append(self.TargetX, [self.ROIlist[i]["ROIx"] for i in self.TargetIdx])
        self.TargetY = np.append(self.TargetY, [self.ROIlist[i]["ROIy"] for i in self.TargetIdx])
        self.numTargets = self.TargetX.shape[0]+p['ExtraTargetX'].shape[0] 
        self.updateTargets()    
        self.updateTargetROIs()
    
    def updateTargets(self):
        # save current targets to p
        p['currentTargetX']  = np.append(self.TargetX,p['ExtraTargetX'])
        p['currentTargetY']  = np.append(self.TargetY,p['ExtraTargetY'])   
        self.updateTargetROIs()
            
    def updateTargetROIs(self):
        # redraw targeted rois in imaging window - to do: draw what is saved in p
        self.Targetcontour_item.clear()
        self.Targetcontour_item.addPoints(x = p['currentTargetX'], y = p['currentTargetY'], pen= self.Targetpen, size = self.RoiRadius*2+5)   
        
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
        self.photoProtoInx =self.photoProto_comboBox.currentIndex()
        self.updateStatusBar('Photostim protocol changed to' + str(self.photoProtoInx))
        p['photoProtoInx'] = self.photoProtoInx
          
    def loadMoviePath(self):
        movie_path = str(QFileDialog.getOpenFileName(self, 'Load movie', '', 'MPTIFF (*.tif);;All files (*.*)')[0])
        self.movie_path = movie_path
        p['moviePath'] = movie_path
        if movie_path:
            self.moviePath_lineEdit.setText(movie_path)
            movie_ext = os.path.splitext(p['moviePath'])[1]
            if movie_ext == '.tif':
                print('Correct video found')
                
    def plotSTAonMasks(self,sta_amp):
        # show cells detected by caiman
        # make a image with sta levels
        cnm = self.proc_cnm
        print(cnm.gnb)
        A, b = cnm.Ab[:, cnm.gnb:], cnm.Ab[:, :cnm.gnb].toarray()

        if issparse(A):
            A = np.array(A.todense())
        else:
            A = np.array(A)
 
        d, nr = np.shape(A)
        
        # use sta value, otherwise use one
        if sta_amp is None: # will show scaled amplitude of A
            sta_amp = np.ones((nr,))*255           
        else: # normalise within component before multiply with sta
            for i in range(nr):
                A[:,i] = A[:,i]/sum(A[:,i])

        # put value into mask
        cellMask = np.zeros((cnm.dims[1]*cnm.dims[0],))

        for i in range(np.minimum(len(sta_amp),nr)):
            if not np.isnan(sta_amp[i]):
                cellMask+=A[:,i].flatten()*sta_amp[i]
                print(max(A[:,i]))
                print('sum =' + str(sum(A[:,i])))

        cellMask2D = np.reshape(cellMask,cnm.dims,order='F')
        cellMask2D = cellMask2D/max(cellMask)*255

        norm = plt.Normalize(0,1)
        im = plt.imshow(norm(cellMask2D),aspect = 'equal',cmap = 'Greys')
        plt.colorbar(im, orientation='horizontal')
        plt.show()
        
        cellMask2D =  np.repeat(cellMask2D[:,:,None],3,axis=-1)
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
        sta_tiral_avg_amp = np.nanmean(sta_trial_avg[:,p['staAvgStartFrame']+p['staPreFrame']:p['staAvgStopFrame']+p['staPreFrame']],1)

        
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
            axs[i].set_title('roi'+str(i)+' amp = '+str("%.2f"%sta_tiral_avg_amp[i]))
            axs[i].set_aspect('auto')
            axs[i].set_frame_on(False)
            axs[i].set_adjustable('box')
            
        self.sta_amp = sta_tiral_avg_amp
        print(self.sta_amp.shape)
        
#        plt.show()
        #  add figure to layout
        self.staTracesFig = fig
        self.staTracesCanvas = FigureCanvas(self.staTracesFig)
        self.staPlot_gridLayout.addWidget(self.staTracesCanvas)
        self.updateStatusBar('STA traces plotted on plot tab')
    
    
    def resetFigure(self):
        plt.close('all')
        
        self.ax1.cla()
        self.ax2.cla()
        self.ax3.cla()
        self.axcomp.cla()
        
#            self.fig.clear()
        self.figCanvas.draw()
#            self.figCanvas.update()
        self.figCanvas.flush_events()
        
        g = gc.collect()
        print('gc: ', g)

        
    def showOnacidResults(self, cnm_struct, t, img=None, plot=True):
        if plot:
            self.plotOnacidTraces(cnm_struct, t, img)
        self.c['cnm2'] = cnm_struct
        self.checkOpsin()
        
        
    def plotOnacidTraces(self, cnm_struct, t, img=None):

        try:
            self.resetFigure()
        except:
            pass
            
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
        
        plt.ion()
        if 'csc_matrix' not in str(type(A)):
            A = csc_matrix(A)
        if 'array' not in str(type(b)):
            b = b.toarray()
    
        nr, T = C.shape
        nb = f.shape[0]
    
        Y_r = YrA + C
        d1, d2 = cnm_struct.dims
    
        if img is None:
            img = np.reshape(np.array(A.mean(axis=1)), (d1, d2), order='F')
        
        self.axcomp = self.fig.add_axes([0.05, 0.05, 0.9, 0.03])
        self.ax1 = self.fig.add_axes([0.05, 0.55, 0.4, 0.4])
        self.ax3 = self.fig.add_axes([0.55, 0.55, 0.4, 0.4])
        self.ax2 = self.fig.add_axes([0.05, 0.1, 0.9, 0.4])
        self.s_comp = Slider(self.axcomp, 'Component', 0, nr + nb - 1, valinit=0)
        
        vmax = np.percentile(img, 95)
        
        def update(val):
            i = np.int(np.round(self.s_comp.val))
            print(('Component:' + str(i)))
    
            if i < nr:
    
                self.ax1.cla()
                imgtmp = np.reshape(A[:, i].toarray(), (d1, d2), order='F')
                self.ax1.imshow(imgtmp, interpolation='None', cmap=pl.cm.gray, vmax=np.max(imgtmp)*0.5)
                self.ax1.set_title('Spatial component ' + str(i + 1))
                self.ax1.axis('off')
    
    
                self.ax2.cla()
                self.ax2.plot(np.arange(T), Y_r[i], 'c', linewidth=3)
                self.ax2.plot(np.arange(T), C[i], 'r', linewidth=2)
                self.ax2.set_title('Temporal component ' + str(i + 1))
                self.ax2.legend(labels=['Filtered raw data', 'Inferred trace'])

    
                self.ax3.cla()
                self.ax3.imshow(img, interpolation='None', cmap=pl.cm.gray, vmax=vmax)
                imgtmp2 = imgtmp.copy()
                imgtmp2[imgtmp2 == 0] = np.nan
                self.ax3.imshow(imgtmp2, interpolation='None',
                           alpha=0.5, cmap=pl.cm.hot)
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
                
                
        def arrow_key_image_control(self, event):  # this doesn't work now
#            print(event.key())
            if event.key() == Qt.Key_Left: # event.key == 'left':
                new_val = np.round(self.s_comp.val - 1)
                if new_val < 0:
                    new_val = 0
                self.s_comp.set_val(new_val)
    
            elif event.key() == Qt.Key_Right: # event.key == 'right':
                new_val = np.round(self.s_comp.val + 1)
                if new_val > nr + nb:
                    new_val = nr + nb
                self.s_comp.set_val(new_val)
            else:
                pass
            
        
        self.s_comp.on_changed(update)
        self.s_comp.set_val(0)
        self.figCanvas.mpl_connect('key_release_event', arrow_key_image_control)  # TODO: implement this. event should be input to function
        self.figCanvas.draw()
#            self.figCanvas.update()
        self.figCanvas.flush_events()
        self.updateStatusBar('OnACID traces plotted in plot tab')


    def loadRefMoviePath(self):
        if self.UseOpsinMask_radioButton.isChecked() or self.UseOnacidMask_radioButton.isChecked():
            ref_movie_path = str(QFileDialog.getOpenFileName(self, 'Load reference', '', 'MPTIFF (*.tif);;All files (*.*)')[0])
        else:
            ref_movie_path = str(QFileDialog.getOpenFileName(self, 'Load reference', '', 'PKL(*.pkl);;MPTIFF (*.tif);;All files (*.*)')[0])
        
        self.ref_movie_path = ref_movie_path
        self.refMoviePath_lineEdit.setText(self.ref_movie_path)
        print('Reference movie selected')
        #self.initialiseCaiman()
        
        
    def loadMask(self):
        mask_path = str(QFileDialog.getOpenFileName(self, 'Load reference', '', 'NPZ (*.npz);;PKL (*.pkl);;MPTIFF (*.tif);;All files (*.*)')[0])
        
        self.mask_path = mask_path
        if self.mask_path:
            loaded_mask_ext = os.path.splitext(self.mask_path)[1]
        
        if loaded_mask_ext == '.tif' or loaded_mask_ext == '.pkl':
            print('Correct file with a mask selected')
        
        if loaded_mask_ext == '.npz': # faster than pkl for loading mask
            data_loaded = np.load(self.mask_path)
            self.A_loaded = data_loaded['A']
            dims = data_loaded['dims']
            
        elif loaded_mask_ext == '.pkl':
            data_loaded = load_object(self.mask_path)
            A_loaded = data_loaded['cnm_init'].A
            self.A_loaded = np.array(A_loaded.todense())
            dims = data_loaded['cnm_init'].dims

        
        if self.A_loaded.size:
            ds_factor = self.dsFactor_doubleSpinBox.value()
            expected_dim = int(512/ds_factor)
            equal_size = list(dims) == [expected_dim for i in range(2)]
            if not equal_size:
                self.updateStatusBar('Mask dimension mismatch. Change downsampling to ' + str(512/dims[0])[:4] +
                                     ' or select a different mask')
                return
                
            self.opsinMaskOn = False
            self.dims = dims
            loaded_mask = np.reshape(np.array(self.A_loaded.mean(axis=1)), dims, order='F')
            self.updateImage(cv2.resize(np.squeeze(loaded_mask), (512, 512), interpolation=cv2.INTER_CUBIC))
            self.updateStatusBar('Mask with ' + str(self.A_loaded.shape[-1]) + ' cells was loaded')    
            self.UseOnacidMask_radioButton.setEnabled(True)
            self.UseOnacidMask_radioButton.setChecked(True)
        else:
            self.updateStatusBar('No mask was not found in the file selected')
            self.UseOnacidMask_radioButton.setChecked(False)
    
    
    def initialiseCaiman(self):
        if self.ref_movie_path:
            movie_ext = os.path.splitext(p['refMoviePath'])[1]
            self.updateStatusBar('Initialising...')
        else:
            self.updateStatusBar('No reference movie selected!')
            return
            
        K = self.InitNumROIs_spinBox.value()
        ds_factor = self.dsFactor_doubleSpinBox.value()
        
        # init method to be used
        opsin_seeded = self.UseOpsinMask_radioButton.isChecked()
        mask_seeded = self.UseOnacidMask_radioButton.isChecked()
        cnmf_init = not (opsin_seeded or mask_seeded)
                
        # process image or load from pkl file directly
        if cnmf_init:
            if movie_ext == '.tif':
                K = self.InitNumROIs_spinBox.value()
    
                print('Starting CNMF initialisation with tiff file')
                lframe, init_values = initialise(self.ref_movie_path, init_method='cnmf', K=K, initbatch=500, ds_factor=ds_factor,
                                                 rval_thr=0.85, thresh_overlap=0, save_init=True, mot_corr=True, merge_thresh=0.85,
                                                 decay_time=0.2, min_SNR=2.5)
                                                     # rval_thr for component quality
                                                 
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
                    lframe, init_values = initialise(self.ref_movie_path, init_method='seeded', Ain=self.A_opsin,
                                                     ds_factor=ds_factor, initbatch=500, rval_thr=0.85,
                                                     thresh_overlap=0.1, save_init=True, mot_corr=True,
                                                     merge_thresh=0.85, CNN_filter=True)
                    
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
        
        # Extract number of cells detected
        K = cnm2.N
        self.InitNumROIs_spinBox.setValue(K)
        print('Number of components initialised: ' + str(K))
        
        if self.MaxNumROIs_spinBox.value() < K+5:
            self.MaxNumROIs_spinBox.setValue(K+10)
        
        print(cnm2.N)
        self.ds_factor = ds_factor
        self.c = init_values
        self.c['cnm2'] = cnm2
        
        if idx_components is not None:
            coms = self.c['coms_init']
            self.c['coms_init'] = coms[idx_components]
            
#        if opsin_seeded: self.c['cnm2'].opsin = [True]*len(init_values['idx_components'])
        
        self.proc_cnm = init_values['cnm2'] #['cnm_init']
        self.dims = init_values['cnm_init'].dims
        self.InitNumROIs = K
        self.opsinMaskOn = False
        self.imageItem.setImage(cv2.resize(self.c['img_norm'],(512,512),interpolation=cv2.INTER_CUBIC)) # display last frame of the initiialisation movie

        if self.A_opsin.size:
            print('Checking opsin overlap')
            self.checkOpsin()
        
        ### add a step to discard unwanted components from initialisation?
        
        self.initialiseROI()
        for i in range(K):
            y, x = init_values['coms_init'][i]  # reversed
            self.getROImask(thisROIx = x, thisROIy = y)
    
        self.updateStatusBar('Initialision completed')
        show_traces = 1
        if show_traces:
            cnm_init = init_values['cnm_init']
            self.plotOnacidTraces(cnm_struct=cnm2,t=cnm_init.initbatch,img=init_values['Cn_init'])


# THIS ALLOWS FOR DRAWING and initialises ROIs -- unnecessary now. unless you want to press reset initially #
#    def drawROI(self, thisROI):
#        # self.initialiseROI()
#        def addROI(self):
#            #  if self.drawOn:
#            self.graphicsView.addItem(thisROI)#(self.thisROI)
#        #self.drawOn = not self.drawOn
#        if self.drawOn:
#            addROI(self)

    
    def loadOpsinImg(self):
        opsin_img_path = str(QFileDialog.getOpenFileName(self, 'Load a C1V1 image', '', 'MPTIFF (*.tif)')[0])
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
            radius = self.cellRadius_spinBox.value()
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
        else:
            A = Ain
            
        dims = self.dims # tuple([np.sqrt(A.shape[0])]*2)
        overlap = []
        opsin = []
        self.opsin_overlap_ratio = self.overlapRatio_doubleSpinBox.value()
        print(self.opsin_overlap_ratio)
        
        # Convert A to numpy array
        if issparse(A):
            A = np.array(A.todense())
        else:
            A = np.array(A)
        
        onacid_mask = (deepcopy(A)>0).astype('int')  # binarise the onacid output mask
        
        for cell in range(onacid_mask.shape[-1]):
            cell_mask = (np.reshape(onacid_mask[:,cell], dims, order='F'))
            cell_pix = sum(sum(cell_mask == 1))
            
            inter = cv2.bitwise_and(self.opsin_mask, cell_mask)
            inter_pix = sum(sum(inter))
            cell_overlap = inter_pix/cell_pix
            overlap.append(cell_overlap)
            
            if cell_overlap <= self.opsin_overlap_ratio:
                onacid_mask[:,cell][onacid_mask[:,cell] == 1] = -3
            else:
                onacid_mask[:,cell][onacid_mask[:,cell] == 1] = 3
                
            opsin.append(cell_overlap > self.opsin_overlap_ratio)
        
        if Ain is None:
            # store info on opsin
            self.c['cnm2'].opsin = opsin
            self.c['overlap'] = overlap
            self.c['opsin_mask'] = self.opsin_mask
            self.c['opsin_overlap_ratio'] = self.opsin_overlap_ratio
            self.onacid_mask = onacid_mask
        return onacid_mask
        
            
    def showCellsOnMask(self):
        if self.c != {}:
            self.checkOpsin()  # for overlap update

            # visualise all comps
            try:
                summed_A = np.hstack((self.A_opsin, self.onacid_mask))
                summed_mask = np.reshape(np.array(summed_A.sum(axis=1)), self.dims, order='F')
                pl.figure();pl.imshow(summed_mask)
                pl.colorbar()
            except Exception as e:
                print(e)
            
        elif self.A_loaded.size:
            loaded_mask = self.checkOpsin(self.A_loaded)
            
            summed_A = np.hstack((self.A_opsin, loaded_mask))
            summed_mask = np.reshape(np.array(summed_A.sum(axis=1)), self.dims, order='F')
            pl.figure();pl.imshow(summed_mask)
            pl.colorbar()
            
        else:
            self.updateStatusBar('No cells to be shown!')


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
                self.ROIcontour_item.addPoints(x = [thisROIx], y = [thisROIy], pen = pg.mkPen(qcolor, width=2), size = self.RoiRadius*2) #TODO: how is this plotted??

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
        if(self.selectAll_checkBox.isChecked):
            self.selectAllROIs()
            
    def autoAddClicked(self):
        p['FLAG_AUTO_ADD_TARGETS'] = self.addNewROIsToTarget_checkBox.isChecked
        print('FLAG_AUTO_ADD_TARGETS ='+str(p['FLAG_AUTO_ADD_TARGETS']))
        
    def enterEvent(self,event):
        self.graphicsView.setCursor(Qt.CrossCursor)
        

    def clickRun(self):
        self.resetFigure()
        
        self.getValues()
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
        self.workerObject.sendTraces_signal.connect(self.showOnacidResults)
        
        # start and finish
        self.workerObject.transDataToMain_signal.connect(self.transDataToMain)
        self.workerThread.started.connect(self.workerObject.work)
        self.workerObject.finished_signal.connect(self.workerThread.exit)
        
        # triggers
        self.workerObject.sendTTLTrigger_signal.connect(self.sendTTLTrigger)
        self.workerObject.sendCoords_signal.connect(self.sendCoords)
        self.workerObject.sendPhotoStimTrig_signal.connect(self.sendPhotoStimTrigger)
        self.stop_pushButton.clicked.connect(self.stop_thread)
        self.workerThread.start()
        self.updateStatusBar('Worker started')
        
        
        # start stream thread
        kwargs = {"prairie": self.pl}
        self.streamObject = DataStream(**kwargs)
        self.streamObject.moveToThread(self.streamThread)
        self.streamThread.started.connect(self.streamObject.stream)
        self.streamObject.stream_finished_signal.connect(self.streamThread.exit)
        self.streamThread.start()
        self.updateStatusBar('stream started')


    def refreshPlot(self,arr):
        self.plotItem.clear()
        for i in range(self.thisROIIdx):
            self.plotItem.plot(arr[i,:], antialias=True, pen=pg.mkPen(self.ROIlist[i]["ROIcolor"], width=1))#self.ROIpen)   #roipen width for this as well?
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
        p['plotOn'] = self.plotOn_checkBox.isChecked

    def switch_displayOn(self):
        p['displayOn'] = self.displayOn_checkBox.isChecked

    def switch_denoiseOn(self):
        p['denoiseOn'] = self.denoiseOn_checkBox.isChecked
        
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
        
    def stop_thread(self):
        self.workerObject.stop()
        self.workerThread.wait(1)
        self.workerThread.quit()
        # self.workerThread.wait()
        
    def transDataToMain(self, sta_amp,proc_cnm): # to do - 
        self.sta_amp = sta_amp
        self.proc_cnm = proc_cnm # processed cnm
                
        
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
            self.workerObject.stop()
            self.workerThread.wait(1)
            self.workerThread.quit()
            
        if self.streamObject:
            self.streamObject.stop()
            self.streamThread.wait(1)
            self.streamThread.quit()
        
        # reset AO to zero
        self.niStimWriter.write_one_sample(0,10.0)
        self.niStimTask.stop()
        self.niStimTask.close()
        del self.niStimTask
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
    
   # global data buffer
    global qbuffer 
    qbuffer = Queue(maxsize = 0)

    # launch program
    try:
        main(sys.argv)
    except Exception as e:
        logger.exception(e)
        print(str(e))
