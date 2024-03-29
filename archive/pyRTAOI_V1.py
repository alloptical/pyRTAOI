''''
 Dependences:
 follow instructions in these links to install all dependences
 caiman  : https://github.com/flatironinstitute/CaImAn 
 opencv  (pip install .whl) : https://www.lfd.uci.edu/~gohlke/pythonlibs/#opencv
 PyDAQmx (pip) : https://github.com/clade/PyDAQmx




 TO DO:
1. log - photostimulation targets and frame idx; save online analysis result
    offline analysis - use cnmf in matlab?
 

2. deal with zoom - check
3. mark point setting for changing photostim power; if not fast enough use dummy targets in zero blocker
 
 
 
4. Plot STA dF/F - send triggers - check
5. threshold sometimes becomes Inf - check
6. auto initialisation - take reference movie and load 
7. delete target by right clicking 

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
from skimage.external import tifffile
import numpy as np
import pyqtgraph as pg
import time
import pylab as pl
import cv2
from copy import deepcopy
import matplotlib.pyplot as plt
from scipy.sparse import issparse, spdiags, coo_matrix, csc_matrix
import pickle
import logging
import json

# prairie link
import pvlink
from pvlink import *
import ctypes 
from ctypes import *

# PYCUDA - imorted on the worker thread
#import pycuda.driver as cuda
#from pycuda.tools import make_default_context
#import pycuda.gpuarray as gpuarray
#from pycuda.compiler import SourceModule

# blink
from bLink import bLink

if sys.platform == 'win32':
    time_ = time.clock
else:
    time_ = time.time

# caiman libraries
import caiman as cm
from caiman_func.initialisation import initialise
from caiman.motion_correction import motion_correct_iteration_fast
from caiman.base.rois import com
from caiman.utils.utils import load_object
# plot libs
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib

# ni max
# pydaqmx badly documented, use nidaqmx instead
import nidaqmx
import nidaqmx.task as daqtask
from nidaqmx import stream_writers


# configure logging
logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)
handler = logging.FileHandler('errors.log')  # create a file handler
handler.setLevel(logging.DEBUG)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')  # create a logging format
handler.setFormatter(formatter)
logger.addHandler(handler)  # add the handlers to the logger
logger.info('Started application')

    
# Ensure using PyQt5 backend
matplotlib.use('QT5Agg')

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


    def __init__(self, **kwargs ):
        super(Worker,self).__init__()
        
        # get parameters from Main thread
        self.c = kwargs.get('caiman',{})
        self.pl = kwargs.get('prairie',[])
        print('worker pl empty? ' + str(self.pl == []))
        
        # sta traces
        self.sta_traces = np.empty([p['MaxNumROIs'],p['numberStims'],p['staPreFrame']+p['staPostFrame']])
        self.sta_traces.fill(np.nan)
        p['stimStopFrame'] = p['stimStartFrame']+p['numberStims']*p['interStimInterval']+p['staPreFrame']
        self.stim_frames = np.arange(p['stimStartFrame'],p['stimStopFrame']+1,p['interStimInterval'])
#

    def work(self):     
        self.STOP_JOB_FLAG = None
        self.UseONACID = p['UseONACID']
        self.FLAG_PV_CONNECTED = p['FLAG_PV_CONNECTED']
        self.FLAG_OFFLINE = p['FLAG_OFFLINE']
        print('pv connected? ' + str(self.FLAG_PV_CONNECTED))
        print('offline analysis' + str(self.FLAG_OFFLINE))
        
        if self.c == {}:
            self.UseONACID = False
            print('No reference movie provided')
        if self.pl == []:
            print('No Prairie object provided')
        
        # initialise sta 
        sta_start_frames = self.stim_frames -p['staPreFrame'] 
        photo_stim_frames = self.stim_frames + p['photoWaitFrames']
        flag_sta_recording = False
        sta_stim_idx = 0
        sta_frame_idx = 0
        sta_trace_length = p['staPreFrame']+p['staPostFrame']
        num_stims = p['numberStims']
        
        # get current parameters
        BufferLength = p['BufferLength']
        RoiMeanBuffer = p['RoiMeanBuffer']
        BufferPointer = p['BufferPointer']
        NumROIs = p['MaxNumROIs']
        ROIlist_threshold = np.zeros(NumROIs)

        LastPlot = 0       
        bufferFrames = p['pvBufferFrames']
        refresh = p['refreshFrame']
        
#%% load onACID pre-settings
        if self.UseONACID:
            self.status_signal.emit('Using OnACID')

            # Extract initialisation parameters
            # mot_corr = self.c['mot_corr']
            mot_corr = True # hard coded now; using motion correction from on acid
            img_norm = self.c['img_norm'].copy().astype(np.float32)
            img_min = self.c['img_min'].copy().astype(np.float32)
            T1 = self.c['T1']
            ds_factor = self.c['ds_factor']
            Cn_init = self.c['Cn_init']
            cnm_init = self.c['cnm_init']
            #gnb =  self.c['cnm_init'].gnb
            dims =  self.c['cnm_init'].dims

            # cells detected during initialisation
            K = self.c['K']
            coms_init = self.c['coms_init']

            # Define OnACID parameters
            max_shift = np.ceil(10./ds_factor).astype('int')  # max shift allowed
            cnm2 = deepcopy(cnm_init)
            t = cnm2.initbatch
            tottime = []
            Cn = Cn_init.copy()
            shifts = []
            coms = coms_init.copy()
            com_count = K
            rejected = 0
            accepted = list(range(1,K+1))
            expect_components = True
            dist_thresh = 21
            display_shift = False  # option to display motion correction shift
            
        
#%% movie streaming from prairie
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
            BUFFER_SIZE = bufferFrames*samplesPerFrame
            thisID = os.getpid()
            buf = np.ndarray(BUFFER_SIZE, dtype = c_int16)
            buf_add, read_only_flag = buf.__array_interface__['data']
            _rrd = make_command('-rrd',str(thisID),str(c_int64(buf_add).value),str(BUFFER_SIZE))
            self.status_signal.emit('cpu buffer created')
            
            # assign destination buffer in gpu
            dest = np.ndarray(pixelsPerLine*linesPerFrame, dtype = 'int32')
            dest_gpu = gpuarray.to_gpu(dest.astype(np.int32))
            self.status_signal.emit('gpu initiated')
            
            # initialise before start scanning
            self.FLAG_SCANNING  = False
            loop = 0
            framesRecvd = 0
            recvTime = 0
            convtTime = 0
            incorrectCount = 0
            incorrectLimit = 1000
            
            # start imaging
            num_samples = self.pl.send_recv(_rrd) # flush
            if not self.pl.send_done(self.pl._ts): # t-series
                self.FLAG_SCANNING = True
                self.status_signal.emit('Started scanning')
                
            while self.FLAG_SCANNING:
                loop += 1
                # print('This is loop '+ str(loop))
                loop_start_time = time.time()
                num_samples = self.pl.send_recv(_rrd) # kernel dies here on the second read; sometimes first read 20180322
            
                if not num_samples:
                    self.FLAG_SCANNING  = False
                    
                elif int(num_samples) == samplesPerFrame:    

                    print('recv time: '+str("%.4f"%(time.time()-loop_start_time)))
                    recvTime += time.time()-loop_start_time
                    recv_time = time.time()
                    
                    # -----  use cuda -------
                    # try:
                    buf_gpu = gpuarray.to_gpu(buf.astype(np.int16))
                    print('from cpu to gpu:' + str("%.4f"%(time.time()-recv_time)))
                    transfer_time = time.time()
                    sample_mean( buf_gpu, pixelsPerLine, linesPerFrame, samplesPerPixel, p['flipEvenRows'], dest_gpu,
                        block=(1024,1,1), grid = (int(samplesPerFrame/1024),1,1) )
                    print('convert time:' + str("%.4f"%(time.time()-transfer_time)))
                    transfer_time = time.time()
                    dest_gpu = dest_gpu.reshape(linesPerFrame,pixelsPerLine)
                    thisFrame = dest_gpu.get()
                    print('gpu to cpu:' + str("%.4f"%(time.time()-transfer_time)))
            
                    # except Exception as e:
                    #     print('error gpu convert: '+str(e))
                    #     break
                    convtTime += time.time()-recv_time
                    
                    temp_time = time.time()
                    frame_ = thisFrame.copy().astype(np.float32)
                    framesRecvd += 1
                    self.frame_signal.emit(framesRecvd)
                    
                    if self.UseONACID: # HAVENT TESTED ONLINE YET
                        # move trace buffer pointer
                        if BufferPointer==BufferLength-1:
                            BufferPointer = 0
                        else:
                            BufferPointer +=1
                            
                        # process current frame
                        if ds_factor > 1:
                            frame_ = cv2.resize(frame_, img_norm.shape[::-1])   # downsampling
        
                        frame_ -= img_min                                       # make data non-negative
        
                        if mot_corr:                                            # motion correct
                            templ = cnm2.Ab.dot(cnm2.C_on[:cnm2.M, t - 1]).reshape(cnm2.dims, order='F') * img_norm
                            frame_cor, shift = motion_correct_iteration_fast(frame_, templ, max_shift, max_shift)
                            shifts.append(shift)
                        else:
                            templ = None
                            frame_cor = frame_
        
                        frame_cor = frame_cor / img_norm                        # normalize data-frame
                        cnm2.fit_next(t, frame_cor.reshape(-1, order='F'))      # run OnACID on this frame
                        
                        # detect new compunents 
                        if expect_components:
                            if cnm2.N - (com_count+rejected) == 1:
                                new_coms = com(cnm2.Ab[:, -1], dims[0], dims[1])[0]
        
                                # Check for repeated components
                                close = abs(coms - new_coms) < dist_thresh/ds_factor
                                repeated = any(np.all(close,axis=1))
        
                                if repeated == False: # add component to ROI
                                    coms = np.vstack((coms, new_coms))
                                    y, x = new_coms   # reversed
        
                                    com_count += 1
                                    accepted.append(cnm2.N)
                                    print('New cell detected (' + str(cnm2.N-rejected) + ')')
                                    
                                    # add new ROI to photostim target, if required
                                    if p['FLAG_BLINK_CONNECTED'] and p['FLAG_AUTO_ADD_TARGETS']:
                                        p['currentTargetX'] = np.append(p['currentTargetX'],x*ds_factor)
                                        p['currentTargetY'] = np.append(p['currentTargetY'],y*ds_factor)
                                        self.sendCoords_signal.emit()
                                        self.updateTargetROIs_signal.emit()
                                    
                                    self.getROImask_signal.emit(x,y)
        
                                    if com_count == NumROIs:
                                        expect_components = False
                                        print('Not accepting more components')
                                else:
                                    print('Repeated component found')
                                    rejected += 1
                                    
                        
                        # display current frame
                        if p['denoiseOn']:
                            A, b = cnm2.Ab[:, cnm2.gnb:], cnm2.Ab[:, :cnm2.gnb].toarray()
                            C_t, f_t = cnm2.C_on[cnm2.gnb:cnm2.M, t], cnm2.C_on[:cnm2.gnb, t]
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
                        else:
                            if ds_factor == 1:
                                self.roi_signal.emit(frame_cor)
                            else:
                                res_frame_cor = cv2.resize(frame_cor, (512, 512), interpolation=cv2.INTER_CUBIC)
                                self.roi_signal.emit(res_frame_cor)
                                
        
                        shift_ = [round(shift[0]), round(shift[1])]
                        if shift_ != [0,0]:
                            if display_shift == True: #denoiseOn == False:
                                self.refreshScatter_signal.emit(shift_)
                                
                        # add data to buffer
                        try:
                            RoiMeanBuffer[:com_count, BufferPointer] = cnm2.C_on[accepted,t] #[gnb:gnb+com_count, t]
                            ROIlist_threshold[:com_count] = np.nanmean(RoiMeanBuffer[:com_count,:], axis=1) + 3*np.nanstd(RoiMeanBuffer[:com_count,:], axis=1)
                        except Exception as e:
                            print(e)
                            logger.exception(e)
                            print(RoiMeanBuffer[:com_count,:])
                            
                        # sta recording 
                        if sta_stim_idx <num_stims:
                            if p['FLAG_STIM_TRIG_ENABLED'] and framesRecvd == self.stim_frames[sta_stim_idx]: # send TTL 
                                self.sendTTLTrigger_signal.emit()
                            if p['photoProto'] == CONSTANTS.PHOTO_FIX_FRAMES and framesRecvd == photo_stim_frames[sta_stim_idx]:
                                self.sendPhotoStimTrig_signal.emit()
                            
                            if not flag_sta_recording:
                                if framesRecvd == sta_start_frames[sta_stim_idx]:
                                    flag_sta_recording = True
                                    sta_frame_idx = 0
                            else:
                                self.sta_traces[:com_count,sta_stim_idx,sta_frame_idx] = cnm2.C_on[accepted,t]
                                sta_frame_idx +=1
                                if sta_frame_idx == sta_trace_length:
                                    flag_sta_recording = False
                                    sta_stim_idx += 1
                        
                        if framesRecvd>refresh-1: #frame_count>BufferLength-1:
                            if LastPlot == refresh:
                                if p['plotOn']:
                                    self.refreshPlot_signal.emit(RoiMeanBuffer[:com_count,:])
                                if p['displayOn']:
                                    self.thresh_signal.emit(ROIlist_threshold[:com_count])
                                LastPlot = 0
                            else:
                                LastPlot += 1
        
                        elif framesRecvd == refresh-1 and p['plotOn']: #frame_count == BufferLength-1
                            self.status_signal.emit('Plotting trace')
                            
                    # else not using onACID        
                    else: 
                        self.roi_signal.emit(frame_)
                        print('display time'+str("%.4f"%(time.time()-temp_time)))
                        print('sample received:'+num_samples[0:-2])
                        
                        
                    # trigger photostim
                
                    QApplication.processEvents()
                    lp_time = time.time()-loop_start_time
                    print('loop:'+str(loop)+' time:' + str("%.4f"%(lp_time))+"\n\r")
                    incorrectCount = 0
                    
                else: 
                    # print('sample number incorrect, samples recved =' + num_samples)
                    incorrectCount +=1
                if incorrectCount>incorrectLimit: 
                    print('too many incorrect frames, quit loop')
                    self.FLAG_SCANNING   = False
                    self.context.pop()
                    del self.context
                    from pycuda.tools import clear_context_caches
                    clear_context_caches()
                    del buf_gpu
                    del dest_gpu
                    del sample_mean
            # ADD FINISHING HERE
            

#%% OFFLINE
        if self.FLAG_OFFLINE and self.UseONACID:
            self.status_signal.emit('Running OnACID')
            
            # load movie
            print('Loading video')
            Y_ = cm.load(p['moviePath'], subindices=slice(0,T1,None))
            print('Video loaded')

            # old_comps = cnm2.N                              # number of existing components

            # MAIN LOOP
            for frame_count, frame in enumerate(Y_):        # now process each file
                t0 = time_()
                
                if self.STOP_JOB_FLAG:
                    STOP_JOB_FLAG  = False
                    break

                if BufferPointer==BufferLength-1:
                    BufferPointer = 0
                else:
                    BufferPointer +=1
                self.frame_signal.emit(frame_count)

                if np.isnan(np.sum(frame)):
                    raise Exception('Frame ' + str(frame_count) + ' contains nan')
#                    if t % 100 == 0:
#                        print('Epoch: ' + str(iter+1) + '. ' + str(t-cnm_init.initbatch)+' frames have beeen processed in total. '+str(cnm2.N - old_comps)+' new components were added. Total number of components is '+str(cnm2.Ab.shape[-1]-gnb))
#                        old_comps = cnm2.N

                frame_ = frame.copy().astype(np.float32)
                
                # process current frame
                if ds_factor > 1:
                    frame_ = cv2.resize(frame_, img_norm.shape[::-1])   # downsampling

                frame_ -= img_min                                       # make data non-negative

                if mot_corr:                                            # motion correct
                    templ = cnm2.Ab.dot(cnm2.C_on[:cnm2.M, t - 1]).reshape(cnm2.dims, order='F') * img_norm
                    frame_cor, shift = motion_correct_iteration_fast(frame_, templ, max_shift, max_shift)
                    shifts.append(shift)
                else:
                    templ = None
                    frame_cor = frame_
                frame_cor = frame_cor / img_norm                        # normalize data-frame
                cnm2.fit_next(t, frame_cor.reshape(-1, order='F'))      # run OnACID on this frame

                if expect_components:
                    if cnm2.N - (com_count+rejected) == 1:
                        new_coms = com(cnm2.Ab[:, -1], dims[0], dims[1])[0]

                        # Check for repeated components
                        close = abs(coms - new_coms) < dist_thresh/ds_factor
                        repeated = any(np.all(close,axis=1))

                        if repeated == False:
                            coms = np.vstack((coms, new_coms))
                            y, x = new_coms   # reversed

                            com_count += 1
                            accepted.append(cnm2.N)
                            print('New cell detected (' + str(cnm2.N-rejected) + ')')
                                                        
                            # add new ROI to photostim target, if required
                            if p['FLAG_BLINK_CONNECTED'] and p['FLAG_AUTO_ADD_TARGETS']:
									 # sending coords
                                print('updating targets')
                                p['currentTargetX'] = np.append(p['currentTargetX'],x*ds_factor)
                                p['currentTargetY'] = np.append(p['currentTargetY'],y*ds_factor)
								
                                self.sendCoords_signal.emit()
                                self.updateTargetROIs_signal.emit()
								
                            self.getROImask_signal.emit(x,y)

                            if com_count == NumROIs:
                                expect_components = False
                                print('Not accepting more components')
                        else:
                            print('Repeated component found')
                            rejected += 1
                            
                
                # display current frame
                if p['denoiseOn']:
                    A, b = cnm2.Ab[:, cnm2.gnb:], cnm2.Ab[:, :cnm2.gnb].toarray()
                    C_t, f_t = cnm2.C_on[cnm2.gnb:cnm2.M, t], cnm2.C_on[:cnm2.gnb, t]
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
                else:
                    if ds_factor == 1:
                        self.roi_signal.emit(frame_cor)
                    else:
                        res_frame_cor = cv2.resize(frame_cor, (512, 512), interpolation=cv2.INTER_CUBIC)
                        self.roi_signal.emit(res_frame_cor)
                        

                shift_ = [round(shift[0]), round(shift[1])]
                if shift_ != [0,0]:
                    if display_shift == True: #denoiseOn == False:
                        self.refreshScatter_signal.emit(shift_)
                        
                # add data to buffer
                # TODO: every now and then there is an overflow here
                try:
                    RoiMeanBuffer[:com_count, BufferPointer] = cnm2.C_on[accepted,t] #[gnb:gnb+com_count, t]
                    ROIlist_threshold[:com_count] = np.nanmean(RoiMeanBuffer[:com_count,:], axis=1) + 3*np.nanstd(RoiMeanBuffer[:com_count,:], axis=1)
                except Exception as e:
                    print(e)
                    logger.exception(e)
                    print(RoiMeanBuffer[:com_count,:])
                    
                # sta recording
                if sta_stim_idx <num_stims:
                    # send triggers, if required
                    if p['FLAG_STIM_TRIG_ENABLED'] and framesRecvd == self.stim_frames[sta_stim_idx]: # send TTL 
                        self.sendTTLTrigger_signal.emit()
                    if p['photoProto'] == CONSTANTS.PHOTO_FIX_FRAMES and framesRecvd == photo_stim_frames[sta_stim_idx]:
                        self.sendPhotoStimTrig_signal.emit()
                    
                    if flag_sta_recording == False:
                        if frame_count == sta_start_frames[sta_stim_idx]:
                            flag_sta_recording = True
                            sta_frame_idx = 0
                    else:
                        self.sta_traces[:com_count,sta_stim_idx,sta_frame_idx] = cnm2.C_on[accepted,t]
                        sta_frame_idx +=1
                        if sta_frame_idx == sta_trace_length:
                            flag_sta_recording = False
                            sta_stim_idx += 1
                
#                    if displayOn:
#                        self.thresh_signal.emit(ROIlist_threshold)  #TODO: maybe display every X frames, with refresh plot?
                if frame_count>refresh-1: #frame_count>BufferLength-1:
                    if LastPlot == refresh:
                        if p['plotOn']:
                            self.refreshPlot_signal.emit(RoiMeanBuffer[:com_count,:])
                        if p['displayOn']:
                            self.thresh_signal.emit(ROIlist_threshold[:com_count])
                        LastPlot = 0
                    else:
                        LastPlot += 1

                elif frame_count == refresh-1 and p['plotOn']: #frame_count == BufferLength-1
                    self.status_signal.emit('Plotting trace')
        
                QApplication.processEvents()
                t += 1
                tottime.append(time_() - t0)                             # store time for each frame

            print('Cumulative processing speed is ' + str(t/np.sum(tottime))[:5] + ' frames per second.')
            print(t/np.sum(tottime))
            
            # show indices on viewbox
            self.showROIIdx_signal.emit()
            
            # save sta traces to Temp folder
            sta_trial_avg = np.nanmean(self.sta_traces,1)
            if(p['staAvgStopFrame']>p['staAvgStartFrame']):
                sta_tiral_avg_amp = np.nanmean(sta_trial_avg[:,p['staAvgStartFrame']+p['staPreFrame']:p['staAvgStopFrame']+p['staPreFrame']],1)
                self.sta_amp_signal.emit(sta_tiral_avg_amp[:com_count])
            else:
                self.status_signal.emit('Check sta average range')
#            save_name = os.path.join(os.path.dirname(__file__),'/Temp/sta_traces.npy')
            save_name = '/Temp/sta_traces.npy'
            print('sta trace save name: '+save_name)
            np.save(save_name,self.sta_traces)
            
            
            #% Optionally save results
            save_results = False
            if save_results:
                np.savez('results_analysis_online_MOT_CORR.npz',
                         Cn=Cn, Ab=cnm2.Ab, Cf=cnm2.C_on, b=cnm2.b, f=cnm2.f,
                         dims=cnm2.dims, tottime=tottime, noisyC=cnm2.noisyC, shifts=shifts)

            BufferPointer = 0
            frame_count = 0
            self.finished_signal.emit()
            
#%% FINISHING WORK
        # save sta traces to Temp folder - havent testing in live stream mode
        sta_trial_avg = np.nanmean(self.sta_traces,1)
        if(p['staAvgStopFrame']>p['staAvgStartFrame']):
            sta_tiral_avg_amp = np.nanmean(sta_trial_avg[:,p['staAvgStartFrame']+p['staPreFrame']:p['staAvgStopFrame']+p['staPreFrame']],1)
            self.sta_amp_signal.emit(sta_tiral_avg_amp[:com_count])
        else:
            self.status_signal.emit('Check sta average range')
#            save_name = os.path.join(os.path.dirname(__file__),'/Temp/sta_traces.npy')
        save_name = '/Temp/sta_traces.npy'
        print('sta trace save name: '+save_name)
        np.save(save_name,self.sta_traces)
            
        # transfer data to main
        self.transDataToMain_signal.emit(sta_trial_avg, cnm2)
        
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
        self.pvBufferFrames = 1 # worked for 25
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
        self.RoiMeanBuffer = np.empty([self.MaxNumROIs,self.BufferLength])
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
        self.ds_factor = 1.5

        # display traces
        self.plotItem = self.plotArea_graphicsLayoutWidget.addPlot(antialias=True)
        self.plotArea_graphicsLayoutWidget.setAntialiasing(True)
        self.plotItem.disableAutoRange()
        self.plotItem.setXRange(0,self.BufferLength-1)

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
            
        # make worker thread
        self.workerThread = MyThread()
        kwargs = {"caiman": self.c,"prairie": self.pl}
        self.workerObject = Worker(**kwargs)
        
        # signal/slot connections
        self.setConnects()
        

        
    def random_color(self):   # TODO: maybe make it non random but varied depending on expected comps
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
        self.RoiMeanBuffer = np.empty([NumROIs,self.BufferLength])

        # dictionary of ROIs        
        self.resetROIlist()

        # threshold table in GUI
        self.updateTable()

        # clear scatter plot
        self.ROIcontour_item.clear()
        
    def updateTable(self):
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
        self.Targetcontour_item.addPoints(x = [x], y = [y], pen= self.Targetpen, size = pointSize)   
        p['ExtraTargetX'] = np.append(p['ExtraTargetX'],x)
        p['ExtraTargetY'] = np.append(p['ExtraTargetY'],y)
        print('selected x ='+str(x))
        print('selected y ='+ str(y))

    def showMousePosition(self,event):
        x = event.x()
        y = event.y()  
        self.xPos_label.setText(str(x))
        self.yPos_label.setText(str(y))
        
    def tempTest(self):
        print('button clicked')
        self.plotSTAonMasks()       

    def setConnects(self):
        # load and save
        self.loadMoviePath_pushButton.clicked.connect(self.loadMoviePath)
        self.loadRefMoviePath_pushButton.clicked.connect(self.loadRefMoviePath)
        self.saveConfig_pushButton.clicked.connect(self.saveConfig)
        self.loadConfig_pushButton.clicked.connect(self.loadConfig)
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
                        # save parameters
                        
        # parameters not in widgets                
        p['BufferLength'] = self.BufferLength
        p['RoiMeanBuffer'] = self.RoiMeanBuffer
        p['BufferPointer'] = self.BufferPointer
        p['FLAG_PV_CONNECTED'] = self.FLAG_PV_CONNECTED
        p['FLAG_OFFLINE'] = self.IsOffline
        p['photoProtoInx'] = self.photoProtoInx
        
 
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
        self.niPhotoStimWriter.write_many_sample(self.daq_array,10.0)
        print('photostim trigger sent')
		  
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


    def loadRefMoviePath(self):
        ref_movie_path = str(QFileDialog.getOpenFileName(self, 'Load reference', '', 'PKL(*.pkl);;MPTIFF (*.tif);;All files (*.*)')[0])
        self.ref_movie_path = ref_movie_path
        self.initialiseCaiman()
    
    def initialiseCaiman(self):
        if self.ref_movie_path:
            self.refMoviePath_lineEdit.setText(self.ref_movie_path)
            movie_ext = os.path.splitext(p['refMoviePath'])[1]
            
        K = self.InitNumROIs_spinBox.value()
        print('Number of components to initialise: ' + str(K))
        ds_factor = self.dsFactor_doubleSpinBox.value()
        self.ds_factor = ds_factor
        
        # process image or load from pkl file directly
        if movie_ext == '.tif':
            print('Starting initialisation with tiff file')
            lframe, init_values = initialise(self.ref_movie_path, init_method='cnmf', K=K, initbatch=500, ds_factor = ds_factor,
                                             rval_thr = 0.85, thresh_overlap = 0,save_init = True)
                                             # rval_thr for component quality
        if movie_ext == '.pkl': # something was not saved/loaded properly; did not work for actual processing
            print('loading initialisation file')
            init_values = load_object(self.ref_movie_path)
                # Prepare object for OnACID
            init_values['cnm_init']._prepare_object(np.asarray(init_values['Yr']), init_values['T1'], init_values['expected_comps'], idx_components=None,
                             min_num_trial = 2, N_samples_exceptionality = int(init_values['N_samples']))
                             
        self.c = init_values
        self.proc_cnm = init_values['cnm_init']
        self.InitNumROIs = K
        self.imageItem.setImage(cv2.resize(self.c['img_norm'],(512,512),interpolation=cv2.INTER_CUBIC)) # display last frame of the initiialisation movie

        ### add a step to discard unwanted components from initialisation? and move below elsewhere/button pressed ###
        self.initialiseROI()
        for i in range(K):
            y, x = init_values['coms_init'][i]  # reversed
            self.getROImask(thisROIx = x, thisROIy = y)
        
############ THIS ALLOWS FOR DRAWING and initialises ROIs -- unnecessary now. unless you want to press reset initially ##############
    def drawROI(self, thisROI):
        # self.initialiseROI()
        def addROI(self):
            #  if self.drawOn:
            self.graphicsView.addItem(thisROI)#(self.thisROI)
        #self.drawOn = not self.drawOn
        if self.drawOn:
            addROI(self)


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
                #self.Thresh_tableWidget.setItem(0,self.thisROIIdx,self.Thresh_labels[self.thisROIIdx]) # threshold table
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
        for i in range(self.thisROIIdx):
                # TODO: uses Qpoints so bear the radius/pos offset in mind
            self.ROIcontour_item.addPoints(x = [self.ROIlist[i]["ROIx"] + shift[0]],
                                           y = [self.ROIlist[i]["ROIy"] + shift[1]],
                                           pen = pg.mkPen(self.ROIlist[i]["ROIcolor"], width=2),
                                           size = self.RoiRadius*2) #self.ROIlist[i]["ROIradius"][0]*2)

                #TODO: can add multiple points at once??
                                               #but what about their colors then. + is it gonna save much time?//is it worth it?
                                               # potentially can add list of pens??
            # BUT THIS REQUIRES CHANGING THE WAY THE VALUES ARE STORED. CURRENTLY CAN'T EASILY ACCESS ALL XS AND YS

        #self.ROIcontour_item.addPoints(x = [x], y = [y], pen = color, size = radius*2)
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
        self.getValues()
        kwargs = {"caiman": self.c,"prairie": self.pl}
        self.workerObject = Worker(**kwargs)
        self.updateStatusBar('Worker created')        
        self.workerObject.moveToThread(self.workerThread)
        
        self.workerObject.status_signal.connect(self.updateStatusBar)
        
        # update display
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
        
        # triggers
        self.workerObject.sendTTLTrigger_signal.connect(self.sendTTLTrigger)
        self.workerObject.sendCoords_signal.connect(self.sendCoords)
        self.workerObject.sendPhotoStimTrig_signal.connect(self.sendPhotoStimTrigger)
        self.stop_pushButton.clicked.connect(self.stop_thread)
        self.workerThread.start()
        self.updateStatusBar('Worker started')
        

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
        #self.Thresh_labels[ROIIdx].setText(str('{0:.1f}'.format(thresh))) #OLD
    def updateSTA(self, sta_arr):
        for i in range(0,len(sta_arr)):
            print('updating sta value')
#            self.STA_labels[i].setForeground(QBrush(self.ROIlist[i]["ROIcolor"]))
            self.STA_labels[i].setText(str('{0:.1f}'.format(sta_arr[i])))
            self.Thresh_tableWidget.setItem(i,3,self.STA_labels[i]) 

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
            self.workObject.stop()
        if self.workerThread.isAlive:
            self.workerThread.wait(1)
            self.workerThread.quit()
        
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
    # launch program
    try:
        main(sys.argv)
    except Exception as e:
        logger.exception(e)
        print(str(e))

