# Dependences:
# follow instructions in the links to install all dependences
# caiman  : https://github.com/flatironinstitute/CaImAn 
# opencv  (pip install .whl) : https://www.lfd.uci.edu/~gohlke/pythonlibs/#opencv
# PyDAQmx (pip) : https://github.com/clade/PyDAQmx




# TO DO:
# record cell traces once they are found
# Plot STA dF/F - send triggers




import sys
import random
import os

# Qt
import GUI
from PyQt5.QtCore import Qt,QObject, pyqtSignal, QThread, QTimer, QRectF, QUrl,QPoint, QRect, QSize,QPointF,QSizeF, QSettings, QCoreApplication
from PyQt5.QtWidgets import (QComboBox, QCheckBox, QLineEdit, QSpinBox, QLabel,
                             QDoubleSpinBox, QFileDialog, QApplication,
                             QDesktopWidget, QMainWindow, QMessageBox, QTableWidgetItem)
from PyQt5.QtGui import QColor, QIcon, QPalette, QDesktopServices, QImage, QPixmap, QPainter,QGraphicsScene,QGraphicsEllipseItem,QPen,QBrush
from PyQt5 import QtGui, QtCore
# utilities
from skimage.external import tifffile
import numpy as np
import pyqtgraph as pg
import time
import pylab as pl
import cv2
from copy import deepcopy
import matplotlib.pyplot as plt

# prairie link
import pvlink
from pvlink import *
import ctypes 
from ctypes import *

# PYCUDA
#import pycuda.driver as cuda
#from pycuda.tools import make_default_context
#import pycuda.gpuarray as gpuarray
#from pycuda.compiler import SourceModule

# blink
import bLink

if sys.platform == 'win32':
    time_ = time.clock
else:
    time_ = time.time

# caiman libraries
import caiman as cm
from caiman_func.initialisation import initialise
from caiman.motion_correction import motion_correct_iteration_fast
from caiman.base.rois import com

# plot libs
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib

# Ensure using PyQt5 backend
matplotlib.use('QT5Agg')



# Interpret image data as row-major instead of col-major
pg.setConfigOptions(imageAxisOrder = 'row-major')
class MyThread(QThread):
    def run(self):
        self.exec_()

class Worker(QObject):
    status_signal            = pyqtSignal(str,name ='statusSignal')
    roi_signal               = pyqtSignal(object,name = 'roiSignal')
    refreshPlot_signal       = pyqtSignal(object,name = 'refreshPlot')
    thresh_signal            = pyqtSignal(object,name = 'threshSignal')
    sta_amp_signal           = pyqtSignal(object,name = 'staAmpSignal')
    frame_signal             = pyqtSignal(object,name = 'frameSignal')
    refreshScatter_signal    = pyqtSignal(object)

    #drawROI_signal = pyqtSignal(object, name = 'drawROISignal')
    getROImask_signal        = pyqtSignal(object, object, name = 'getROImaskSignal')
    finished_signal          = pyqtSignal()


    def __init__(self, **kwargs ):
        super(Worker,self).__init__()
        self.c = kwargs.get('caiman',{})
        self.pl = kwargs.get('prairie',[])
        print('worker pl not empty? ' + str(self.pl == []))
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
        print('pv connected? ' + str(self.FLAG_PV_CONNECTED))
        
        if self.c == {}:
            self.UseONACID = False
            print('No reference movie provided')
        if self.pl == []:
            print('No Prairie object provided')
        
        # initialise sta 
        sta_start_frames = self.stim_frames -p['staPreFrame']
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
        # print('Maximum number of ROIs to be found: ' + str(NumROIs))
        ROIlist_threshold = np.zeros(NumROIs)

        plotOn = p['plotOn']
        displayOn = p['displayOn']
        denoiseOn = p['denoiseOn']
        LastPlot = 0
        
        bufferFrames = p['pvBufferFrames']
        flipEvenRows = p['flipEvenRows']
        
        # test movie streaming from prairie
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
                    sample_mean( buf_gpu, pixelsPerLine, linesPerFrame, samplesPerPixel, flipEvenRows, dest_gpu,
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
                    self.roi_signal.emit(frame_)
                    print('display time'+str("%.4f"%(time.time()-temp_time)))
                    print('sample received:'+num_samples[0:-2])
                    framesRecvd += 1
            
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


        if self.UseONACID:
            self.status_signal.emit('Running OnACID')

            # Extract initialisation parameters
            mot_corr = self.c['mot_corr']
            img_norm = self.c['img_norm']
            img_min = self.c['img_min']
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
                            self.getROImask_signal.emit(x,y)

                            if com_count == NumROIs:
                                expect_components = False
                                print('Not accepting more components')
                        else:
                            print('Repeated component found')
                            rejected += 1
                            
                
                # display current frame
                if denoiseOn:
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

                # TODO: every now and then there is an overflow here
                try:
                    RoiMeanBuffer[:com_count, BufferPointer] = cnm2.C_on[accepted,t] #[gnb:gnb+com_count, t]
                    ROIlist_threshold[:com_count] = np.nanmean(RoiMeanBuffer[:com_count,:], axis=1) + 3*np.nanstd(RoiMeanBuffer[:com_count,:], axis=1)
                except Exception as e:
                    print(e)
                    print(RoiMeanBuffer[:com_count,:])
                    
                # sta recording
                if sta_stim_idx <num_stims:
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
                refresh = 5 # define this somewhere in the gui
                if frame_count>refresh-1: #frame_count>BufferLength-1:
                    if LastPlot == refresh:
                        if plotOn:
                            self.refreshPlot_signal.emit(RoiMeanBuffer[:com_count,:])
                        if displayOn:
                            self.thresh_signal.emit(ROIlist_threshold[:com_count])
                        LastPlot = 0
                    else:
                        LastPlot += 1

                elif frame_count == refresh-1 and plotOn: #frame_count == BufferLength-1
                    self.status_signal.emit('Plotting trace')
        
                #time.sleep(0.33)  # 1: delay for debugging
                QApplication.processEvents()

                t += 1
                tottime.append(time_() - t0)                             # store time for each frame

            print('Cumulative processing speed is ' + str(t/np.sum(tottime))[:5] + ' frames per second.')
            print(t/np.sum(tottime))
            
            # save sta traces to Temp folder
            sta_trial_avg = np.nanmean(self.sta_traces,1)
            if(p['staAvgStopFrame']>p['staAvgStartFrame']):
                sta_tiral_avg_amp = np.nanmean(sta_trial_avg[:,p['staAvgStartFrame']+p['staPreFrame']:p['staAvgStopFrame']+p['staPreFrame']],1)
                self.sta_amp_signal.emit(sta_tiral_avg_amp[:com_count])
            else:
                self.status_signal.emit('Check sta average range')
            np.save(os.path.join(os.path.dirname(__file__),'/Temp/sta_traces.npy'),self.sta_traces)
            #% Optionally save results
            save_results = False
            if save_results:
                np.savez('results_analysis_online_MOT_CORR.npz',
                         Cn=Cn, Ab=cnm2.Ab, Cf=cnm2.C_on, b=cnm2.b, f=cnm2.f,
                         dims=cnm2.dims, tottime=tottime, noisyC=cnm2.noisyC, shifts=shifts)

            BufferPointer = 0
            frame_count = 0
            self.finished_signal.emit()

    def stop(self):
        self.status_signal.emit('Stop clicked')
        self.STOP_JOB_FLAG = True

class myGraphicsScene(QtGui.QGraphicsScene):
    def __init__(self,parent = None):
        QtGui.QGraphicsScene.__init__(self)
    
        self.Targetpen = QPen(QColor(255,112,75))
        self.Targetpen.setStyle(Qt.DashDotLine)
        self.Targetpen.setWidth(2)
        
        def mouseClickEvent(self, event): #UNFINISHED!        
            pointSize = self.RoiRadius*2+5
    #        globalPoint = QPoint(event.pos().x(),event.pos().y())
           # print('xy: '+str(event.pos().x())+str(event.pos().y()))
    #        localPoint = self.graphicsView.mapFromGlobal(event.pos())
    #        print('localPoint pos = '+str(localPoint.pos()))
            print('mouse clicked')
            x = event.scenePos().x()
            y = 512-event.scenePos().y()
            print(str(x)+' '+str(y))
            self.addEllipse(x-pointSize/2, y -pointSize/2, pointSize,pointSize,pen= self.Targetpen )   

    
class MainWindow(QMainWindow, GUI.Ui_MainWindow):
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
            

        # make worker thread
        self.workerThread = MyThread()

        # signal/slot connections
        self.setConnects()
        
        
        # prairie default settings - maybe move to a xml and load from there 
        self.PV_IP = '128.40.156.161'  # bruker 2: TCP_IP = '128.40.202.220'
        self.PV_PORT = 1236
        self.flipEvenRows = np.int32(1)
        self.pvBufferFrames = 1 # worked for 25
        self.FLAG_PV_CONNECTED = False
        self.pl = []
        
        # Blink settings
        self.BLINK_IP = '128.40.156.163'
        self.BLINK_PORT = 8888
        self.FLAG_BLINK_CONNECTED = False
        self.bl = []

        # setup image window
        self.ImageWindow.ci.layout.setContentsMargins(0, 0, 0, 0)
        self.ImageWindow.ci.layout.setSpacing(0)
        self.graphicsView = self.ImageWindow.addViewBox(enableMouse=False)
        print(self.graphicsView.screenGeometry())
        
        brush = QBrush()
        
        # add scene to view -- should add items to the scene
        self.graphicScene = myGraphicsScene(self)
        self.graphicScene.setSceneRect(QtCore.QRectF(QPointF(0.0,0.0),QSizeF(512,512)))
        self.imageGraphicsView.setScene(self.graphicScene)

        
        self.imageItem = pg.ImageItem(border='w', invertY=False)
        self.graphicScene.addItem(self.imageItem)
#        self.imageGraphicsView.setAspectLocked(True)  # lock the aspect ratio so pixels are always square
#        self.graphicScene.setRange(QRectF(0, 0, 512, 512), padding=0)  # Set initial view bounds
        

        
        self.Targetpen = pg.mkPen(color = (255,112,75), width = 2, stype = Qt.DotLine)
        self.MouseOffsetX = self.ImageWindow.geometry().getCoords()[0] # this does not give the actual offset ~ 4 pixels away
        self.MouseOffsetY = self.ImageWindow.geometry().getCoords()[1]+self.iconSize().height()
        self.BlankImage = np.zeros(shape = (512,512))
        self.imageItem.setImage(self.BlankImage)
        self.movie_path = 'U:/simulate_movie/20170823.tif'
        self.ref_movie_path = 'U:/simulate_movie/20170823_ref.tif'
        self.movie = np.atleast_3d(self.BlankImage)
        self.movie = np.swapaxes(self.movie,0,2)   # TODO: SWAPPED axes of the movie before. not anymore. should it?
        self.FrameIdx = 0;

        # ROI selection
        self.ROIlist = dict()
        self.drawOn = False
        self.RoiCenter = QPoint(244,244)  #TODO: UNNECESSARY? but maybe it is
        self.RoiRadius = 10
        self.thisROI = pg.CircleROI(self.RoiCenter,self.RoiRadius*2)  # TODO: UNNECESSARY?

        # initialise ROI list
        self.InitNumROIs = 3
        self.NumROIs = 30  # desired number of cells
        self.thisROIIdx = 0
        self.ALL_ROI_SELECTED = False

        # mean intensity of ROI
        self.BufferLength = 200
        self.RoiMeanBuffer = np.empty([self.NumROIs,self.BufferLength])
        self.BufferPointer = 0

        # coordinates of ROI pixels
       # self.ROIxcoor = list()
       # self.ROIycoor = list()
        self.ROIpen = QPen()
        self.ROIpen.setWidthF = 0.5
#        self.ROIscn = QGraphicsScene()
#        self.ROIscnBg = QRectF(QPointF(0.0,0.0),QSizeF(256,256))
#        self.ROImask_graphicsView.scale(1,-1)
#        self.ROImask_graphicsView.setScene(self.ROIscn)
#        self.ROImask_graphicsView.ensureVisible(self.ROIscnBg, xMargin = 0, yMargin = 0)
        self.ROIcontour_item = pg.ScatterPlotItem() # brush for cells on current frame
        self.ROIcontour_item.setBrush(brush.setStyle(0))  # this gives open circles
        self.graphicScene.addItem(self.ROIcontour_item)
        
        self.Targetcontour_item = pg.ScatterPlotItem()
        self.Targetcontour_item.setBrush(brush.setStyle(0))
        self.graphicScene.addItem(self.Targetcontour_item)

        self.plotOn = True
        self.displayOn = True
        self.denoiseOn = True
        self.ds_factor = 1.5

        # display traces
        self.plotItem = self.plotArea_graphicsLayoutWidget.addPlot(antialias=True)
        self.plotArea_graphicsLayoutWidget.setAntialiasing(True)
        self.plotItem.disableAutoRange()
        self.plotItem.setXRange(0,self.BufferLength-1)

        # get gui elements
        self.getValues()

        # caiman values
        self.c = {}
        
        # buffer for sta traces
        self.sta_traces = np.array([])
        
        # photostim target list
        self.TargetIdx = np.array([],dtype ='uint16')# indices in ROIlist
        self.TargetX = np.array([]) # all target coords
        self.TargetY = np.array([])
        self.ExtraTargetX = np.array([]) #selected targets (not in the ROIlist) -- TODO
        self.ExtraTargetY = np.array([])
        self.numTargets = 0
        
        

    def initialiseROI(self):
        self.ALL_ROI_SELECTED = False
        self.thisROIIdx = 0
#        self.NumROIs = self.MaxNumROIs_spinBox.value()  # TODO: move elsewhere?
        NumROIs =p['MaxNumROIs']
        print('Number of ROIs to be tracked: ' + str(NumROIs))
        self.RoiMeanBuffer = np.empty([NumROIs,self.BufferLength])

        # dictionary of ROIs
        self.ROIlist = [dict() for i in range(NumROIs)]

        def color():   # TODO: maybe make it non random but varied depending on expected comps
            r = random.randrange(0, 255)
            g = random.randrange(0, 255)
            b = random.randrange(0, 255)
            return QColor(r, g, b)

        for ROIIdx in range(0,NumROIs):
            self.ROIlist[ROIIdx]["ROIx"]= list()
            self.ROIlist[ROIIdx]["ROIy"] = list()
#            self.ROIlist[ROIIdx]["ROIcenter"] = list()  # or use ROIcenter
#            self.ROIlist[self.thisROIIdx]["ROIradius"]  # needed or not?
            self.ROIlist[ROIIdx]["ROIcolor"] = list()
            self.ROIlist[ROIIdx]["ROIcolor"] = color()
            self.ROIlist[ROIIdx]["threshold"] = list()   # ADDED THRESHOLD
            self.ROIlist[ROIIdx]["STA"] = list()

            
            # currently radius always assumed to be 10

        # threshold table in GUI
        self.Thresh_tableWidget.setRowCount(NumROIs)
        self.Thresh_tableWidget.setColumnCount(4)
        self.Thresh_tableWidget.setHorizontalHeaderLabels(['X','Y','Threshold','STA'])
        self.Thresh_labels = [QTableWidgetItem() for i in range(NumROIs)]
        self.STA_labels = [QTableWidgetItem() for i in range(NumROIs)]
        self.X_labels = [QTableWidgetItem() for i in range(NumROIs)]
        self.Y_labels = [QTableWidgetItem() for i in range(NumROIs)]


        # clear scatter plot
        self.ROIcontour_item.clear()
    

    def setConnects(self):
        self.loadMoviePath_pushButton.clicked.connect(self.loadMoviePath)
        self.loadRefMoviePath_pushButton.clicked.connect(self.loadRefMoviePath)
       # self.reset_pushButton.clicked.connect(self.resetAll)  # can be implemented in the future
        self.run_pushButton.clicked.connect(self.clickRun)
        self.plotOn_checkBox.clicked.connect(self.switch_plotOn)
        self.displayOn_checkBox.clicked.connect(self.switch_displayOn)
        self.denoiseOn_checkBox.clicked.connect(self.switch_denoiseOn)
        self.connectPV_pushButton.clicked.connect(self.connectPV)
        self.selectAll_checkBox.clicked.connect(self.selectAllROIs)
        self.flipRow_comboBox.currentIndexChanged.connect(self.flipRowChanged)
        self.plotSTA_pushButton.clicked.connect(self.plotSTA)
        self.connectBlink_pushButton.clicked.connect(self.connectBlink)
        self.SendCoords_pushButton.clicked.connect(self.sendCoords)

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
        p['BufferLength'] = self.BufferLength
        p['RoiMeanBuffer'] = self.RoiMeanBuffer
        p['BufferPointer'] = self.BufferPointer
        p['FLAG_PV_CONNECTED'] = self.FLAG_PV_CONNECTED
        p['flipEvenRows'] = self.flipEvenRows
        self.NumROIs = p['MaxNumROIs']
        
        #p['InitNumROIs'] = K  #TODO: keep one copy only
    def connectBlink(self):
        self.updateStatusBar('connecting Blink')
        self.bl = bLink(self.BLINK_IP,self.BLINK_PORT)
        print('bl created')
        self.FLAG_BLINK_CONNECTED = True
        
        
    def sendCoords(self):
        if(self.bl.CONNECTED):
            if not self.bl.send_coords(self.ROIlist[:]["ROIx"],self.ROIlist[:]["ROIy"]):
                self.updateStatusBar('Phase mask updated')
            else:
                self.updateStatusBar('Update phase mask ERROR!')
        else:
            self.updateStatusBar('Send coords faild, check blink connection')
                
        
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
        if self.flipRow_comboBox.currentText() is 'Even':
            self.flipEvenRows = np.int32(1)
            print('flip even rows')
        elif  self.flipRow_comboBox.currentText() is 'Odd':
            self.flipEvenRows = np.int32(0)
            print('flip odd rows')
            
            
    def loadMoviePath(self):
        movie_path = str(QFileDialog.getOpenFileName(self, 'Load movie', '', 'MPTIFF (*.tif);;All files (*.*)')[0])
        self.movie_path = movie_path
        p['moviePath'] = movie_path
        if movie_path:
            self.moviePath_lineEdit.setText(movie_path)
            movie_ext = os.path.splitext(p['moviePath'])[1]
            if movie_ext == '.tif':
                print('Correct video found')
                
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
            
        self.staTracesFig = fig
#        plt.show()
        #  add figure to layout
        self.staTracesCanvas = FigureCanvas(self.staTracesFig)
        self.staPlot_gridLayout.addWidget(self.staTracesCanvas)
        self.updateStatusBar('STA traces plotted on plot tab')

    def resetAll(self):
        pass

    def loadRefMoviePath(self):
        K = self.InitNumROIs_spinBox.value()
        ds_factor = self.dsFactor_doubleSpinBox.value()
        self.ds_factor = ds_factor
        print('Number of components to initialise: ' + str(K))

        ref_movie_path = str(QFileDialog.getOpenFileName(self, 'Load movie', '', 'MPTIFF (*.tif);;All files (*.*)')[0])
        self.ref_movie_path = ref_movie_path
        if ref_movie_path:
            self.refMoviePath_lineEdit.setText(ref_movie_path)
            movie_ext = os.path.splitext(p['refMoviePath'])[1]

        if movie_ext == '.tif':
            print('Starting initialisation')
            lframe, init_values = initialise(ref_movie_path, init_method='cnmf', K=K, initbatch=500, ds_factor = ds_factor,
                                             rval_thr = 0.85, thresh_overlap = 0)
                                             # rval_thr for component quality

        self.c = init_values
        self.InitNumROIs = K

        self.imageItem.setImage(lframe) # display last frame of the initiialisation movie
        self.FrameIdx = 0;

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
            self.graphicScene.addItem(thisROI)#(self.thisROI)
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

            def create_pixmap(x,y,radius,qcolor): # added qcolor. changed center to x and y
                center = QPoint(x, y)

                pixmap = QPixmap(512, 512)
                painter = QPainter()
                painter.begin(pixmap)
                thiscolor = qcolor
                painter.setBrush(thiscolor)
                painter.drawEllipse(center,radius,radius)
                #TODO: check. ELLIPSE WITH RADIUSES FOR DRAWING?? do i even need the qpoint?
                painter.end()
#                thisRgb = thiscolor.rgb()
                thisRgb = qcolor.rgb()

                return pixmap, thisRgb

            # used to retrieve original central values. pos different from center. size different from radius (2*radius)
            # now: all values fine so not needed?
#            currentRadius = thisROI.size()/2 #self.thisROI.size()/2
            ##print('confirm correct radius:', currentRadius == 10)  #TODO: MAKE SURE THE RADIUS IS BOTH ROUNDED UP AND DOWN
#            currentCenter = thisROI.pos() #+ currentRadius -- TODO: figure out this position stuff!
            #TODO: currentRadius still needed for plotting casue pos different than center. but here?

            qcolor = self.ROIlist[self.thisROIIdx]["ROIcolor"]


            try:
                pixelmap,color = create_pixmap(thisROIx, thisROIy, self.RoiRadius, qcolor)
            except Exception as e:
                print(e)

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
#            self.ROIlist[self.thisROIIdx]["ROIcenter"] = [thisROIx, thisROIy]
#            self.ROIlist[self.thisROIIdx]["ROIradius"] = currentRadius
            self.ROIlist[self.thisROIIdx]["ROIx"] = thisROIx
            self.ROIlist[self.thisROIIdx]["ROIy"] = thisROIy
            self.thisROIIdx += 1
            if self.thisROIIdx == self.NumROIs:
                self.ALL_ROI_SELECTED = True

        else:
            print('All ROI selected')

        # except:
        #      print('ROI out of FOV')
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
    
    def updateROI(self,img):
        # import pdb
        # pdb.set_trace()
        self.updateImage(img)
        
    def selectAllROIs(self):
        # add all ROIlist to the targetlist
        if(self.selectAll_checkBox.checkState):
            self.TargetIdx =np.append(self.TargetIdx, range(0,self.thisROIIdx))
            self.TargetIdx = np.unique(self.TargetIdx)
            self.TargetX = np.append(self.ExtraTargetX, [self.ROIlist[i]["ROIx"] for i in self.TargetIdx])
            self.TargetY = np.append(self.ExtraTargetY, [self.ROIlist[i]["ROIy"] for i in self.TargetIdx])
            self.numTargets = self.TargetX.shape[0]
            
        self.updateTargetROIs()
            
    def updateTargetROIs(self):
        # redraw targeted rois in imaging window
        self.Targetcontour_item.clear()
        self.Targetcontour_item.addPoints(x = self.TargetX, y = self.TargetY, pen= self.Targetpen, size = self.RoiRadius*2+5)   
    
        # select and draw ROIs
    def enterEvent(self,event):
        self.imageGraphicsView.setCursor(Qt.CrossCursor)
        
    def clickRun(self):
        self.getValues()
        kwargs = {"caiman": self.c,"prairie": self.pl}
        self.workerObject = Worker(**kwargs)
        self.updateStatusBar('Worker created')        
        self.workerObject.moveToThread(self.workerThread)
        self.workerObject.status_signal.connect(self.updateStatusBar)
        self.workerObject.roi_signal.connect(self.updateROI)
        self.workerObject.refreshPlot_signal.connect(self.refreshPlot)
        self.workerObject.frame_signal.connect(self.updateFrameInfo)
        self.workerObject.thresh_signal.connect(self.updateThresh)
        self.workerObject.sta_amp_signal.connect(self.updateSTA)
        self.workerObject.refreshScatter_signal.connect(self.addScatterROI)
        #self.workerObject.drawROI_signal.connect(self.drawROI)  #TODO: REMOVE/CHANGE
        self.workerObject.getROImask_signal.connect(self.getROImask)
        self.workerThread.started.connect(self.workerObject.work)
        self.workerObject.finished_signal.connect(self.workerThread.exit)
        self.stop_pushButton.clicked.connect(self.stop_thread)
        self.workerThread.start()
        self.updateStatusBar('Worker started')

#        self.workerObject.stop()
#        self.workerThread.stop()
#        self.workerThread.wait()

    def refreshPlot(self,arr):
       # t0 = time_()
        self.plotItem.clear()
        for i in range(self.thisROIIdx):
            self.plotItem.plot(arr[i,:], antialias=True, pen=pg.mkPen(self.ROIlist[i]["ROIcolor"], width=1))#self.ROIpen)   #roipen width for this as well?

        try:
            self.plotItem.setYRange(np.round(arr.min())-1,np.round(arr.max())+1)
            
        except:
            self.plotItem.setYRange(0,30)
       # print('plotting took ',time_()-t0)

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
        for i in range(self.thisROIIdx):
            self.STA_labels[i].setText(str('{0:.1f}'.format(sta_arr[i])))
    

    def updateFrameInfo(self,FrameIdx):
        self.CurrentFrameIdx_label.setText(str(FrameIdx))

    def switch_plotOn(self):
        self.plotOn = self.plotOn_checkBox.checkState

    def switch_displayOn(self):
        self.displayOn = self.displayOn_checkBox.checkState

    def switch_denoiseOn(self):
        self.denoiseOn = self.denoiseOn_checkBox.checkState
        
 
    def stop_thread(self):
        self.workerObject.stop()
        self.workerThread.wait(1)
        self.workerThread.quit()
        # self.workerThread.wait()
        
    def closeEvent(self,event):
        # override closeEvent method in Qt
        print('form closing, PV connection = ' + str(self.FLAG_PV_CONNECTED))
        if self.FLAG_PV_CONNECTED:
            del self.pl
            print('pv disconnected')
        if self.bl.CONNECTED:
            self.bl.abort()
            del self.bl
            print('blink disconnected')
            
        self.workObject.stop()
        self.workerThread.quit()
        event.accept()




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
    # launch program
    global p
    p = {}
    main(sys.argv)

