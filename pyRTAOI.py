''''
 Dependences:
 follow instructions in these links to install all dependences
 caiman  : https://github.com/flatironinstitute/CaImAn
 opencv  (pip install .whl) : https://www.lfd.uci.edu/~gohlke/pythonlibs/#opencv (caiman installs it)
 nidaqwx (pip install)
 pyqtgraph (pip install)

SPYDER HAS NO AUTOSAVE FUNCTION!  (only autosaves when script is run)
CAREFUL WHEN USING 'REPLACE ALL' - IT WILL QUIT, WITHOUT SAVING!!

UPDATE REV40 PV SETTINGS MAX VOLTAGE FOR PHOTOSTIM CONTROL!!!

TO DO    log this to Notes ToDo_log when it's done and tested

-0. multiply selected cells by weight and compare with threshold - done, test on rig
-1. test motion correction on gpu using caiman code - using onacid template matching, working well
-2. disable adding components - done
3. load and save out fixed targets centroid image - done
4. temporally save a sequence of phase masks in holoblink, display on slm on command
5. load threshold and weights file - done
5.5 check trialon trigger - current stim trigger is for sensory stim directly - added 'offsetFrames', done
6. configure stim types from pybehavior output file - done
7. seed cell detection by cell masks from previous recordings
8. set different targets for stim types
9. photoexcitability check (stim one by one after detection) - done, test on rig
test timing on rig
log data for offline anaysis - keep checking if everything is saved out

'''
#%% imports
import sys
import random
import os
import glob
import math
import scipy.io as sio

# Qt
import GUI
from PyQt5.QtCore import Qt,QObject, pyqtSignal, QThread, QTimer, QRectF, QUrl,QPoint, QRect, QSize,QPointF,QSizeF, QSettings, QCoreApplication, QDir, QEventLoop
from PyQt5.QtWidgets import (QComboBox, QCheckBox, QLineEdit, QSpinBox, QLabel,
							 QDoubleSpinBox, QFileDialog, QApplication, QWidget,
							 QDesktopWidget, QMainWindow, QMessageBox, QTableWidgetItem,QAbstractScrollArea)
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
from skimage.morphology import remove_small_objects
import pickle
import logging
import json
from itertools import compress
from pkl2mat import pkl2mat

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
from caiman.base.rois import extract_binary_masks_from_structural_channel as extract_mask
from caiman.utils.visualization import get_contours
from caiman.utils.utils import load_object, save_object
from caiman.paths import caiman_datadir
from caiman.components_evaluation import evaluate_components_CNN
from caiman.source_extraction.cnmf.online_cnmf import RingBuffer
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
from loadMatFile import get_power_params, get_triggertargets_params, get_stimOrder

# sta
from utils import STATraceMaker

#socket
import socket

#%% configure logging
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
	PHOTO_REPLAY_TRIAL = 5 # stimulate all responsive cells at once
	PHOTO_REPLAY_TRAIN = 6 # stim in a sequence
	PHOTO_CLAMP_DOWN = 7 # photostim after sensory stim until cell goes back to baseline (opto inhibition)
	PHOTO_WEIGH_SUM = 8 # compare weighted average of ROI values with tresh
	PHOTO_SEQUENCE = 9 # stimulate loaded targets one by one for photoexcitability test
	PHOTO_FIX_SEQUENCE = 10 # stimulate loaded sequence of targets

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

#%% prepare streaming data from PV
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
						result[result_idx] = 0;
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

#%% start imaging loop
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
					qbuffer.put(thisFrame.copy().astype(np.float32))
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

#% offline, read frames from tiff
		elif p['FLAG_OFFLINE']:
			# load movie
			if p['moviePath'] != 'U:/simulate_movie/20170823.tif':
				print('Loading video')
				Y_ = cm.load(p['moviePath'])
				print(Y_.shape)
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

#%% photostim thread
class Photostimer(QObject):
	# send coordinates to blink and trigger photostimulation
	status_signal            = pyqtSignal(str,name ='statusSignal')
	finished_signal          = pyqtSignal()
	totalNumPhotostim_signal = pyqtSignal(int,name = 'totalNumPhotostimSignal')

	def __init__(self, **kwargs ):
		super(Photostimer,self).__init__()
		self.bl = kwargs.get('blink',[])
		self.NumStimSent = 0
		self.photoPowerPerCell = p['photoPowerPerCell']
		self.power_polyfit_p = p['power_polyfit_p']
		p['FLAG_SKIP_FRAMES'] = False

		if self.bl.send_duration(p['photoDuration']):
			print('duration setting error')

		# setup a trigger for measure timing
#		self.niTimingTask = daqtask.Task() # TTL trigger for measure timing
#		self.niTimingTask.ao_channels.add_ao_voltage_chan('Dev5/ao7','timing_trig')
#		self.niTimingWriter= stream_writers.AnalogSingleChannelWriter(self.niTimingTask.out_stream,True)
#		self.niTimingWriter.write_one_sample(0,10.0)
#		print('timing trigger set')

	def photostimNewTargets(self):
		print('this is photostimer!')
		try:
			if (not qtarget.empty()):
				[thisX,thisY] = qtarget.get(timeout = 10)
				self.sendCoordsAndTriggerStim(thisX,thisY,IF_PHOTOSTIM = True)
				print('coords and stim sent')
				self.NumStimSent +=1
				self.totalNumPhotostim_signal.emit(self.NumStimSent)
				qtarget.task_done()
		except Exception as e:
			print('photostim error')
			print(e)

	def photostimCurrentTargets(self,num_stim_targets):
		ERROR = False
		p['FLAG_SKIP_FRAMES'] = True
		this_volt = np.polyval(self.power_polyfit_p,self.photoPowerPerCell*num_stim_targets)
		if not self.bl.send_trigger_power(this_volt):
			print('Photostimer msg: photostim sent from blink')
		else:
			print('photostimCurrentTargets: msg to blink ERROR!')
			ERROR = True
		p['FLAG_SKIP_FRAMES'] = False
		return ERROR

	def updateNewTargets(self):
		ERROR = False
		if (not qtarget.empty()):
			[thisX,thisY] = qtarget.get(timeout = 10)
			ERROR = self.sendCoordsAndTriggerStim(thisX,thisY,IF_PHOTOSTIM = False)
			qtarget.task_done()
			self.status_signal.emit('Photostimer msg: updated phasemask')
			return ERROR
		else:
			return ERROR

	def sendCoordsAndTriggerStim(self,thisX,thisY,IF_PHOTOSTIM = False):
		ERROR = False
		t0 = time.time()
		p['FLAG_SKIP_FRAMES'] = True
		this_volt = np.polyval(self.power_polyfit_p,self.photoPowerPerCell*len(thisX))
		print('sending coords')
			# scale targets coordinates (to refZoom with which the SLM transform matrix was computed);
		currentTargetX = [int((item-255.0)*p['targetScaleFactor']+255.0) for item in thisX]
		currentTargetY = [int((item-255.0)*p['targetScaleFactor']+255.0) for item in thisY]

		if IF_PHOTOSTIM:
			if self.bl.send_coords_power(currentTargetX, currentTargetY,this_volt):
				print('sendCoordsAndTriggerStim: msg to blink ERROR!')
				ERROR = True
		else:
			if self.bl.send_coords(currentTargetX, currentTargetY):
				print('sendCoordsAndTriggerStim: msg to blink ERROR!')
				ERROR = True

#		self.niTimingWriter.write_many_sample( p['NI_1D_ARRAY'],10.0) # timing test - this trigger didnt work,, dont know why
		print('echo received')
		print(time.time() - t0) # almost same delay as TimerWoeker to BlinkSpiral when no photostim sent - tcp is not the major delay?
		p['FLAG_SKIP_FRAMES'] = False
		return ERROR


	def updatePhotostimParameters(self):
		ERROR = False
		self.photoPowerPerCell = p['photoPowerPerCell']
		if self.bl.send_duration(p['photoDuration']):
			print('duration setting error')
			ERROR = True
		else:
			print('photostim parameters updated')
		return ERROR

#%% processing thread
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
	sendTrialStartTrigger_signal = pyqtSignal()
	sendCoords_signal        = pyqtSignal(name = 'sendCoordsSignal')
	sendPhotoStimTrig_signal = pyqtSignal()
	getROImask_signal        = pyqtSignal(object, object,name = 'getROImaskSignal')
	transDataToMain_signal   = pyqtSignal(object,object,object,object,object,object,object,name = 'transDataToMain')
	updateTargetROIs_signal  = pyqtSignal()
	photostimNewTargets_signal = pyqtSignal()
	updateNewTargets_signal  = pyqtSignal()
	photostimCurrentTargets_signal = pyqtSignal(int)
	finished_signal          = pyqtSignal()


	def __init__(self, **kwargs ):
		super(Worker,self).__init__()

		# get parameters from Main thread
		try:
			self.c = kwargs.get('caiman',{})
			self.pl = kwargs.get('prairie',[])
			self.ROIlist = kwargs.get('ROIlist',dict())
			print('worker pl empty? ' + str(self.pl == []))

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

			# initialise stim frames
			frame_rate = 30 # hard-coded
			self.photo_stim_frames = p['photo_stim_frames']
			self.stim_frames = p['sen_stim_frames'] # trialOn triggers
			self.flag_sta_recording = False
			self.sta_frame_idx = 0
			self.sta_trace_length = p['staPreFrame']+p['staPostFrame']
			self.photo_duration= math.ceil(p['photoDuration']*0.001*frame_rate) # in frames;  frame rate 30hz

			# replay sequence
			self.replay_idx = []
			self.replay_frames = []

			# stim sequence
			self.sequence_idx = p['photo_sequence_idx']

			# count
			self.tot_num_photostim = len(self.photo_stim_frames) # only used for replay or fixed-photostim protocols
			self.tot_num_senstim = len(self.stim_frames)

			print('sensory stim frames:')
			print(self.stim_frames)
			print('num photostim (start) frames:'+str(self.tot_num_photostim))

			# get current parameters
			self.BufferLength = p['BufferLength']
			self.RoiBuffer = p['RoiBuffer']
			self.BufferPointer = p['BufferPointer']
			self.ROIlist_threshold = np.zeros(p['MaxNumROIs'])
	#
	#		self.T1 = self.c['T1']
			self.online_C = np.zeros_like(self.c['cnm2'].C_on)

	#        self.dist_thresh = 21 # reject new rois too close to existing rois
			self.display_shift = False  # option to display motion correction shift
			self.tottime = []

			# lateral motion
			self.shifts = []

			# make Temp folder
			self.save_dir = p['temp_save_dir']

	#		# setup a trigger for measure timing - using this as trial start trigger
	#		self.niTimingTask = daqtask.Task() # TTL trigger for measure timing
	#		self.niTimingTask.ao_channels.add_ao_voltage_chan('Dev5/ao7','timing_trig')
	#		self.niTimingWriter= stream_writers.AnalogSingleChannelWriter(self.niTimingTask.out_stream,True)
	#		self.niTimingWriter.write_one_sample(0,10.0)
	#		print('timing trigger set')

			global STOP_JOB_FLAG
			STOP_JOB_FLAG = False

			p['FLAG_END_LOADING'] = False
			print('worker initialised')
		except Exception as e:
			print('Worker init error:' + e)

	def sortTargets(self,target_idx,target_frames):
		# generate sequence of phase masks
		# bin by photostim duration
		self.replay_idx = []
		self.replay_frames = []
		bin_size = self.photo_duration

		for idx in range(0,len(target_frames),bin_size):
			self.replay_idx.append(np.unique(target_idx[idx:idx+bin_size-1]))
			self.replay_frames.append(target_frames[idx])
		print('replay idx and frames:')
		print(self.replay_idx)
		print(self.replay_frames)

	def work(self):
		print('this is work')
		try:
			# local counts and copy of params
			blockPostFrames = p['blockPostFrames']
			last_com_count = 0
			sens_stim_idx = 0
			LastPlot = 0
			framesProc = self.framesProc # all online frames
			frames_skipped = []
			framesCaiman = 0 # frames count passed to caiman
			photo_stim_frames_caiman = []
			stim_frames_caiman = []
			refreshFrame = p['refreshFrame']
			num_photostim = 0
			current_target_idx = []
			current_target_frames = []
			replay_count = 0
			photostim_train_count = 0
			num_stim_targets = len(p['currentTargetX'])
			last_target_idx = []
			frames_post_photostim = 0
			frames_post_stim = 0
			photo_duration = self.photo_duration
			monitor_frames = p['staPostFrame']+p['offsetFrames'] # offsetFramaes are frame to start monitoring wrt stim start trigger frame (use it when TTL are trialOn triggers for pybehavior rather than direct sensory stim triggers )
			stim_duration = p['staPostFrame']
			wait_frames = p['photoWaitFrames']
			baseline_frames = p['staPreFrame']
			ROIsumThresh = p['ROIsumThresh']
			photostim_flag = 0
			tot_num_photostim = self.tot_num_photostim
			tot_num_senstim = self.tot_num_senstim
			print('duration per photostim = ' + str(photo_duration)+' frames')
			print('tot num photostim = ' + str(tot_num_photostim))
			buflength = self.BufferLength

			FLAG_USE_ONACID = self.UseONACID
			FLAG_SAVE_TIFF = p['saveAsTiff']
			FLAG_SAVE_PROC_TIFF = p['saveProcTiff']
			BufferPointer = self.BufferPointer
			photo_stim_frames = self.photo_stim_frames # predefined fixed frames
			print('photostim frames:')
			print(photo_stim_frames)
			stim_frames = self.stim_frames # stim at fixed intervals

			# flags for photostim and slm signals
			FLAG_TRIG_PHOTOSTIM = False
			FLAG_SEND_COORDS = False

			# control trials
			numExtraPhoto = len(p['extraPhotoFrames'])
			numExtraSensory = len(p['extraSensoryFrames'])
			numExtraSensoryPhoto = len(p['extraSensoryPhotoFrames'])
			extraPhotoFrames = p['extraPhotoFrames']
			extraSensoryFrames = p['extraSensoryFrames']
			extraSensoryPhotoFrames = p['extraSensoryPhotoFrames']
			countExtraPhoto = 0
			countExtraSensory = 0
			countExtraSensoryPhoto = 0
			extraSensoryPhoto_caiman = []
			extraSensory_caiman = []
			extraPhoto_caiman = []
			FLAG_TRIAL_START = False

			# Local buffer for recording online protocol
			online_photo_frames = []
			online_photo_targets = []
			online_photo_x = []
			online_photo_y = []
			online_thresh = []
			online_clamp_level = []
			print('local params in work set')

		except Exception as e:
			print('work local params setting error:')
			print(e)

		# get power control prameters if provided
		try:
			power_polyfit_p = p['power_polyfit_p']
			photoPowerPerCell = p['photoPowerPerCell']
			NI_UNIT_POWER_ARRAY = p['NI_UNIT_POWER_ARRAY']
		except:
			print('power setting error')
			pass

		if FLAG_USE_ONACID:
			# Use locally scoped variables to speed up
			self.status_signal.emit('Using OnACID')

			# Extract initialisation parameters
			cnm2 = self.c['cnm2'] #deepcopy(cnm_init)
			print('cnm2 t = ')
			print(self.c['cnm2'].t)

			img_norm = self.c['img_norm'].copy().astype(np.float32)
			img_min = self.c['img_min'].copy().astype(np.float32)
			ds_factor = self.c['ds_factor']
			dims = self.c['cnm_init'].dims
			gnb = cnm2.gnb
			keep_prev = self.c['keep_prev']
			cell_radius = cnm2.gSig[0]     # dowsampled radius value
			dist_thresh = 2*cell_radius    # reject new rois too close to existing rois - based on cell radius
			reject_mask = self.c['reject_mask']

			# Define OnACID parameters
			max_shift = np.ceil(50./ds_factor).astype('int')  # max shift allowed
			accepted = [ix+gnb for ix in cnm2.accepted] # count from 1 for easy C_on access
			print('worker accepted index:')
			print(accepted)
			rejected_count = cnm2.N - len(accepted)
			expect_components = True
			cnm2.update_num_comps = True
			com_count = len(accepted) #self.c['cnm2'].N
			if p['addNewROIs']==False or com_count == p['MaxNumROIs']:
				expect_components = False
				cnm2.update_num_comps = False
			init_t = cnm2.initbatch
			K_init = self.c['K_init']


			# Local record of all roi coords
			print('number of initial rois '+str(len(accepted)))
			ROIx = np.asarray([item.get("ROIx") for item in self.ROIlist[:com_count]])   # could also predefine ROIx/y as arrays of bigger size (e.g. 100) and fill them in instad of append
			ROIy = np.asarray([item.get("ROIy") for item in self.ROIlist[:com_count]])
			ROIw = np.asarray([item.get("weight") for item in self.ROIlist[:com_count]])
			# Keep or discard previous online cells if rerunning OnACID
			keep_full_traces = True  # keep full trace history of the previous session (may be necessary for deconv of last minibatch, otherwise errors possible)
			if keep_prev:
				coms_init = self.c['coms']
				if keep_full_traces:
					t_cnm = self.c['cnm2'].t
					print(t_cnm)
#                else:
##                    mbs = cnm2.minibatch_shape
#                    init_shape = cnm2.initbatch
#                    t_cnm = init_shape # + mbs
##                    cnm2.C_on[:cnm2.N+cnm2.gnb,init_shape:t_cnm] = cnm2.C_on[:cnm2.N+cnm2.gnb,-mbs:]
#                    cnm2.C_on[:cnm2.N+cnm2.gnb,t_cnm:] = 0  # clears previous results -- should keep at least the last minibatch if clearing
			else:
				coms_init = self.c['coms_init']
				t_cnm = cnm2.initbatch

			init_t = t_cnm
			coms = coms_init.copy()

			store_all_online = True # if True store online traces of both accepted and rejected cells --> important for full plotting in plot window; keep True
			frame_extra_detected = []

			 # fill the init_t trace frames with init results
			if store_all_online:
				self.online_C[:cnm2.M,:init_t] = cnm2.C_on[:cnm2.M,:init_t]
			else:
				accepted = [0] + accepted # include gnb
				self.online_C[:com_count+gnb,:init_t] = cnm2.C_on[accepted,:init_t]

			# Extract opsin mask info constructed from c1v1 image
			try:
				opsin_mask = self.c['opsin_mask']
				opsin_thresh = self.c['opsin_thresh']
				opsin = self.c['cnm2'].opsin
				overlap = self.c['overlap']

			except Exception as e:
				opsin_mask = np.array([])
				opsin_positive = True  # default for when no opsin mask
				opsin = [True for item in range(com_count)]
				print(e)

		# check new cells for opsin
		if opsin_mask.size:
			check_opsin = True
			offset = 2  # for use coms option
			print('Checking opsin online')
		else:
			check_opsin = False

		# check new cells against reject_mask
		if reject_mask.size:
			check_reject = True
			print('Checking reject mask online')
		else:
			check_reject = False

		# store info on why cell was not included
		repeated_idx = []
		rejected_idx = []
		repeated = 0
		rejected = 0

		# prepare to save frames to tiff file
		if FLAG_SAVE_TIFF or FLAG_SAVE_PROC_TIFF:
			try:
				MyTiffwriter =  tifffile.TiffWriter(p['currentMoviePath'], bigtiff=True)
				print('movie will be saved as '+p['currentMoviePath'])
			except Exception as e:
				print('tiffwriter error:' + e)
				FLAG_SAVE_TIFF = False
				FLAG_SAVE_PROC_TIFF = False


		if p['FLAG_OFFLINE']:
			max_wait = None  # wait till video loaded and buffer starts filling
		else:
			max_wait = 10 # timeout
#%% Worker main Loop
		print('Woker entering loop..')
		# keep processing frames in qbuffer
		while ((not p['FLAG_END_LOADING']) or (not qbuffer.empty())) and not STOP_JOB_FLAG:
			app.processEvents(QEventLoop.ExcludeUserInputEvents)
			# skip frames contaminated by photostim but save them to tiff
			if frames_post_photostim >0:
				if p['FLAG_SKIP_FRAMES'] or (frames_post_photostim>0 and frames_post_photostim<photo_duration):
					frames_post_photostim+=1
					frame_in = qbuffer.get(timeout = max_wait)
					qbuffer.task_done()
					frames_skipped.append(framesProc)
					framesProc += 1
					if FLAG_SAVE_TIFF:
						MyTiffwriter.save(frame_in)
					print('framesProc'+str(framesProc)+'skipped')
					continue
				elif frames_post_photostim == photo_duration:
					frames_post_photostim = 0
			elif frames_post_photostim == -1: # let one more frame in
				frames_post_photostim = 2

			# skip frames during (sensory) stimulation in case photostim is triggered during these frames outside of RTAOI
			if blockPostFrames:
				if frames_post_stim>0 and frames_post_stim<stim_duration:
					frames_post_stim +=1
					frame_in = qbuffer.get(timeout = max_wait)
					qbuffer.task_done()
					frames_skipped.append(framesProc)
					framesProc += 1
					if FLAG_SAVE_TIFF:
						MyTiffwriter.save(frame_in)
					print('framesProc'+str(framesProc)+'skipped')
					continue
				elif frames_post_stim == stim_duration:
					frames_post_stim = 0

			# get data from queue
			framesCaiman+=1
			try:
				frame_in = qbuffer.get(timeout = max_wait)
				if FLAG_SAVE_TIFF:
					MyTiffwriter.save(frame_in)
			except Exception: # timeout exception
				print('Timeout Exception: empty qbuffer')
				break # this may be a problem for online too? stops reading the stream completely if timeout happens. continue?

			# timer for processing time
			t0 = time.time()
			framesProc = framesProc+1

			if FLAG_USE_ONACID:
				try:
					# move trace buffer pointer
					if BufferPointer==self.BufferLength-1:
						BufferPointer = 0
					else:
						BufferPointer +=1

					# process current frame
					if ds_factor > 1:
						frame_in = cv2.resize(frame_in, img_norm.shape[::-1])   # downsampling
					frame_in -= img_min                                       # make data non-negative

					# motion correction: template matching with denoised frame
#                   mot_corr_start = time.time()
					templ = cnm2.Ab.dot(cnm2.C_on[:cnm2.M, t_cnm - 1]).reshape(cnm2.dims, order='F') * img_norm
					frame_cor, shift = motion_correct_iteration_fast(frame_in, templ, max_shift, max_shift)
					self.shifts.append(shift)
	#               print('caiman motion correction time:' + str("%.4f"%(time.time()-mot_corr_start)))

					# save to tiff (normalised but not motion corrected)
					if FLAG_SAVE_PROC_TIFF:
						MyTiffwriter.save(frame_cor)

					frame_cor = frame_cor / img_norm                            # normalize data-frame
					cnm2.fit_next(t_cnm, frame_cor.reshape(-1, order='F'))      # run OnACID on this frame

#%% detect new compunents
					if expect_components:
						update_comp_time = time.time()
						if cnm2.N - (com_count+rejected_count) == 1:
							frame_extra_detected.append(t_cnm)
							new_coms = com(cnm2.Ab[:, -1], dims[0], dims[1])[0]

							# Check for repeated components
							close = abs(coms - new_coms) < dist_thresh
							repeated = any(np.all(close,axis=1))

							if check_reject:  # TODO: could also check it by looking at multiple cells around com (like opsin check)
								x = round(new_coms[0])  # no need to inverse here
								y = round(new_coms[1])

								rejected = reject_mask[int(x)][int(y)]


							if not (repeated or rejected): #repeated == False: # add component to ROI
								coms = np.vstack((coms, new_coms))
								y, x = new_coms   # reversed
								ROIx = np.append(ROIx,x*ds_factor)  # ~0.11 ms // filling in empty array: ~0.07 ms
								ROIy = np.append(ROIy,y*ds_factor)

								com_count += 1
								accepted.append(cnm2.N)
								print('New cell detected (' + str(cnm2.N-rejected_count) + ')')

								# Check cell for c1v1
	#                            tt = time_()
								if check_opsin:
									if p['use_mask']:
										cell_A = np.array(cnm2.Ab[:,-1].todense())
										cell_mask = (np.reshape(cell_A, dims, order='F') > 0).astype('int')
										cell_pix = sum(sum(cell_mask == 1))

										inter = cv2.bitwise_and(opsin_mask, cell_mask)
										inter_pix = sum(sum(inter))
										cell_overlap = inter_pix/cell_pix

										overlap.append(cell_overlap)
										opsin_positive = cell_overlap > opsin_thresh
										opsin.append(opsin_positive)

									else:                                    # TODO: could optimise the below
										x_ = int(round(y))
										y_ = int(round(x))
										x_all = [x_, x_+offset, x_-offset]
										y_all = [y_, y_+offset, y_-offset]

										opsin_count = 0
										count = 0
										for i in range(len(x_all)):
											for j in range(len(y_all)):
												op = opsin_mask[x_all[i]][y_all[j]]
												opsin_count += op
												count += 1

										opsin_positive = opsin_count/count >= opsin_thresh
										opsin.append(opsin_positive)
								else:
									opsin.append(True)

								# add new ROI to photostim target, if required
								if p['FLAG_BLINK_CONNECTED'] and p['addNewROIsToTarget']:
									if opsin_positive:  # add target only if opsin present
										p['currentTargetX'].append(x*ds_factor)
										p['currentTargetY'].append(y*ds_factor)
										current_target_idx.append(cnm2.N)
										num_stim_targets = len(p['currentTargetX'])

										if p['stimFromBlink']:
											qtarget.put([p['currentTargetX'].copy(),p['currentTargetY'].copy()])
											self.updateNewTargets_signal.emit()
											# ---  test timing --- DELETE THIS LATER
#												self.niTimingWriter.write_many_sample( p['NI_1D_ARRAY'],10.0)
#												self.photostimNewTargets_signal.emit()
										else:
											self.sendCoords_signal.emit()


										self.updateTargetROIs_signal.emit()
										FLAG_SEND_COORDS = False

								self.getROImask_signal.emit(x,y) # add roi coords to list in main
								print('add new component time:' + str("%.4f"%(time.time()-update_comp_time)))

							else:
								if repeated:
									repeated_idx.append(cnm2.N)
									print('Repeated component found!')
								if rejected:
									rejected_idx.append(cnm2.N)
									print('Rejected component found!')

								rejected_count += 1
						if com_count == p['MaxNumROIs']:
							expect_components = False
							cnm2.update_num_comps = False


					# add data to buffer
					try:
						self.RoiBuffer[:com_count, BufferPointer] = cnm2.C_on[accepted, t_cnm]
#						self.ROIlist_threshold[:com_count] = np.nanmean(self.RoiBuffer[:com_count,:], axis=1) + 2*np.nanstd(self.RoiBuffer[:com_count,:], axis=1)  # changed 3 to 2
						# use noisy C to estimate noise:
						if not (p['photoProtoInx'] == CONSTANTS.PHOTO_CLAMP_DOWN):
							self.ROIlist_threshold[:com_count] = np.nanmean(cnm2.C_on[accepted, t_cnm - buflength:t_cnm], axis=1) + 2*np.nanstd(cnm2.noisyC[accepted, t_cnm - buflength:t_cnm], axis=1)  # changed 3 to 2
							online_thresh.append([items for items in self.ROIlist_threshold[0:com_count]])

						 # record the buffer values for offline analysis
						if store_all_online:
							self.online_C[:cnm2.M, t_cnm] = cnm2.C_on[:cnm2.M, t_cnm]  # storing also background signal (gnb)
						else:
							self.online_C[:com_count, t_cnm] = cnm2.C_on[accepted, t_cnm]
					except Exception as e:
						print('add data to buffer error')
						print(e)
						logger.exception(e)

#%% photostim protocols:
					if p['photoProtoInx'] == CONSTANTS.PHOTO_WEIGH_SUM:
						# for closed-loop trials
						if sens_stim_idx < tot_num_senstim and framesProc == stim_frames[sens_stim_idx]-1:
							# get baseline level
							current_bs_level= np.nanmean(cnm2.C_on[accepted, t_cnm - baseline_frames:t_cnm], axis=1)

						if sens_stim_idx>0 and framesProc < stim_frames[sens_stim_idx-1]+ monitor_frames and framesProc > stim_frames[sens_stim_idx-1]+ wait_frames:
							if p['ROIsumAbove']:
								photostim_flag = np.sum(np.multiply(self.RoiBuffer[:com_count, BufferPointer],ROIw))-current_bs_level- ROIsumThresh
							elif p['ROIsumBelow']:
								photostim_flag = ROIsumThresh - (np.sum(np.multiply(self.RoiBuffer[:com_count, BufferPointer],ROIw))-current_bs_level)
							if photostim_flag>0:
								num_stim_targets = len(p['fixedTargetX'])
								FLAG_TRIG_PHOTOSTIM = True

					elif p['photoProtoInx'] == CONSTANTS.PHOTO_FIX_SEQUENCE:
						if num_photostim < tot_num_photostim:
							if framesProc == photo_stim_frames[num_photostim]:
								FLAG_TRIG_PHOTOSTIM = True

					elif p['photoProtoInx'] == CONSTANTS.PHOTO_CLAMP_DOWN:
						# for closed-loop trials
						if sens_stim_idx < tot_num_senstim and framesProc == stim_frames[sens_stim_idx]-1:
							# get baseline level
							current_clamp_level= np.nanmean(cnm2.C_on[accepted, t_cnm - baseline_frames:t_cnm], axis=1) + np.nanstd(cnm2.noisyC[accepted, t_cnm - baseline_frames:t_cnm], axis=1)
							online_clamp_level.append(current_clamp_level.copy())
							self.ROIlist_threshold[:com_count] = current_clamp_level.copy()
							last_com_count = com_count
							print('updated clamp level')

						if sens_stim_idx>0 and framesProc < stim_frames[sens_stim_idx-1]+ monitor_frames and framesProc > stim_frames[sens_stim_idx-1]+ wait_frames:
							# stim rois: 1.detected before this sensory stim and passed clamp level 2 detected during this sensory stim
							photostim_flag = self.RoiBuffer[:last_com_count, BufferPointer]-current_clamp_level
							photostim_flag = np.concatenate((photostim_flag,np.ones(com_count - last_com_count)))
							above_thresh = np.array(photostim_flag>0)
							opsin_ok = np.array(opsin)
							current_target_idx = np.where(above_thresh & opsin_ok == True)
							num_stim_targets = len(current_target_idx[0])
							if (num_stim_targets>0):
								FLAG_TRIG_PHOTOSTIM = True
								print('photostim on:')
								print(current_target_idx)
								if not np.array_equal(last_target_idx,current_target_idx):
									p['currentTargetX'] = list(ROIx[current_target_idx])
									p['currentTargetY'] = list(ROIy[current_target_idx])
									FLAG_SEND_COORDS = True
									print('different target found:')
									print(last_target_idx)
									print(current_target_idx)


						elif countExtraSensoryPhoto>0 and framesProc < extraSensoryPhotoFrames[countExtraSensoryPhoto-1]+ monitor_frames and framesProc > extraSensoryPhotoFrames[countExtraSensoryPhoto-1]+ wait_frames:
							num_stim_targets = len(p['fixedTargetX'])
							FLAG_TRIG_PHOTOSTIM = True

						elif countExtraPhoto>0 and framesProc < extraPhotoFrames[countExtraPhoto-1]+ monitor_frames and framesProc > extraPhotoFrames[countExtraPhoto-1]+ wait_frames:
							num_stim_targets = len(p['fixedTargetX'])
							FLAG_TRIG_PHOTOSTIM = True
					elif p['photoProtoInx'] == CONSTANTS.PHOTO_SEQUENCE:
						if num_photostim < tot_num_photostim:
							if framesProc == photo_stim_frames[num_photostim]:
								current_target_idx = self.sequence_idx[num_photostim]
								p['currentTargetX'] = [p['fixedTargetX'][current_target_idx]]
								p['currentTargetY'] = [p['fixedTargetY'][current_target_idx]]
								FLAG_TRIG_PHOTOSTIM = True
								FLAG_SEND_COORDS = True
								current_target_idx = [current_target_idx] # make it a list
								print('photostim triggered')

					elif p['photoProtoInx'] == CONSTANTS.PHOTO_REPLAY_TRIAL and num_photostim < tot_num_photostim:
						if (sens_stim_idx>0 and framesProc < stim_frames[sens_stim_idx-1]+ monitor_frames and framesProc > stim_frames[sens_stim_idx-1]):
							# add responsive cells to target group
							photostim_flag = self.RoiBuffer[:com_count, BufferPointer]-self.ROIlist_threshold[:com_count]
							above_thresh = np.array(photostim_flag>0)
							opsin_ok = np.array(opsin)
							current_target_idx = np.append(np.array(current_target_idx),np.where(above_thresh & opsin_ok == True))

						elif framesProc == stim_frames[sens_stim_idx-1] + monitor_frames:
							# update phase mask
							current_target_idx = np.unique(current_target_idx)
							current_target_idx = current_target_idx.astype(int)
							p['currentTargetX'] = list(ROIx[current_target_idx])
							p['currentTargetY'] = list(ROIy[current_target_idx])

							if p['stimFromBlink']:
								qtarget.put([p['currentTargetX'].copy(),p['currentTargetY'].copy()])
								self.updateNewTargets_signal.emit()
							else:
								self.sendCoords_signal.emit()

							self.updateTargetROIs_signal.emit()
							num_stim_targets = len(current_target_idx)

						if framesProc == photo_stim_frames[num_photostim]:
							FLAG_TRIG_PHOTOSTIM = True

					elif  p['photoProtoInx'] == CONSTANTS.PHOTO_REPLAY_TRAIN:
						if (sens_stim_idx>0 and framesProc < stim_frames[sens_stim_idx-1]+ monitor_frames and framesProc > stim_frames[sens_stim_idx-1]):
							# take notes of targets and frames
							photostim_flag = self.RoiBuffer[:com_count, BufferPointer]-self.ROIlist_threshold[:com_count]
							above_thresh = np.array(photostim_flag>0)
							opsin_ok = np.array(opsin)
							current_target_idx = np.append(np.array(current_target_idx),np.where(above_thresh & opsin_ok == True))
							current_target_frames.append(framesProc - stim_frames[sens_stim_idx-1])

						elif sens_stim_idx>0 and framesProc == stim_frames[sens_stim_idx-1] + monitor_frames:
							# prepare for stim sequence
							self.sortTargets(current_target_idx, current_target_frames)
							replay_count = 0
							photostim_train_count +=1

						elif framesProc <= photo_stim_frames[photostim_train_count-1]+ monitor_frames and framesProc >= photo_stim_frames[photostim_train_count-1]:
							# start replay previous trial
							if framesProc == (photo_stim_frames[photostim_train_count-1] + self.replay_frames[replay_count]):
								replay_count+=1
								current_target_idx = self.replay_idx[replay_count].astype(int)
								p['currentTargetX'] = list(ROIx[current_target_idx])
								p['currentTargetY'] = list(ROIy[current_target_idx])
								num_stim_targets = len(current_target_idx)
								if (num_stim_targets>0):
									FLAG_TRIG_PHOTOSTIM = True
									FLAG_SEND_COORDS = True
									print('photostim triggered')

					elif p['photoProtoInx'] == CONSTANTS.PHOTO_FIX_FRAMES:
						if num_photostim < tot_num_photostim:
							if framesProc == photo_stim_frames[num_photostim]:
								FLAG_TRIG_PHOTOSTIM = True

						if countExtraSensoryPhoto < numExtraSensoryPhoto:
							if framesProc == extraSensoryPhotoFrames[countExtraSensoryPhoto]:
								FLAG_TRIAL_START = True
								FLAG_TRIG_PHOTOSTIM = True
						if countExtraPhoto < numExtraPhoto:
							if framesProc == extraPhotoFrames[countExtraPhoto]:
								p['currentTargetX'] = p['fixedTargetX'].copy()
								p['currentTargetY'] = p['fixedTargetY'].copy()
								num_stim_targets = len(p['fixedTargetX'])
								current_target_idx = [] # reset target
								FLAG_TRIAL_START = True
								FLAG_SEND_COORDS = True
								FLAG_TRIG_PHOTOSTIM = True

					elif p['photoProtoInx'] == CONSTANTS.PHOTO_ABOVE_THRESH: # TODO: test the check for opsin
						photostim_flag = self.RoiBuffer[:com_count, BufferPointer]-self.ROIlist_threshold[:com_count]
						above_thresh = np.array(photostim_flag>0)
						opsin_ok = np.array(opsin)
						print('above thresh = ')
						print(above_thresh)
						current_target_idx = np.where(above_thresh & opsin_ok == True)
						num_stim_targets = len(current_target_idx[0])

						if (num_stim_targets>0):
							FLAG_TRIG_PHOTOSTIM = True
							print('photostim trigger true, current target idx')
							print(current_target_idx)
							if not np.array_equal(last_target_idx,current_target_idx):
								p['currentTargetX'] = list(ROIx[current_target_idx])
								p['currentTargetY'] = list(ROIy[current_target_idx])
								FLAG_SEND_COORDS = True

						print('ROIx of current targets:')
						print(ROIx[current_target_idx])

						if countExtraSensoryPhoto < numExtraSensoryPhoto:
							if framesProc == extraSensoryPhotoFrames[countExtraSensoryPhoto]:
								FLAG_TRIAL_START = True
								FLAG_TRIG_PHOTOSTIM = True
						if countExtraPhoto < numExtraPhoto:
							if framesProc == extraPhotoFrames[countExtraPhoto]:
								p['currentTargetX'] = p['fixedTargetX'].copy()
								p['currentTargetY'] = p['fixedTargetY'].copy()
								num_stim_targets = len(p['fixedTargetX'])
								current_target_idx = [] # reset target
								FLAG_TRIAL_START = True
								FLAG_SEND_COORDS = True
								FLAG_TRIG_PHOTOSTIM = True

					elif p['photoProtoInx'] == CONSTANTS.PHOTO_BELOW_THRESH:
						photostim_flag = self.ROIlist_threshold[:com_count] - self.RoiBuffer[:com_count, BufferPointer]
						above_thresh = np.array(photostim_flag>0)
						opsin_ok = np.array(opsin)
						current_target_idx = np.where(above_thresh & opsin_ok == True)
						num_stim_targets = len(current_target_idx[0])

						if (num_stim_targets>0):
							FLAG_TRIG_PHOTOSTIM = True
							if not np.array_equal(last_target_idx,current_target_idx):
								p['currentTargetX'] = list(ROIx[current_target_idx])
								p['currentTargetY'] = list(ROIy[current_target_idx])
								FLAG_SEND_COORDS = True

						print('ROIx of current targets:')
						print(ROIx[current_target_idx])

						if countExtraSensoryPhoto < numExtraSensoryPhoto:
							if framesProc == extraSensoryPhotoFrames[countExtraSensoryPhoto]:
								FLAG_TRIAL_START = True
								FLAG_TRIG_PHOTOSTIM = True
						if countExtraPhoto < numExtraPhoto:
							if framesProc == extraPhotoFrames[countExtraPhoto]:
								p['currentTargetX'] = p['fixedTargetX'].copy()
								p['currentTargetY'] = p['fixedTargetY'].copy()
								num_stim_targets = len(p['fixedTargetX'])
								current_target_idx = [] # reset target
								FLAG_TRIAL_START = True
								FLAG_SEND_COORDS = True
								FLAG_TRIG_PHOTOSTIM = True

					# Trigger sensory stimulation
					if sens_stim_idx < tot_num_senstim:
						if p['enableStimTrigger'] and framesProc == stim_frames[sens_stim_idx]: # send TTL
							self.sendTTLTrigger_signal.emit()
							if p['photoProtoInx'] == CONSTANTS.PHOTO_FIX_SEQUENCE: # change phasemask at sensory trial onsets
								p['currentTargetX'] = p['fixedTargetSeqX'][sens_stim_idx].copy()
								p['currentTargetY'] = p['fixedTargetSeqY'][sens_stim_idx].copy()
								num_stim_targets = len(p['currentTargetX'])
								if num_stim_targets>0:
									FLAG_SEND_COORDS = True
								else:
									FLAG_SEND_COORDS = False

							sens_stim_idx += 1
							stim_frames_caiman.append(framesCaiman)
							# reset target indices and frames (for replay protocols)
							current_target_idx = [] # reset target
							frames_post_stim = 1
							FLAG_TRIAL_START = True

					if countExtraSensoryPhoto <numExtraSensoryPhoto:
						if p['enableStimTrigger'] and framesProc == extraSensoryPhotoFrames[countExtraSensoryPhoto]: # send TTL
							self.sendTTLTrigger_signal.emit()
							p['currentTargetX'] = p['fixedTargetX'].copy()
							p['currentTargetY'] = p['fixedTargetY'].copy()
							num_stim_targets = len(p['fixedTargetX'])
							current_target_idx = [] # reset target
							FLAG_SEND_COORDS = True
							countExtraSensoryPhoto += 1
							extraSensoryPhoto_caiman.append(framesCaiman)
							frames_post_stim = 1
							FLAG_TRIAL_START = True

					if countExtraSensory <numExtraSensory:
						if p['enableStimTrigger'] and framesProc == extraSensoryFrames[countExtraSensory]: # send TTL
							self.sendTTLTrigger_signal.emit()
							countExtraSensory += 1
							extraSensory_caiman.append(framesCaiman)
							frames_post_stim = 1
							FLAG_TRIAL_START = True

					if countExtraPhoto <numExtraPhoto:
						if framesProc == extraPhotoFrames[countExtraPhoto]: # send TTL
							p['currentTargetX'] = p['fixedTargetX'].copy()
							p['currentTargetY'] = p['fixedTargetY'].copy()
							num_stim_targets = len(p['fixedTargetX'])
							current_target_idx = [] # reset target
							FLAG_SEND_COORDS = True
							countExtraPhoto += 1
							extraPhoto_caiman.append(framesCaiman)
							frames_post_stim = 1
							FLAG_TRIAL_START = True

					# Send out trial start trigger
					if FLAG_TRIAL_START:
						self.sendTrialStartTrigger_signal.emit()
						print('trial start trigger sent')
						FLAG_TRIAL_START = False

					# Trigger photostimulation
					if p['FLAG_SKIP_FRAMES']:
						FLAG_TRIG_PHOTOSTIM = False
						FLAG_SEND_COORDS = False

#					FLAG_SEND_COORDS = False # delete this later!
					if p['photoProtoInx'] == CONSTANTS.PHOTO_FIX_SEQUENCE and num_stim_targets==0:
						FLAG_TRIG_PHOTOSTIM = False

					if FLAG_TRIG_PHOTOSTIM:
						frames_post_photostim = 1
						if p['stimFromBlink']: # send trigger from photostimer
							if FLAG_SEND_COORDS:
								qtarget.put([p['currentTargetX'].copy(),p['currentTargetY'].copy()])
								self.photostimNewTargets_signal.emit()
								frames_post_photostim = -1 # accept next frame during phasemask update

							else:
								self.photostimCurrentTargets_signal.emit(num_stim_targets)

						else:
							# update phase mask
							if FLAG_SEND_COORDS:
								self.sendCoords_signal.emit()
							# trigger spiral
							print('sending photostim trigger')
							p['NI_2D_ARRAY'][1,:] = NI_UNIT_POWER_ARRAY *np.polyval(power_polyfit_p,photoPowerPerCell*num_stim_targets)
							self.sendPhotoStimTrig_signal.emit()

						last_target_idx = np.copy(current_target_idx)
						FLAG_TRIG_PHOTOSTIM = False
						if FLAG_SEND_COORDS:
							self.updateTargetROIs_signal.emit() # update display
							FLAG_SEND_COORDS = False

						# take notes and update counts
						online_photo_frames.append(framesProc)
						online_photo_targets.append(current_target_idx[:])
						online_photo_x.append(p['currentTargetX'][:])
						online_photo_y.append(p['currentTargetY'][:])
						num_photostim +=1
						photo_stim_frames_caiman.append(framesCaiman)
						print('Number photostim,',num_photostim)
					elif FLAG_SEND_COORDS: # update phase mask without photostim
						self.sendCoords_signal.emit()
						FLAG_SEND_COORDS = False


					# Update GUI display
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

				except Exception as e:
					print('this is main loop error')
					print(e)
					exc_type, exc_obj, exc_tb = sys.exc_info()
					fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
					print(exc_type, fname, exc_tb.tb_lineno)

			# else not using onACID
			else:
				temp_time = time.time()
				self.roi_signal.emit(frame_in)
				print('display time '+str("%.4f"%(time.time()-temp_time)))


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

		if FLAG_SAVE_TIFF or FLAG_SAVE_PROC_TIFF:
			MyTiffwriter.close()

		if self.UseONACID:
			try:
				self.BufferPointer = BufferPointer
				self.BufferPointer = 0
				# show indices on viewbox
				self.showROIIdx_signal.emit()

				accepted = [idx-1 for idx in accepted] # count from 0
				cnm2.accepted = accepted # for easier access in MainWindow
				cnm2.t = t_cnm
				self.c['cnm2'].t = t_cnm
				self.c['coms'] = coms
				if opsin_mask.size:
					cnm2.opsin = opsin

				frame_detected_init = np.ones(K_init).astype('int')*cnm2.initbatch
				frame_detected = np.concatenate((frame_detected_init, frame_extra_detected))

				# save results to the results folder
				self.movie_name = p['currentMoviePath']

				save_dict = dict()
				save_dict['cnm2'] = deepcopy(cnm2)  # opsin info a part of cnm struct for now
				save_dict['online_C']  = self.online_C
				save_dict['init_com_count'] = K_init # com_count from init file (in case any cells removed from init file)
				save_dict['online_com_count'] = com_count
				save_dict['online_photo_x'] = online_photo_x
				save_dict['online_photo_y'] = online_photo_y
				save_dict['online_thresh'] = online_thresh
				save_dict['online_clamp_level'] = online_clamp_level

				save_dict['accepted_idx'] = accepted  # accepted currently stored inside cnm2 as well
				save_dict['repeated_idx'] = repeated_idx
				save_dict['rejected_idx'] = rejected_idx

				save_dict['frame_detected'] = frame_detected
				save_dict['t_cnm'] = t_cnm
				save_dict['tottime'] = self.tottime
				save_dict['shifts'] = self.shifts
				save_dict['coms'] = coms
				save_dict['ds_factor'] = ds_factor
				save_dict['t_init'] = init_t

				# for easier offline check
#				save_dict['Cn'] = self.c['Cn_init']
				save_dict['img_norm'] = img_norm
				save_dict['img_min'] = img_min

				save_dict['online_photo_frames'] = online_photo_frames
				save_dict['online_photo_targets'] = online_photo_targets
				save_dict['photo_stim_frames_caiman'] = photo_stim_frames_caiman
				save_dict['stim_frames_caiman'] = stim_frames_caiman
				save_dict['frames_skipped'] = frames_skipped
				save_dict['framesProc'] = framesProc

				save_dict['extraSensoryPhoto_caiman'] = extraSensoryPhoto_caiman
				save_dict['extraPhoto_caiman'] = extraPhoto_caiman
				save_dict['extraSensory_caiman'] = extraSensory_caiman

				save_dict['all_trial_types'] = p['all_trial_types']
				save_dict['all_trial_frames'] = p['all_trial_frames']

				save_dict['keep_prev'] = keep_prev
				save_dict['reject_mask'] = reject_mask
			# keep notes of initialisation file used
			except Exception as e:
				print('this is making save dict error')
				print(e)
				exc_type, exc_obj, exc_tb = sys.exc_info()
				fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
				print(exc_type, fname, exc_tb.tb_lineno)
			try:
				save_dict['coms_init'] = self.c['coms_init']
				save_dict['cnm_init'] = self.c['cnm_init']
				save_dict['coms_init_orig'] = self.c['coms_init_orig']
				save_dict['cnm_init_orig'] = self.c['cnm_init_orig']
			except:
				print('init fields not found')

			# save out mat file for STAMovieMaker
			mat_dict = dict()
			mat_dict['oris'] = p['all_trial_types']
			mat_dict['frames'] = p['all_trial_frames']

			if p['enableStimTrigger'] :
				save_dict['sensory_stim_frames'] = self.stim_frames
			else:
				save_dict['sensory_stim_frames'] = []

			try:
				save_dict['opsin_thresh'] = opsin_thresh
				save_dict['opsin_mask'] = opsin_mask
				save_dict['overlap'] = overlap
				# opsin currently as a property of cnm2
			except: pass

			try:
				save_dict['p'] = p
			except:
				print('p not saved')

			try:
				print('Saving onacid output')
#                save_separately = True # temp flag to save in a new folder inside movie folder - default
				daytimestr = time.strftime("%Y%m%d-%H%M%S")
				timestr = daytimestr[-6:]

				movie_name = os.path.basename(self.movie_name)
				i = movie_name.find('Cycle')
				if i>-1:
					movie_string = movie_name[:i] + 'rtaoi'
				else:
					movie_string = movie_name[:-4]

				filename = movie_string + '_OnlineProc_' + timestr + '.pkl'
				matname =  movie_string + '_OnlineProc_' + timestr + '_TrialType' + '.mat'

				movie_folder = os.path.dirname(self.movie_name)
				save_folder = os.path.join(movie_folder, 'pyrtaoi_results')  # save init result in a separate folder

				if not os.path.exists(save_folder):
					os.makedirs(save_folder)
				save_path = os.path.join(save_folder, filename)
				mat_save_path = os.path.join(save_folder, matname)
#                else:
#                    save_path = self.movie_name[:-4] + '_OnlineProc_' + timestr + '.pkl'   # save results in the same folder
				print('OnACID result is being saved as: ' + save_path)
				save_object(save_dict, save_path)
				print('pkl file saved')

				sio.savemat(mat_save_path,mat_dict)
				print('tiral type fle saved')

				pkl2mat(file_full_name = save_path)
				print('mat file saved')

			except Exception as e:
				print(e)


			# transfer data to main and show traces in plot tab
			self.transDataToMain_signal.emit(cnm2, self.online_C, coms, accepted, t_cnm,[],[])

			# delete big variables
			del self.c
			del ROIx
			del ROIy

#		try:
#			self.niTimingWriter.write_one_sample(0,10.0)
#			self.niTimingTask.stop()
#			self.niTimingTask.close()
#			del self.niTimingWriter
#			print('closed timing trigger task')
#		except Exception as e:
#			print(e)

		# finishing
		self.finished_signal.emit()
		print('worker finished')

#    def photostimTargets(self, coords):

	def stop(self):
		self.status_signal.emit('Stop clicked')
		print('worker stop')
		global STOP_JOB_FLAG
		STOP_JOB_FLAG = True


#%% Main GUI
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

		# trial type
		p['trialOrder'] = []
		p['all_trial_types'] = []
		p['photo_sequence_idx'] = []

		# camain init
		p['FLAG_USING_RESULT_FILE'] = False

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
		self.use_mask = False
		self.opsinMaskOn = False
		self.showROIsOn = False
		self.A_loaded = np.array([])
		self.mask_path = ''
		self.reject_mask = np.array([])

		# ROI selection
#        self.drawOn = False
		self.RoiCenter = QPoint(244,244)
		self.RoiRadius = 10
		self.thisROI = pg.CircleROI(self.RoiCenter,self.RoiRadius*2)
		self.removeModeOn = False
		self.removeIdx = []
		self.rejectModeOn = False
		self.rejectIdx = []
		self.rejectMaskOn = False

		# initialise ROI list
		self.InitNumROIs = 3
		self.MaxNumROIs = 30  # desired number of cells
		self.thisROIIdx = 0  # current number of ROIs
		self.ALL_ROI_SELECTED = False
		self.ROIlist = dict()
		self.resetROIlist()

		# mean intensity of ROI
		self.BufferLength = 100 # make it same as minibatch
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
		self.RejectPen = pg.mkPen(color = (255, 0, 0), width = 3, style = Qt.SolidLine)
		self.TriggerPen = pg.mkPen(color = (137,182,255), width = 3, style = Qt.DotLine)
		self.Targetcontour_item = pg.ScatterPlotItem()
		self.Targetcontour_item.setBrush(self.scatterbrush)
		self.graphicsView.addItem(self.Targetcontour_item)
		self.Triggercontour_item = pg.ScatterPlotItem()
		self.Triggercontour_item.setBrush(self.scatterbrush)
		self.graphicsView.addItem(self.Triggercontour_item)

		# dictionary of ROIs
		self.ROIlist = [dict() for i in range(self.MaxNumROIs)] # this only include accepted rois
		self.Thresh_tableWidget.setColumnCount(5)
		self.Thresh_tableWidget.setHorizontalHeaderLabels(['X','Y','Thresh','W','STA'])
		self.Thresh_tableWidget.setSizeAdjustPolicy(QAbstractScrollArea.AdjustToContents)
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

		# trigger list
		self.TriggerIdx = []
		self.TriggerWeights = []
		self.TriggerFrames = []


		# photostim target list
		self.numTargets = 0
		self.bl = []
		self.TargetIdx = [] # indices in ROIlist
		self.TargetX = [] # all target coords
		self.TargetY = []
		p['ExtraTargetX'] = [] # selected targets (not in the ROIlist)
		p['ExtraTargetY'] = []
		p['currentTargetX'] = [] # keep track of what is currently on SLM
		p['currentTargetY'] = []
		p['FLAG_BLINK_CONNECTED'] = False
		p['stimFromBlink'] = False
		p['targetScaleFactor']  = 1 # currentZoom/refZoom
		p['FLAG_SKIP_FRAMES'] = False
		p['fixedTargetX'] = []
		p['fixedTargetY'] = []
		p['fixedTargetSeqX'] = []
		p['fixedTargetSeqY'] = []

		# reject list
		p['rejectedX'] = []
		p['rejectedY'] = []

		# groupbox list
		self.groupBoxes = [self.Prairie_groupBox,self.caiman_groupBox,
						   self.Blink_groupBox, self.opsinMask_groupBox,
						   self.StimOptions_groupBox, self.OfflineAnalysis_groupBox,
						   self.photostim_groupBox,
						   self.groupBox, self.config_groupBox, self.run_groupBox] # removed self.DisplayOptions_groupBox

		# offline
		self.IsOffline = False

		# initialise FLAGs
		p['enableStimTrigger'] = False
		self.photoProtoInx = CONSTANTS.PHOTO_NONE

		# make dir to save data temporally
		p['temp_save_dir'] = os.getcwd()+'/Temp/'
		if not os.path.exists(p['temp_save_dir']):
			os.makedirs(p['temp_save_dir'])
		p['saveResultPath'] = p['temp_save_dir'] +str('temp.pkl')
		self.saveResultPath_lineEdit.setText(p['saveResultPath'])

		# power control
		self.POWER_CONTROL_READY = False

		# get gui elements
		self.getValues()
		self.updateTrialType()

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
		self.photostimObject = None
		self.niStimWriter = None
		self.niPhotoStimWriter = None
		self.niPhotostimFullWriter = None

		# daq config
		NI_SAMPLE_RATE = 10000
		NI_TTL_NUM_SAMPLES = int(0.01*NI_SAMPLE_RATE) #10 ms
		NI_STIM_NUM_SAMPLES = int(p['photoDuration']*0.001*NI_SAMPLE_RATE)
		p['NI_SAMPLE_RATE'] = NI_SAMPLE_RATE
		p['NI_TTL_NUM_SAMPLES'] = NI_TTL_NUM_SAMPLES
		self.StimTask_READY = False
		self.TrialTask_READY = False

		try:
			self.daq_array = np.ones((NI_TTL_NUM_SAMPLES-1), dtype = np.float)*5 # to write to single output channel
			self.daq_array = np.append(self.daq_array,0)
			self.daq_1s_array = np.ones((int(NI_SAMPLE_RATE)-1), dtype = np.float)*5
			self.daq_1s_array = np.append(self.daq_1s_array,0)
			self.niStimTask = daqtask.Task()
			self.niPhotoStimTask = daqtask.Task() #  TTL trigger only
			self.niPhotoStimFullTask = daqtask.Task() # TTL trigger plus power control
			p['NI_1D_ARRAY'] = self.daq_array
		except:
			print('setting daq tasks error')

		try: # sensory stim
			self.niStimTask.ao_channels.add_ao_voltage_chan(p['stimDaqDevice']) # ao3
			self.niStimWriter = stream_writers.AnalogSingleChannelWriter(self.niStimTask.out_stream,True)
			self.niStimWriter.write_one_sample(0,10.0)
			print('sensory stim channel created')
			self.StimTask_READY = True

		except Exception as e:
			print(e)
			logger.exception(e)
			time.sleep(2)

		try: # trial start trigger
			self.niTrialStartTask = daqtask.Task()
			self.niTrialStartTask.ao_channels.add_ao_voltage_chan('Dev5/ao0','trialStart')
			self.niTrialStartWriter = stream_writers.AnalogSingleChannelWriter(self.niTrialStartTask.out_stream,True)
			self.niTrialStartWriter.write_one_sample(0,10.0)
			self.niTrialStartWriter.write_one_sample(1,10.0)
			time.sleep(1)
			self.niTrialStartWriter.write_one_sample(0,10.0)
			print('trial start channel created')
			self.TrialTask_READY = True

		except Exception as e:
			print(e)
			logger.exception(e)
			time.sleep(2)

		try: # photo stim
			init_output = np.zeros([2, max(NI_TTL_NUM_SAMPLES+1,NI_STIM_NUM_SAMPLES+1)])
			init_output.ravel()[:len(self.daq_array)] = self.daq_array  # use ravel to data between array (fastest way)
			init_output[1,:-1] = 1
			p['NI_2D_ARRAY'] = np.copy(init_output)
			p['NI_UNIT_POWER_ARRAY'] = np.copy(init_output[1,:])
			p['NI_2D_ARRAY'][1,:] = p['NI_UNIT_POWER_ARRAY'] *np.polyval(p['power_polyfit_p'],p['photoPowerPerCell']) # single cell power

			# single-channel writer - only send triggers, used when power control is disabled
			self.niPhotoStimTask.ao_channels.add_ao_voltage_chan('Dev5/ao2','photostim_simple_trig') # hard coded
			self.niPhotoStimWriter= stream_writers.AnalogSingleChannelWriter(self.niPhotoStimTask.out_stream,True)
			self.niPhotoStimWriter.write_one_sample(0,10.0)
			print('single channel photostim trigger set')

			# multi-channel writer - send triggers and voltage to power control device  - worked when using ao0 and ao1
			self.niPhotoStimFullTask.ao_channels.add_ao_voltage_chan(p['photostimDaqDevice'],'photostim_trig')
			self.niPhotoStimFullTask.ao_channels.add_ao_voltage_chan(p['powerDaqDevice'],'photostim_power')
			self.niPhotostimFullWriter = stream_writers.AnalogMultiChannelWriter(self.niPhotoStimFullTask.out_stream,auto_start = True)
			self.niPhotoStimFullTask.timing.cfg_samp_clk_timing(rate = NI_SAMPLE_RATE,sample_mode= nidaqmx.constants.AcquisitionType.FINITE, samps_per_chan = max(NI_TTL_NUM_SAMPLES+1,NI_STIM_NUM_SAMPLES+1))

			# send a trigger - for debug
#			self.niPhotostimFullWriter.write_many_sample(p['NI_2D_ARRAY'])
#			while(not self.niPhotoStimFullTask.is_task_done()):
#				pass
#			self.niPhotoStimFullTask.stop()
#			self.niPhotostimFullWriter = stream_writers.AnalogMultiChannelWriter(self.niPhotoStimFullTask.out_stream,auto_start = True)

			print('multi channel photostim trigger set')
			self.POWER_CONTROL_READY = True

		except Exception as e:
			print('power control initialisation error')
			print(str(e))
			self.updateStatusBar(str(e))


		# sensory TCP config
		self.SENSTIM_IP = '128.40.156.162' #OPTILEX IP
		self.SENSTIM_PORT = 8087
		self.FlAG_stimSOCK_READY = False


		# make threads
		self.workerThread = MyThread()
		self.streamThread = MyThread()
		self.imageSaverThread = MyThread()
		self.photostimThread = MyThread()

		# flag reading/streaming data
		p['FLAG_END_LOADING'] = False

		# allow drag and drop
		self.setAcceptDrops(True)

		# signal/slot connections
		self.setConnects()

		# methods for elements within MainWindow
		self.ImageWindow.enterEvent = self.displayImg.__get__(self.ImageWindow)  # ImageWindow when e.g. opsinMaskOn
		self.ImageWindow.leaveEvent = self.displayMask.__get__(self.ImageWindow)
		self.figCanvas.keyPressEvent = self.arrow_key_image_control.__get__(self.figCanvas) # changing cell display in the plot tab with arrows

	def tempTest(self):
		print('test button clicked')
		try:
			self.sendTrialStartTrigger()
		except Exception as e:
			print(e)
#%% connect GUI elements
	def setConnects(self):
		# load and save
		self.loadMoviePath_pushButton.clicked.connect(self.loadMoviePath)
		self.loadRefMoviePath_pushButton.clicked.connect(self.loadRefMoviePath)
		self.loadOpsinImgPath_pushButton.clicked.connect(self.loadOpsinImg)
		self.loadMask_pushButton.clicked.connect(self.loadMask)
		self.initialise_pushButton.clicked.connect(self.initialiseCaiman)
		self.saveConfig_pushButton.clicked.connect(self.saveConfig)
		self.loadConfig_pushButton.clicked.connect(self.loadConfig)
		self.saveResult_pushButton.clicked.connect(self.saveResults)
		self.browseResultPath_pushButton.clicked.connect(self.browseResultPath)
		self.loadCalciumImg_pushButton.clicked.connect(self.loadCalciumImg)
		self.loadResultPath_pushButton.clicked.connect(self.loadResultPath)
		self.loadTargetCentroid_pushButton.clicked.connect(self.loadTargetCentroid)
		self.loadTriggerConfig_pushButton.clicked.connect(self.loadTriggerConfig)
		self.saveCurrentTargetCentroids_pushButton.clicked.connect(self.saveCurrentTargetCentroids)
#        self.createRejectMask_pushButton.clicked.connect(self.autoCreateRejectMask)

		# rois
		self.createMask_pushButton.clicked.connect(self.createOpsinMask)
		self.showCellsOnMask_pushButton.clicked.connect(self.showCellsOnMask)
		self.removeMode_pushButton.clicked.connect(self.removeModeController) # remove manually selected or CNN filtered ROIs
		self.rejectMode_pushButton.clicked.connect(self.rejectModeController) # reject filled cells in gcamp image
		self.resetRejectMask_pushButton.clicked.connect(self.resetRejectMask)
		self.applyRejectMask_pushButton.clicked.connect(self.applyRejectMask)
		self.rejectMaskDisplay_pushButton.clicked.connect(self.rejectMaskDisplay)
		self.showCellsOnRejectedMask_pushButton.clicked.connect(self.showCellsOnRejectedMask)
		self.CNNFilter_pushButton.clicked.connect(self.filterResults)
		self.removeSelected_pushButton.clicked.connect(self.removeSelected)
		self.keepSelected_pushButton.clicked.connect(self.keepSelected)

		# start worker
		self.run_pushButton.clicked.connect(self.clickRun)

		# start photostimer
		self.stimFromBlink_checkBox.stateChanged.connect(self.outsourcePhotostim)
		self.updatePhotostimParameters_pushButton.clicked.connect(self.updatePhotostimParameters)

		# display
		self.plotOn_checkBox.clicked.connect(self.switch_plotOn)
		self.displayOn_checkBox.clicked.connect(self.switch_displayOn)
		self.denoiseOn_checkBox.clicked.connect(self.switch_denoiseOn)
		self.plotSTA_pushButton.clicked.connect(self.plotSTA)
		self.plotSTAonMasks_pushButton.clicked.connect(lambda: self.plotSTAonMasks(self.sta_amp))
		self.showComponents_pushButton.clicked.connect(lambda: self.plotSTAonMasks(None))
		self.showFOV_pushButton.clicked.connect(self.showFOV)
		self.showROIIdx_pushButton.clicked.connect(self.showIdx)
		self.loadProcData_pushButton.clicked.connect(self.loadProcData)

		# triggers
		self.enableStimTrigger_checkBox.clicked.connect(self.enableStimTrigger)
		self.testTTLTrigger_pushButton.clicked.connect(self.testTTLTrigger)
		self.testPhotoStimTrigger_pushButton.clicked.connect(self.testPhotoStimTrigger)
		self.connectSenStimTCP_pushButton.clicked.connect(self.connectSenStimTCP)
		self.removeButTrigger_pushButton.clicked.connect(self.removeButTrigger)

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
		self.selectAll_checkBox.clicked.connect(self.selectAllROIsClicked)
		self.clearTargets_pushButton.clicked.connect(self.clearTargets)
		self.updateTargets_pushButton.clicked.connect(self.updateTargets)
		self.updateFixTargets_pushButton.clicked.connect(self.updateFixTargets)
		self.loadStimOrder_pushButton.clicked.connect(self.loadStimOrder)
		self.previewStimOrder_pushButton.clicked.connect(self.previewStimOrder)

		# mouse events
		self.graphicsView.scene().sigMouseClicked.connect(self.getMousePosition)
		self.graphicsView.scene().sigMouseMoved.connect(self.showMousePosition)

		# reference movie
		self.takeRefMovie_pushButton.clicked.connect(self.takeRefMovie)

		# calibration check
		self.calcheck_pushButton.clicked.connect(self.calCheck)

		# sequence stimulation config
		self.configSeqStim_pushButton.clicked.connect(self.configSeqStim)

		# set weights (in ROIlist table)
		self.Thresh_tableWidget.itemChanged.connect(self.tableItemChanged)

		# others
		self.test_pushButton.clicked.connect(self.tempTest)
		self.reset_pushButton.clicked.connect(self.resetAll)

		# tiral types
		self.plotTrialTypes_pushButton.clicked.connect(self.plotTrialTypes)
		self.updateTrialType_pushButton.clicked.connect(self.updateTrialType)
		self.loadTrialOrder_pushButton.clicked.connect(self.loadTrialOrder)

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

		self.photoProtoInx = self.photoProto_comboBox.currentIndex()

		p['BufferLength'] = self.BufferLength
		p['RoiBuffer'] = self.RoiBuffer
		p['BufferPointer'] = self.BufferPointer
		p['FLAG_PV_CONNECTED'] = self.FLAG_PV_CONNECTED
		p['FLAG_OFFLINE'] = self.IsOffline
		p['photoProtoInx'] = self.photoProtoInx
		p['use_mask'] = self.use_mask

		self.flipRowChanged()
		self.currentZoomChanged()

		if self.POWER_CONTROL_READY:
			self.updatePowerControl()
		p['tot_num_trials'] = p['numberStims']+p['ifExtraSensory']*p['numExtraSensory']+ p['ifExtraPhoto']*p['numExtraPhoto']+p['numExtraSensoryPhoto']*p['ifExtraSensoryPhoto']

		self.framesNeeded_label.setText('Frames needed:'+str(p['stimStartFrame']+p['interStimInterval']*p['tot_num_trials']+p['staPostFrame']))

		try:
			self.c['cnm2'].max_num_added = p['MaxNumROIs']
		except:
			pass


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
#%% set trial types
	def updateFixTargets(self):
		p['fixedTargetX'] = p['currentTargetX'].copy()
		p['fixedTargetY'] = p['currentTargetY'].copy()
		print('fixed targets updated')

	def plotTrialTypes(self):
		type_names1 = ['extraSensory','extraPhoto','extraSensoryPhoto']
		type_names2 = ['ExtraSensory','ExtraPhoto','ExtraSensoryPhoto']
		pl.figure()

		for typeidx in range(0,len(type_names1)):
			if p[type_names1[typeidx]+'Frames']:
				pl.scatter(p[type_names1[typeidx]+'Frames'],[typeidx+1]*p['num'+type_names2[typeidx]])
#			except Exception as e:
#				print('plot '+type_names1[typeidx]+' error')
#				print(e)
#				continue

		if p['sen_stim_frames']:
			pl.scatter(p['sen_stim_frames'],[0]*p['numberStims'])
		pl.yticks(np.arange(len(type_names1)+1), ('CL','ExtraSensory','ExtraPhoto','ExtraSensoryPhoto'))
		pl.ylim(-1,len(type_names1)+1)
		pl.xlabel('Frames')
		pl.show(block=False)

		print('tot num trials:'+str(p['tot_num_trials']))


	def updateTrialType(self):
		frame_rate = 30
		# these protocols require trains of trials:
		# PHOTO_FIX_FRAMES PHOTO_REPLAY_TRIAL PHOTO_REPLAY_TRAIN PHOTO_CLAMP_DOWN

		p['extraSensoryFrames'] = []
		p['extraSensoryPhotoFrames'] = []
		p['extraPhotoFrames'] = []
		photo_stim_frames = []
		tot_num_trials = p['tot_num_trials']
		all_trial_frames = p['stimStartFrame']+p['interStimInterval']*np.arange(0,tot_num_trials)


		# adding control trials
#		if p['photoProtoInx'] in [CONSTANTS.PHOTO_FIX_FRAMES, CONSTANTS.PHOTO_CLAMP_DOWN,CONSTANTS.PHOTO_REPLAY_TRAIN,CONSTANTS.PHOTO_CLAMP_DOWN]:
		if p['photoProtoInx'] != CONSTANTS.PHOTO_SEQUENCE:
			ifCloseloopTrials = 1
		else:
			ifCloseloopTrials = 0

		if not p['useExternalTrialOrder']:
			all_trial_types = []
			if ifCloseloopTrials:
				all_trial_types+=[1]*p['numberStims']
			if p['ifExtraSensory']:
				all_trial_types+=[2]*p['numExtraSensory']
			if p ['ifExtraPhoto']:
				all_trial_types+=[3]*p['numExtraPhoto']
			if p['ifExtraSensoryPhoto']:
				all_trial_types+=[4]*p['numExtraSensoryPhoto']
			print('generated trial types')
		else:
			if len(p['trialOrder']>0):
				# update trial types according to external file
				p['tot_num_trials'] = len(p['trialOrder'])
				tot_num_trials = p['tot_num_trials']
				all_trial_frames = p['stimStartFrame']+p['interStimInterval']*np.arange(0,tot_num_trials)
				trialOrder_types = np.unique(np.array(p['trialOrder']))

				# enable conditioning for a fraction of trials, for specified stimType
				condition_idx = []
				fraction_trials = p['pcConditionTrial']/100
				for i in range(len(trialOrder_types)):
					if p['ConditionOnStimType'] & (trialOrder_types[i]!= p['senStimType']):
						continue
					this_condition_idx = [ii for ii, x in enumerate(p['trialOrder']) if x == trialOrder_types[i]]
					this_num_idx = len(this_condition_idx)
					this_keep_idx= random.sample(range(this_num_idx),round(this_num_idx*fraction_trials))
					this_condition_idx = [x for ii, x in enumerate(this_condition_idx) if ii in this_keep_idx]
					condition_idx += this_condition_idx

				num_condition_trials = len(condition_idx)
				all_trial_types = [2]*p['tot_num_trials']
				for i in condition_idx:
					all_trial_types[i] = 1

				# update GUI and p
				p['all_trial_types'] = all_trial_types
				p['numberStims'] = num_condition_trials
				self.numberStims_spinBox.setValue(num_condition_trials)
				p['numExtraSensory'] = tot_num_trials - num_condition_trials
				self.numExtraSensory_spinBox.setValue(p['numExtraSensory'])
				self.ifExtraPhoto_checkBox.setCheckState(Qt.Unchecked)
				p['ifExtraPhoto'] = False
				self.ifExtraSensoryPhoto_checkBox.setCheckState(Qt.Unchecked)
				p['ifExtraSensoryPhoto'] = False
				if p['numExtraSensory']>0:
					p['ifExtraSensory'] = True
					self.ifExtraSensory_checkBox.setCheckState(Qt.Checked)
				self.shuffleTrials_checkBox.setCheckState(Qt.Unchecked)
				p['shuffleTrials'] = False

				print('using trial types from trialOrder file')
			else:
				print('load trial type file!')
				return

		if p['shuffleTrials']:
			random.shuffle(all_trial_types)

		sen_stim_frames = [all_trial_frames[i] for i,x in enumerate(all_trial_types) if x == 1] # these will be closed-loop trials for clamping down experiments
		p['extraSensoryFrames'] =  [all_trial_frames[i] for i,x in enumerate(all_trial_types) if x == 2]
		p['extraPhotoFrames'] =  [all_trial_frames[i] for i,x in enumerate(all_trial_types) if x == 3]
		p['extraSensoryPhotoFrames'] =  [all_trial_frames[i] for i,x in enumerate(all_trial_types) if x == 4]


			# FLAG_FIX_PHOTO
		# closed-loop photo and sensory are coupled in these protocols
		if p['photoProtoInx'] == CONSTANTS.PHOTO_REPLAY_TRIAL or p['photoProtoInx'] == CONSTANTS.PHOTO_REPLAY_TRAIN:
			photo_stim_frames = sen_stim_frames+p['photoWaitFrames']
			sen_stim_frames = sen_stim_frames[::2]
		if p['photoProtoInx'] == CONSTANTS.PHOTO_FIX_FRAMES or p['photoProtoInx'] == CONSTANTS.PHOTO_FIX_SEQUENCE:
			relative_frames = math.ceil(frame_rate/p['photoFrequency'])*np.arange(0,p['photoNumStimsPerTrial'])+p['photoWaitFrames']
			print(sen_stim_frames)
			photo_stim_frames = np.concatenate( [item+relative_frames for item in sen_stim_frames])
			print(photo_stim_frames)

		# save to global p
		p['sen_stim_frames'] = sen_stim_frames
		p['photo_stim_frames'] = photo_stim_frames
		p['all_trial_types'] = all_trial_types
		p['all_trial_frames'] = all_trial_frames
#%% threading control functions:
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
			p['currentMoviePath'] = self.pl.get_movie_name()+ '_rtaoi_DS_' + str(self.ds_factor) +'.tif'
			print('current movie:' + p['currentMoviePath'])
		else:
			p['currentMoviePath'] = os.path.splitext(p['moviePath'])[0]+ '_rtaoi_DS_' + str(self.ds_factor) +'.tif'

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
						self.removeIdx = np.arange(self.c['K_init'],curr_count)
						self.removeCells(save_init=False,reinitiate = not p['FLAG_USING_RESULT_FILE'])
					else:
						self.c['keep_prev'] = True
				else:
					return

			self.deleteTextItems()

			if self.opsinMaskOn or self.showROIsOn: # or self.calciumMaskOn:
				self.opsinMaskOn = False
#                self.calciumMaskOn = False
				self.showROIsOn = False
				self.updateImage(cv2.resize(self.c['img_norm'],(512,512),interpolation=cv2.INTER_CUBIC))

			# ensure correct param values following init are displayed
			print('setting parameter values in gui')
			self.dsFactor_doubleSpinBox.setValue(self.ds_factor)
			self.cellRadius_spinBox.setValue(self.cell_radius)
			self.minSNR_doubleSpinBox.setValue(self.min_SNR)
			self.rval_thresh_doubleSpinBox.setValue(self.rval_thr)
			self.merge_thresh_doubleSpinBox.setValue(self.merge_thresh)
			self.overlap_thresh_doubleSpinBox.setValue(self.thresh_overlap)

			print('updating ROI list')
			if self.MaxNumROIs_spinBox.value() > self.MaxNumROIs:
				for extraROI in range(self.MaxNumROIs_spinBox.value() - self.MaxNumROIs):  # add extra fields to ROIlist
					self.ROIlist.append(dict())
					ROIIdx = self.MaxNumROIs + extraROI
					self.ROIlist[ROIIdx]["ROIx"]= list()
					self.ROIlist[ROIIdx]["ROIy"] = list()
					self.ROIlist[ROIIdx]["ROIcolor"] = self.random_color()
					self.ROIlist[ROIIdx]["threshold"] = list()
					self.ROIlist[ROIIdx]["STA"] = list()
					self.ROIlist[ROIIdx]["weight"] = 1

				self.MaxNumROIs = self.MaxNumROIs_spinBox.value()

				self.updateTable()
				self.RoiBuffer = np.zeros([self.MaxNumROIs,self.BufferLength])

			print('updating values')
			self.getValues()

			print('setting uo buttons')
			# connect stop_thread to stop button when analysis running
			self.stop_pushButton.clicked.connect(self.stop_thread)

			def resetworkerObject():
				self.workerObject = None

			def resetstreamObject():
				self.streamObject = None

			print('setting up worker object')
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

			# signaling to photostim thread
			if (not self.IsOffline) and self.photostimObject!=None:
				try:
					self.workerObject.photostimNewTargets_signal.connect(self.photostimObject.photostimNewTargets) # update phase mask and send photostim from blink
					self.workerObject.updateNewTargets_signal.connect(self.photostimObject.updateNewTargets) # update phase mask
					self.workerObject.photostimCurrentTargets_signal.connect(self.photostimObject.photostimCurrentTargets)  # photostim only
				except Exception as e:
					print(e)
					return

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
			self.workerObject.sendTTLTrigger_signal.connect(self.sendTTLTrigger) # TTL and/or TCP trigger
			self.workerObject.sendTrialStartTrigger_signal.connect(self.sendTrialStartTrigger)
			self.workerObject.sendCoords_signal.connect(self.sendCoords)
			self.workerObject.sendPhotoStimTrig_signal.connect(self.sendPhotoStimTrigger)
			self.workerThread.start()
			self.updateStatusBar('Worker started')
			print('worker thread started')

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

		else:
			self.updateStatusBar('No initialisation provided or PV not connected')
			return


	def outsourcePhotostim(self):
		def resetPhotostimObject():
			self.photostimObject = None
		# setup photostim object if not present
		if not p['FLAG_BLINK_CONNECTED']:
			print('Check blink connection!')
			self.stimFromBlink_checkBox.setChecked(False)
			return

		p['stimFromBlink'] = self.stimFromBlink_checkBox.isChecked()
		if p['stimFromBlink']:
			try:
				# setup photostim object if not present
				if not qtarget.empty():
					with qtarget.mutex:
						qtarget.queue.clear()
					print('qtarget cleared!')

				print('setting up photostim object...')
				kwargs = {"blink": self.bl,"prairie": self.pl,"niPhotostimFullWriter":self.niPhotostimFullWriter,"niPhotoStimFullTask":self.niPhotoStimFullTask}
				self.photostimObject = Photostimer(**kwargs)
				self.photostimObject.moveToThread(self.photostimThread)
				print('photostim object created')

				self.photostimObject.status_signal.connect(self.updateStatusBar)
				self.photostimObject.totalNumPhotostim_signal.connect(self.updateNumPhotostimLabel)
				self.photostimObject.finished_signal.connect(self.photostimObject.deleteLater)  # delete qt (C++) instance later
				self.photostimObject.finished_signal.connect(resetPhotostimObject) # remove python reference to photostimObject
				self.photostimObject.finished_signal.connect(self.photostimThread.exit)
				self.photostimThread.start()
				self.updateStatusBar('Photostimer started')

			except Exception as e:
				print('Error initialising photostimer:')
				print(e)
				self.stimFromBlink_checkBox.setChecked(False)
				p['stimFromBlink'] = False
		else:
			# stop and clean up photostimer thread
			try:
				if self.photostimObject:
					self.photostimThread.wait(1)
					self.photostimThread.quit()
					self.photostimThread.wait() # ensures worker finished before next run can be started
					print('photostim thread stopped')
			except Exception as e:
				print('photostim thread stop error:')
				print(e)


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

	def enableStimTrigger(self):
		p['enableStimTrigger'] = self.enableStimTrigger_checkBox.isChecked()

	def updatePhotostimParameters(self):
		if self.photostimObject:
			self.photostimObject.updatePhotostimParameters()
			print('updated parameters in photostimer')
		else:
			print('did not find photostim object')
#%% mouse events:
	def enterEvent(self,event):
		self.graphicsView.setCursor(Qt.CrossCursor)

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

				if detected and self.removeModeOn:
					if idx in self.removeIdx:
						# reset all remove list works
#                            self.Targetcontour_item.clear()
#                            self.updateTargets()
#                            self.removeIdx = []

						# unselect individual remove items
						selected = self.Targetcontour_item.data
						if self.rejectMaskOn:
							selected_xy = [[value[0], value[1]] for value in selected][self.numTargets+len(p['rejectedX']):]
						else:
							selected_xy = [[value[0], value[1]] for value in selected][self.numTargets:]
						print(len(selected_xy))

						ROI_xy = [ROI_x, ROI_y]
						selected_idx = selected_xy.index(ROI_xy)

						del selected_xy[selected_idx]
						remove_x = [item[0] for item in selected_xy]
						remove_y = [item[1] for item in selected_xy]

						self.updateAllROIs()
						self.Targetcontour_item.addPoints(x = remove_x, y = remove_y, pen = self.RemovePen, size = pointSize)

						self.removeIdx.remove(idx)
					else:
						self.removeIdx.append(idx)
						self.Targetcontour_item.addPoints(x = [ROI_x], y = [ROI_y], pen = self.RemovePen, size = pointSize)

		elif self.rejectModeOn:
			det_dist = 10  # small is preferred here
			detected = 0

			for idx in range(len(p['rejectedX'])):
				detected = abs(x - p['rejectedX'][idx]) <= det_dist and abs(y - p['rejectedY'][idx]) <= det_dist
				if detected:
					del p['rejectedX'][idx]
					del p['rejectedY'][idx]
					self.updateAllROIs()
					return

			if not detected:
				self.Targetcontour_item.addPoints(x = [x], y = [y], pen = self.RejectPen, size = pointSize)
				p['rejectedX'].append(x)
				p['rejectedY'].append(y)

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


#%% Daq functions:
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

	def updatePowerControl(self):
		NI_TTL_NUM_SAMPLES = p['NI_TTL_NUM_SAMPLES']
		NI_STIM_NUM_SAMPLES = int(p['photoDuration']*0.001*p['NI_SAMPLE_RATE'])
		init_output = np.zeros([2, max(NI_TTL_NUM_SAMPLES+1,NI_STIM_NUM_SAMPLES+1)])
		init_output.ravel()[:len(p['NI_1D_ARRAY'])] = p['NI_1D_ARRAY']
		init_output[1,:-1] = 1
		p['NI_2D_ARRAY'] = np.copy(init_output)
		p['NI_UNIT_POWER_ARRAY'] = np.copy(init_output[1,:])
		p['NI_2D_ARRAY'][1,:] = p['NI_UNIT_POWER_ARRAY'] *np.polyval(p['power_polyfit_p'],p['photoPowerPerCell'])

		# reconfigure sample number on task
		self.niPhotoStimFullTask.timing.cfg_samp_clk_timing(rate = p['NI_SAMPLE_RATE'],sample_mode= nidaqmx.constants.AcquisitionType.FINITE, samps_per_chan = max(NI_TTL_NUM_SAMPLES+1,NI_STIM_NUM_SAMPLES+1))

	def connectSenStimTCP(self):
		self.FlAG_stimSOCK_READY = False
		self.stimSOCK = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
		self.stimSOCK.settimeout(3)
		try:
			self.stimSOCK.connect((self.SENSTIM_IP, self.SENSTIM_PORT))
			self.FlAG_stimSOCK_READY = True
			print('Seneory TCP connected')
		except Exception as e:
			print('StimSOCK connection error:')
			print(e)

	def testTTLTrigger(self):
		try:
			self.sendTTLTrigger()
			self.updateStatusBar('Test trigger sent')
		except Exception as e:
			print(e)

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
		if self.StimTask_READY:
			self.niStimWriter.write_many_sample(self.daq_1s_array,10.0)
			print('stim TTL trigger sent')
			return
		if self.FlAG_stimSOCK_READY:
			self.stimSOCK.sendall(bytes(str(p['senStimType']),'utf-8'))
			print('stim TCP trigger sent')
			return


	def sendTrialStartTrigger(self):
		if self.TrialTask_READY:
			self.niTrialStartWriter.write_many_sample(self.daq_1s_array,10.0)
			print('trial start trigger sent')
		else:
			print('tiral start trigger not ready')

	def sendPhotoStimTrigger(self):
		try:
			if not p['SendPowerVolt']:
				self.niPhotoStimWriter.write_many_sample( p['NI_1D_ARRAY'],10.0)
			elif not p['stimFromBlink']:
				self.niPhotostimFullWriter.write_many_sample(p['NI_2D_ARRAY'],10.0)
				while(not self.niPhotoStimFullTask.is_task_done()):
					pass
				self.niPhotoStimFullTask.stop()
				self.niPhotostimFullWriter = stream_writers.AnalogMultiChannelWriter(self.niPhotoStimFullTask.out_stream,auto_start = True)
			else:

				self.photostimCurrentTargets(len(p['currentTargetX']))
			print('photostim trigger sent')
		except Exception as e:
			print('send photostim trigger error')
			print(e)

	def photostimCurrentTargets(self,num_stim_targets):
		ERROR = False
		p['FLAG_SKIP_FRAMES'] = True
		this_volt = np.polyval(p['power_polyfit_p'],p['photoPowerPerCell']*num_stim_targets)
		if not self.bl.send_trigger_power(this_volt):
			print('Photostimer msg: photostim sent from blink')
		else:
			print('photostimCurrentTargets: msg to blink ERROR!')
			ERROR = True
		p['FLAG_SKIP_FRAMES'] = False
		return ERROR
#%% Add/remove ROIs:
	def selectAllROIsClicked(self):
		if (self.selectAll_checkBox.isChecked):
			self.selectAllROIs()

	def removeCells(self, CNN=False, save_init=True, FLAG_KEEP_SELECTED = False, reinitiate = True):
		if self.rejectIdx:
			rejected = True
			self.removeIdx = self.rejectIdx   # just remove rejected cells
			self.rejectIdx = []
		else:
			rejected = False
		print('Removing ' + str(len(self.removeIdx)) + ' cells')

		accepted = self.c['cnm2'].accepted
		print(accepted)

		if len(accepted) == self.c['coms_init'].shape[0]:
			print('coms_init size and accepted size are equal')
		else:
			print('coms_init size and accepted size are NOT equal!')
			return

		if CNN: # deal with global indexing - to do
			thisROIIdx = self.c['K_init']  # same as self.thisROIIdx in most cases
			glb_idx_keep = sorted(list(set(range(thisROIIdx)) - set(self.removeIdx)))  # removeIdx here could be orig idx // mostly local idx though
			glb_idx_accepted = list(set(glb_idx_keep).intersection(accepted))
			idx_keep = [accepted.index(x) for x in glb_idx_accepted]

		else:  # deal with ROI indexing
			if FLAG_KEEP_SELECTED:
				thisROIIdx = self.thisROIIdx
				idx_keep = sorted(list(set(self.removeIdx)))  # removeIdx here is local GUI idx
			else:
				thisROIIdx = self.thisROIIdx
				idx_keep = sorted(list(set(range(thisROIIdx)) - set(self.removeIdx)))  # removeIdx here is local GUI idx

		glb_idx_accepted = [accepted[x] for x in idx_keep]
		print(glb_idx_accepted)

		self.c['idx_components'] = glb_idx_accepted
		if reinitiate: # initiate again, only work if starts with init file
			self.prepareOnacid(save_init=save_init,idx_components =glb_idx_accepted)
		else: # change 'accepted' in cnm and update ROIlist
			cnm2 = deepcopy(self.c['cnm_init'])
			cnm2.accepted = glb_idx_accepted
			self.c['cnm2'] = cnm2
			self.c['coms_init'] = self.c['coms_init'][idx_keep]
			self.initialiseROI()
			print(self.c['coms_init'].shape)
			for i in range(self.c['coms_init'].shape[0]): # cnm2.accepted: #range(K):
				y, x = self.c['coms_init'][i]  # reversed
				self.getROImask(thisROIx = x, thisROIy = y)
			print('ROI list done')
			self.plotOnacidTraces(t=self.c['t_cnm'],online=True)


		if self.showROIsOn:
			self.plotSTAonMasks(None)
			plt.close()

		self.updateTargets()
		self.showROIIdx()
		self.updateStatusBar('Removed ' + str(len(accepted)-len(glb_idx_accepted)) + ' components')

		# clear removeIdx array
		self.removeIdx = []
	def removeButTrigger(self):
		self.removeIdx = self.TriggerIdx
		self.removeCells(FLAG_KEEP_SELECTED = True,reinitiate = not p['FLAG_USING_RESULT_FILE'] )
		# load weights and threshold
		self.ROIsumThresh_doubleSpinBox.setValue(self.TriggerThresh)
		for idx in range(self.thisROIIdx):
			self.ROIlist[idx]["weight"] = self.TriggerWeights[idx]
		self.updateTable()
		self.showROIIdx()

		if p['FLAG_USING_RESULT_FILE'] :
			self.plotOnacidTraces(t=self.c['t_cnm'],online=True)
		else:
			self.plotOnacidTraces(t=self.c['cnm_init'].initbatch)

		self.getValues()


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
			p['ExtraTargetX'] = [] #selected targets (not in the ROIlist)
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
			self.ROIlist[ROIIdx]["weight"] = 1 # default weight is 1
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

			# reset masks
			self.opsin_img_path = 'U:/simulate_movie/20170823_Ch01.tif'
			self.opsin_img = np.array([])
			self.opsin_mask = np.array([])
			self.A_opsin = np.array([])
			self.opsinMaskOn = False
			self.showROIsOn = False
			self.A_loaded = np.array([])
			self.mask_path = ''
			self.reject_mask = np.array([])
			self.rejectIdx = []
			p['rejectedX'] = []
			p['rejectedY'] = []

			p['FLAG_SKIP_FRAMES'] = False
#            self.calcium_img_path = ''
#            self.calciumMaskOn = False
#            self.calcium_img = np.array([])
			# clear trial type
			p['all_trial_types'] = []

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
			self.c['opsin_mask'] = np.array([])  # clear explicitly
			self.c = {}
			self.cnm2 = {}
#            self.UseONACID_checkBox.setEnabled(False)
			self.UseONACID_checkBox.setChecked(False)
			self.imageItem.setImage(self.BlankImage)

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
						print('An ROI overlaps with an extra target. Removing the extra duplicate.')
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

#        self.updateTargetROIs()
		self.updateAllROIs()

	def updateTargetROIs(self):
		self.Targetcontour_item.clear()
		self.Targetcontour_item.addPoints(x = p['currentTargetX'], y = p['currentTargetY'], pen = self.TargetPen, size = self.RoiRadius*2+5)

	def updateRejectROIs(self, display=False):
		self.Targetcontour_item.addPoints(x = p['rejectedX'], y = p['rejectedY'], pen = self.RejectPen, size = self.RoiRadius*2+5)  # added
		print('reject rois updated')


	def updateAllROIs(self):
		self.Targetcontour_item.clear()
		self.updateTargetROIs()
		if self.rejectMaskOn:
			self.updateRejectROIs()

	def removeModeController(self):
		self.removeModeOn = not self.removeModeOn
		if self.removeModeOn:
			self.startRemoveMode()
		else:
			self.removeMode_pushButton.setText('Start remove')
			self.updateAllROIs()
			self.removeIdx = []
			for item in self.groupBoxes:
				item.setEnabled(True)

	def startRemoveMode(self):
		print('Remove mode on!')

		for item in self.groupBoxes:
			item.setEnabled(False)

#        self.removeMode_pushButton.setEnabled(True)
		self.removeMode_pushButton.setText('Quit remove mode')

	def removeSelected(self):
		if self.removeIdx:
			msg = QMessageBox()
			msg.setText("Do you want to remove selected cells?")
			msg.setWindowTitle('pyRTAOI Message')
			msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
			retval = msg.exec_()
			if(retval==QMessageBox.Ok):
				self.removeCells(reinitiate = not p['FLAG_USING_RESULT_FILE'])
			else:
				self.updateAllROIs()
				self.removeIdx = []

		self.removeMode_pushButton.setText('Start remove')
		self.removeModeOn = False

		for item in self.groupBoxes:
			item.setEnabled(True)

		print('Remove mode off')

	def keepSelected(self):
		if self.removeIdx:
			msg = QMessageBox()
			msg.setText("Do you want to remove all other cells?")
			msg.setWindowTitle('pyRTAOI Message')
			msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
			retval = msg.exec_()
			if(retval==QMessageBox.Ok):
				self.removeCells(FLAG_KEEP_SELECTED = True, reinitiate = not p['FLAG_USING_RESULT_FILE'] )
			else:
				self.updateAllROIs()
				self.removeIdx = []

		self.removeMode_pushButton.setText('Start remove')
		self.removeModeOn = False
		for item in self.groupBoxes:
			item.setEnabled(True)
		print('Remove mode off')

#%% Rejection mask
	def loadCalciumImg(self):
		calcium_img_path = str(QFileDialog.getOpenFileName(self, 'Load an avg calcium image', self.movie_folder, 'MPTIFF (*.tif)')[0])
		calcium_img = cm.load(calcium_img_path, subindices = slice(0,1,None))
		self.updateImage(calcium_img)

	def rejectModeController(self):
		self.rejectModeOn = not self.rejectModeOn

		if self.rejectModeOn:
			self.startRejectMode()
		else:
			self.exitRejectMode()


	def startRejectMode(self):
		print('Reject mode on!')
		self.rejectMaskOn = True

		for item in self.groupBoxes:
			item.setEnabled(False)

		self.rejectMode_pushButton.setText('Mask done')


	def exitRejectMode(self):
		try:
			self.createRejectMask()
			self.rejectMode_pushButton.setText('Create mask')

			for item in self.groupBoxes:
				item.setEnabled(True)

			print('Reject mode off')
		except Exception as e:
			print(e)

	def rejectMaskDisplay(self):
		self.rejectMaskOn = not self.rejectMaskOn

		if not self.rejectMaskOn:
			self.updateAllROIs()
			self.rejectMaskDisplay_pushButton.setText('Show mask')
		else:
			self.updateAllROIs()
			self.rejectMaskDisplay_pushButton.setText('Hide mask')


	def createRejectMask(self):
		if self.reject_mask.size and not self.rejectMaskOn:
			self.rejectMaskDisplay()

		radius = self.ROIRadius_spinBox.value() #self.RoiRadius+2.5
		dims = (512,512)   # first create a full-sized mask, then downsample

		reject_mask = np.zeros(dims)

		for idx in range(len(p['rejectedX'])):
			y = p['rejectedX'][idx]   # TODO: check this is correct. why have to invert? if have to, invert these while adding to p?
			x = p['rejectedY'][idx]

			xcoords = []
			ycoords = []
			ROI_mask = np.zeros(dims)

			for pixel_x in range(dims[0]):
				for pixel_y in range(dims[1]):

					dx = pixel_x - x
					dy = pixel_y - y
					dist_squared = dx**2 + dy**2

					if dist_squared <= radius**2:
						xcoords.append(pixel_x)
						ycoords.append(pixel_y)

			ROI_mask[xcoords,ycoords] = 1
			reject_mask += ROI_mask

		try:
			dims_res = self.dims
		except:
			self.ds_factor = self.dsFactor_doubleSpinBox.value()
			dims_res = tuple([round(dim/self.ds_factor) for dim in dims])

#        self.updateImage(reject_mask.astype('u1'))
		reject_mask = cv2.resize(reject_mask, dims_res)   # shapes more round if downsampled after creation
		reject_mask[reject_mask>0] = 1  # keep only 0s and 1s

		if reject_mask != np.zeros(dims):  # store reject_mask only if non-zero
			self.reject_mask = reject_mask.astype('int')
		else:
			self.reject_mask = np.array([])

		self.updateImage(cv2.resize(reject_mask, (512, 512), interpolation=cv2.INTER_CUBIC))

		self.rejectMaskOn = True
		self.rejectMaskDisplay_pushButton.setEnabled(True)

		if self.c != {}:
			self.checkRejectedMask()
			self.c['reject_mask'] = self.reject_mask

	def resetRejectMask(self):
		p['rejectedX'] = []
		p['rejectedY'] = []
		self.rejectIdx = []
		self.updateAllROIs()
		self.createRejectMask()

		if self.c['cnm2'].N > self.InitNumROIs_spinBox.value():
			self.prepareOnacid()

		if self.rejectMaskOn:
			self.rejectMaskDisplay()


	def checkRejectedMask(self):
		if self.c != {}:
			self.rejectIdx = []

			for idx in range(self.c['coms_init'].shape[0]):
				coms = self.c['coms_init'][idx,:]
				x = round(coms[0])  # no need to inverse here
				y = round(coms[1])

				rejected = self.reject_mask[int(x)][int(y)]
				if rejected:
					self.rejectIdx.append(idx)

			self.updateStatusBar(str(len(self.rejectIdx)) + ' cells to be rejected: ' + str([idx+1 for idx in self.rejectIdx])[1:-1])


	def showCellsOnRejectedMask(self):
		if self.reject_mask.size and self.c != {}:
			self.checkRejectedMask()

			pl.figure();
			pl.imshow(self.reject_mask, cmap='gray')
#            pl.plot(self.c['coms_init'][:,1],self.c['coms_init'][:,0],'w.')

			cnm_struct = self.c['cnm2']
			A = cnm_struct.Ab[:, cnm_struct.gnb:cnm_struct.M]
			coords = get_contours(A, self.dims)

			for idx in range(cnm_struct.N):
				c = coords[idx]
				contours = c['coordinates']

				if idx in self.rejectIdx:
					col1 = 'r'
					col2 = 'm'
				else:
					col1 = 'g'
					col2 = 'c'

				pl.plot(*contours.T, c=col1)
				pl.plot(self.c['coms_init'][idx,1], self.c['coms_init'][idx,0],
						c=col2, marker='.', markersize=3)

			pl.axis('off')

		else:
			 self.updateStatusBar('No cells to be shown or reject mask empty')

	def applyRejectMask(self):
		if self.rejectIdx:
			msg = QMessageBox()
			msg.setText("Do you want to reject selected cells from online analysis?")
			msg.setWindowTitle('pyRTAOI Message')
			msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
			retval = msg.exec_()
			if(retval==QMessageBox.Ok):
#                self.prepareOnacid()  # keep rejected cells in init object
#                self.removeIdx = self.rejectIdx   # just remove rejected cells
#                self.rejectIdx = []
				self.removeCells()
			else:
				self.rejectIdx = []

				msg.setText("Reset reject mask?")
				msg.setWindowTitle('pyRTAOI Message')
				msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
				retval = msg.exec_()
				if(retval==QMessageBox.Ok):
					p['rejectedX'] = []
					p['rejectedY'] = []
					self.createRejectMask()

				self.updateAllROIs()
#%% Blink related
	def connectBlink(self):
		self.updateStatusBar('connecting Blink')
		self.bl = bLink(self.BLINK_IP,self.BLINK_PORT)
		print('bl created')
		if self.bl.CONNECTED:
			p['FLAG_BLINK_CONNECTED'] = True
			self.updateStatusBar('Blink connected')
		print('bl connected = '+str(self.bl.CONNECTED))

	def sendCoords(self):
		print('sending coords')
		if p['FLAG_BLINK_CONNECTED']:
			# scale targets coordinates (to refZoom with which the SLM transform matrix was computed);
			# currentTargetXY should be coordinates in 512x512 frames
			currentTargetX = [int((item-255.0)*p['targetScaleFactor']+255.0) for item in p['currentTargetX']]
			currentTargetY = [int((item-255.0)*p['targetScaleFactor']+255.0) for item in p['currentTargetY']]
			if(self.bl.CONNECTED):
				# send a trigger - for measure timing only;
				self.niStimWriter.write_many_sample(p['NI_1D_ARRAY'],10.0)
				if not self.bl.send_coords(currentTargetX, currentTargetY):
					print('Phase mask updated')
				else:
					print('Update phase mask ERROR!')
					p['FLAG_BLINK_CONNECTED'] = False

			else:
				self.updateStatusBar('Send coords faild, check blink connection')
				p['FLAG_BLINK_CONNECTED'] = False
#%% Praire related
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

#%% load/transfer data
	def loadMoviePath(self):
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

	def loadResultPath(self):
		result_path = str(QFileDialog.getOpenFileName(self, 'Load result', self.movie_folder, 'PKL (*.pkl);;All files (*.*)')[0])  # changed '' (default) to movie_folder
		if result_path:
			# load pkl file to main
			with open(result_path, 'rb') as f:
				file_data = pickle.load(f)
			self.resultPath_lineEdit.setText(result_path)
			cnm2 = file_data['cnm2']
			online_C = file_data['online_C']
			coms = file_data['coms']
			self.transDataToMain(cnm2,online_C,coms,file_data['accepted_idx'],file_data['t_cnm'])
			self.ds_factor = file_data['ds_factor']
			# show fov
			self.c['img_norm'] = file_data['img_norm']
			self.showFOV()

			# mark roi on FOV
			NumROIs = cnm2.N
			self.MaxNumROIs = NumROIs
			self.MaxNumROIs_spinBox.setValue(NumROIs)
			self.ALL_ROI_SELECTED == False
			self.thisROIIdx = 0
			self.ROIlist = [dict() for i in range(NumROIs)]
			for ROIIdx in range(NumROIs):
				self.ROIlist[ROIIdx]["ROIx"]= list()
				self.ROIlist[ROIIdx]["ROIy"] = list()
				self.ROIlist[ROIIdx]["ROIcolor"] = self.random_color()
				self.ROIlist[ROIIdx]["threshold"] = list()
				self.ROIlist[ROIIdx]["STA"] = list()
				self.ROIlist[ROIIdx]["weight"] = 1 # default weight is 1

			self.updateTable()
			accepted_idx = cnm2.accepted
			for i in range(0,len(accepted_idx)):
				self.getROImask(coms[i,1], coms[i,0])
			self.showROIIdx()
			print('finished loading result')


	def loadProcData(self):
		stim_type = self.StaStimType_comboBox.currentIndex()
		if stim_type == 0:
			stim_frames_field = 'photo_stim_frames_caiman'
		elif stim_type == 1:
			stim_frames_field = 'stim_frames_caiman'

		self.sta_traces = STATraceMaker.make_sta_file(stim_frames_field = stim_frames_field)

	def transDataToMain(self, cnm_struct, online_C, coms, accepted, t, sta_amp = [], sta_traces=[], plot=True):
#        print(self.c['test']) # changing c inside worker automatically changes it inside mainwindow -- this is sent/shared somehow
#        self.cnm2 = cnm_struct
		self.c['coms'] = coms
		self.c['cnm2'] = cnm_struct
		self.c['online_C'] = online_C
		self.accepted = accepted
		self.t_cnm = t
		self.sta_amp = sta_amp
		self.sta_traces = sta_traces
#        if self.A_opsin.size:
#            print('Checking opsin overlap')
#            self.checkOpsin()
		if plot:
			self.plotOnacidTraces(t, online=True)

	def loadRefMoviePath(self):
		if self.UseOpsinMask_radioButton.isChecked() or self.UseOnacidMask_radioButton.isChecked():
			ref_movie_path = str(QFileDialog.getOpenFileName(self, 'Load reference', self.movie_folder, 'MPTIFF (*.tif);;All files (*.*)')[0])
		else:
			ref_movie_path = str(QFileDialog.getOpenFileName(self, 'Load reference', self.movie_folder, 'PKL(*.pkl);;MPTIFF (*.tif);;All files (*.*)')[0])

		if ref_movie_path and ref_movie_path != 'U:/simulate_movie/20170823.tif':
			self.ref_movie_path = ref_movie_path
			self.movie_folder = os.path.dirname(os.path.dirname(self.ref_movie_path))

			self.refMoviePath_lineEdit.setText(self.ref_movie_path)
			print('Initialisation file selected')

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

		try:
			idx_components = data_loaded['idx_components']  # only display accepted cells
		except:
			idx_components = None

		if issparse(A_loaded):
			A_loaded = np.array(A_loaded.todense())

		if idx_components:
			A_loaded = A_loaded[:,idx_components]

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

#%% caiman initialisation
	def initialiseCaiman(self,FLAG_USING_RESULT_FILE = False):
		try:
			self.deleteTextItems()
			if self.ref_movie_path:
				movie_ext = os.path.splitext(p['refMoviePath'])[1]
				movie_name = os.path.splitext(p['refMoviePath'])[0]
				print(movie_name)
				self.updateStatusBar('Initialising...')
			else:
				self.updateStatusBar('No reference movie selected!')
				return


			# init parameters
			if movie_ext == '.tif' or movie_ext == '.tiff':
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
				if movie_ext == '.tif'or movie_ext == '.tiff':
					K = self.InitNumROIs_spinBox.value()
					print('Starting CNMF initialisation with tiff file')
					lframe, init_values = initialise(self.ref_movie_path, init_method='cnmf',
													 K=K, gSig=gSig, initbatch=500,
													 ds_factor=ds_factor, rval_thr=rval_thr,
													 thresh_overlap=thresh_overlap,
													 save_init=True, mot_corr=p['initMotionCorr'],
													 merge_thresh=merge_thresh,
													 decay_time=0.2, min_SNR=min_SNR)

				elif movie_ext == '.pkl':
					print('Loading initialisation file...')
					init_values = load_object(self.ref_movie_path)
					print(init_values.keys())
					ds_factor = init_values['ds_factor']
					self.dsFactor_doubleSpinBox.setValue(ds_factor)
					print('...done')


			elif opsin_seeded:
				if self.A_opsin.size:
					if movie_ext == '.tif' or movie_ext == '.tiff':
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
														 save_init=True, mot_corr=p['initMotionCorr'],
														 min_SNR=min_SNR) #, CNN_filter=True)
					else:
						self.updateStatusBar('A non-tif file was provided - select a tif movie to initialise with a mask')
						return
				else:
					self.updateStatusBar('No opsin mask created!')


			elif mask_seeded:
				if movie_ext == '.tif' or movie_ext == '.tiff':
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
													 thresh_overlap=0.1, save_init=True, mot_corr=p['initMotionCorr'],
													 merge_thresh=0.85)

				else:
					self.updateStatusBar('A non-tif file was provided - select a tif movie to initialise with a mask')
					return

			print('initialising c')
			self.c = deepcopy(init_values) # caiman object

			self.c['cnn_removed_idx'] = []
			if (not '_init_' in movie_name) and movie_ext == '.pkl':
				FLAG_USING_RESULT_FILE = True

			if FLAG_USING_RESULT_FILE:
				self.c['coms_init'] = self.c['coms']
				print(self.c['coms_init'].shape)
				self.c['cnm_init'] = self.c['cnm2']
				try:
					self.c['coms_init_orig'] = self.c['coms_init_orig']
					self.c['cnm_init_orig'] = self.c['cnm_init_orig']
				except Exception as e:
					print('init orig error')
					print(e)
				self.dims = self.c['cnm2'].dims
				cnm = init_values['cnm2']
				self.c['cnm2'].t = init_values['t_cnm']
				# make sure everything called from c in worker is also in result
			else:
				self.c['coms_init_orig'] = self.c['coms_init']  # keep original com data
				self.c['cnm_init_orig'] = self.c['cnm_init']
				self.dims = self.c['cnm_init'].dims
				cnm = init_values['cnm_init']



			# set parameters in gui
			print('setting parameters in gui...')
			if movie_ext == '.pkl':
				cell_radius = round(cnm.gSig[0]*ds_factor)
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
			self.min_SNR = min_SNR
			self.rval_thr = rval_thr
			self.merge_thresh = merge_thresh
			self.thresh_overlap = thresh_overlap
			self.ds_factor = ds_factor
			self.rejectIdx = []
			self.c['rejected_idx'] = [] # temp; remove
			print('...done')

			if opsin_seeded:
				thresh_cnn = 0.1 # default thresh for c1v1 mask
				self.CNNFilter_doubleSpinBox.setValue(thresh_cnn)
				self.filterResults(init=True)
			else:
				self.prepareOnacid(save_init=False, using_result_file = FLAG_USING_RESULT_FILE)

			self.updateStatusBar('Initialisation complete')
			p['FLAG_USING_RESULT_FILE'] = FLAG_USING_RESULT_FILE
		except Exception as e:
			print('caiman initialisation error')
			print(e)

	def prepareOnacid(self, show_results=True, save_init=True, using_result_file = False,idx_components = None):
		print('this is prepareOnacid function')
		try:
			# Reset previous settings
			try:
#				self.c['cnm2']
				self.ROIcontour_item.clear()
				self.deleteTextItems()
				self.plotItem.clear()
				self.resetFigure()
				self.thisROIIdx = 0
				self.resetROIlist()
				self.updateTable()
			except: pass

			# Prepare object for OnACID
#			idx_components = self.c['idx_components']
			path_to_cnn_residual = os.path.join(caiman_datadir(), 'model', 'cnn_model_online.h5')

			# check if loaded pkl file is result file or init file
			if using_result_file: # get cnm from result file
				print('Prepare Onacid from result file')
				cnm2 = deepcopy(self.c['cnm_init'])
				coms = self.c['coms_init'].copy()

			else: # get cnm from initialisation file
				print('Prepare Onacid from init file')
				cnm2 = deepcopy(self.c['cnm_init'])
				cnm2._prepare_object(np.asarray(self.c['Yr']), self.c['T1'],
								 self.c['expected_comps'], idx_components=idx_components,
								 min_num_trial=2, N_samples_exceptionality=int(self.c['N_samples']),
								 path_to_model=path_to_cnn_residual, sniper_mode=False)

	#        cnm2.max_num_added = 5  # max number of cells added per onacid frame --> doesn't look like more cells are added per frame
				cnm2.opsin = []
#				coms = self.c['coms_init_orig'].copy()
				# all rejected deleted already
				accepted = list(range(0,cnm2.N))
				# store info on rejected cells and mask
				cnm2.accepted = accepted
#				self.c['accepted'] = accepted

				if self.rejectMaskOn:
					if not self.c['reject_mask'].size:
						self.c['reject_mask'] = self.reject_mask # if no mask provided with the init file, assign the one created in the interface
					else:
						self.reject_mask = self.c['reject_mask'] # if mask provided, change currently stored reject_mask

				cnm2.t = cnm2.initbatch
				cnm2.max_num_added = p['MaxNumROIs']

#			if idx_components is not None:
#				print('discarded components from coms_init')
#				self.c['coms_init'] = coms[idx_components]

			print(cnm2.Yr_buf)
			print(cnm2.Yr_buf.shape)

			# reset buffers
			cnm2.Yr_buf.cur = 0
			cnm2.Yr_buf.max_ = cnm2.minibatch_shape

			cnm2.Yres_buf.cur = 0
			cnm2.Yres_buf.max_ = cnm2.minibatch_shape # this is inibatch shape

			cnm2.rho_buf.cur = 0
			cnm2.rho_buf.max_ = cnm2.minibatch_shape # this is inibatch shape


			K = cnm2.N
			self.c['cnm2'] = cnm2
			self.c['K_init'] = K


			try: # temp; remove soon
				self.c['thresh_cnn']
			except:
				self.c['thresh_cnn'] = 0


			if show_results:
				# Extract number of cells detected
				self.InitNumROIs_spinBox.setValue(K)

				if self.MaxNumROIs_spinBox.value() < K+5:
					self.MaxNumROIs_spinBox.setValue(K+10)

				print('Number of components initialised: ' + str(K))
				self.MaxNumROIs = self.MaxNumROIs_spinBox.value()

				self.InitNumROIs = K
				self.opsinMaskOn = False
				self.imageItem.setImage(cv2.resize(self.c['img_norm'],(512,512),interpolation=cv2.INTER_CUBIC)) # display FOV

				self.initialiseROI()
				print(self.c['coms_init'].shape)
				for i in range(self.c['coms_init'].shape[0]): # cnm2.accepted: #range(K):
					y, x = self.c['coms_init'][i]  # reversed
					self.getROImask(thisROIx = x, thisROIy = y)
				print('ROI list done')

				self.CNNFilter_doubleSpinBox.setValue(self.c['thresh_cnn'])

				self.UseONACID_checkBox.setEnabled(True)
				self.UseONACID_checkBox.setChecked(True)
				if self.selectAll_checkBox.isChecked():
					self.selectAllROIs()
				self.updateStatusBar('Preparing object completed')

				if self.A_opsin.size:
					print('Checking opsin overlap')
					self.checkOpsin()
				if self.rejectMaskOn:
					if self.c['reject_mask'].size:
						print('Checking reject mask')
						self.checkRejectedMask()

				self.showROIIdx()

				if using_result_file:
					self.plotOnacidTraces(t=self.c['t_cnm'],online=True)
				else:
					self.plotOnacidTraces(t=self.c['cnm_init'].initbatch)

			if save_init:
				self.saveInitResults()
			print('done preparing onAcid')
		except Exception as e:
			print('prepareOnacid error:')
			print(e)


	def saveInitResults(self):              # save changed init object
		init_values_new = deepcopy(self.c)  # copy to avoid messing with self.c
		del init_values_new['cnm2']

#        init_values_new['idx_components'] = orig_keep_idx   # will be used during prepare_object to keep only cells of interest
		init_values_new['coms_init'] = self.c['coms_init_orig']  # keep the original

		save_path = self.c['save_path'][:-4] + '_filtered.pkl'
		init_values_new['save_path'] = save_path

		print('Saving new init file as ' + str(save_path))
		save_object(init_values_new, save_path)


	def resetFiltering(self, show_results=True):
#		del self.c['cnm2']  commented out 20190711
		self.c['thresh_cnn'] = 0
		self.c['cnn_removed_idx'] = []

		# do not include manually deleted cells
		idx_components = sorted(list(set(range(self.c['K'])) - set(self.c['removed_idx']) - set(self.c['rejected_idx'])))

		self.c['idx_components'] = idx_components
		self.prepareOnacid(show_results=show_results)
		print('Onacid object prepared')


	def filterResults(self, init=False):
		try:
			thresh_cnn = self.CNNFilter_doubleSpinBox.value()
			self.c['CNN_filter'] = bool(thresh_cnn)
			self.rejectIdx = []  # this is rejectIdx from rejection mask; have to reset
			print('Using CNN filter? '+ str(self.c['CNN_filter']))

			if self.c['thresh_cnn'] == thresh_cnn:
				print('threshold did not change, return')
				return

			if self.c['CNN_filter']:
				if len(self.c['cnn_removed_idx']):
					self.resetFiltering(show_results=False)

				print('Filtering components')
				print('CNN threshold: ', thresh_cnn)

				if init:
					cnm_struct = self.c['cnm_init']
					A = cnm_struct.A
				else:
					cnm_struct = self.c['cnm2']
					A = cnm_struct.Ab[:, cnm_struct.gnb:cnm_struct.M]

				gSig = cnm_struct.gSig #  expected half size of neurons

				print('using model in: ',os.path.join(caiman_datadir(), 'model', 'cnn_model'))
				predictions, final_crops = evaluate_components_CNN(
					A, self.dims, gSig, model_name=os.path.join(caiman_datadir(), 'model', 'cnn_model'))

				predictions[np.isnan(predictions)] = 0  # exclude nan
				idx_keep_CNN = np.where(predictions[:, 1] >= thresh_cnn)[0].tolist()
				idx_rem_CNN = np.where(predictions[:,1] < thresh_cnn)[0].tolist()
				print(str(len(idx_keep_CNN)) + ' rois left post filtering')

				self.c['thresh_cnn'] = thresh_cnn
	#            self.c['idx_components'] = idx_keep_CNN
				print('idx(global) remove', idx_rem_CNN)

				# debug select cells to remove
				self.removeIdx = idx_rem_CNN # this is global index need to match roi idx
	#			self.removeModeOn = True
	#			for idx in self.removeIdx:
	#				ROI_x = self.ROIlist[idx]["ROIx"]
	#				ROI_y = self.ROIlist[idx]["ROIy"]
	#				self.Targetcontour_item.addPoints(x = [ROI_x], y = [ROI_y], pen = self.RemovePen, size =  self.RoiRadius*2+5)

#				self.startRemoveMode()
				# end of debug

	#			self.removeIdx = idx_rem_CNN
				self.removeCells(CNN=True,reinitiate = init )

			# reset previous CNN filter
			else:
				if len(self.c['cnn_removed_idx']):
					print('Resetting filter results')
					self.resetFiltering()
		except Exception as e:
			print(e)
#%% Opsin mask
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

			# downsample params
			radius = round(radius/ds_factor)
			min_area = min_area/ds_factor**2

			if radius % 2 == 0:
				radius += 1

#            self.cellRadius_spinBox.setValue(radius) # can be hidden - rescaling more important

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
		try:
			if Ain is None:
				cnm_struct = self.c['cnm2']
				A = cnm_struct.Ab[:, cnm_struct.gnb:cnm_struct.M]
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
			use_mask = self.useMask_checkBox.isChecked()
			self.use_mask = use_mask

			c = 0 # temp
			for cell in accepted: # range(onacid_mask.shape[-1]):
				if use_mask:
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

				else:
	#                opsin_pos.append(self.opsin_mask[coms[0]][coms[1]])  # an easy way; checking just com
	#                print('opsin', opsin_pos)

					try:
						coms = self.c['coms'][c]   # coms received from worker post onacid
						print('post onacid coms')
						c += 1
					except:
						coms = self.c['coms_init'][cell]

					coms = [int(round(com)) for com in coms]

#                    x = int(round(self.ROIlist[c]["ROIy"]))
#                    y = int(round(self.ROIlist[c]["ROIx"]))

					offset = 2

					x = coms[0]
					y = coms[1]

					x_all = [x, x+offset, x-offset]
					y_all = [y, y+offset, y-offset]

					opsin_count = 0
					count = 0
					for i in range(len(x_all)):
						for j in range(len(y_all)):
							op = self.opsin_mask[x_all[i]][y_all[j]]
							opsin_count += op
							count += 1

					accept = opsin_count/count >= self.opsin_thresh
					opsin.append(accept)


			if Ain is None:
				# store info on opsin
				self.c['cnm2'].opsin = opsin  # True/False based on threshold
				self.c['overlap'] = overlap   # value
				self.c['opsin_mask'] = self.opsin_mask
				self.c['opsin_thresh'] = self.opsin_thresh
				self.onacid_mask = onacid_mask

			if self.selectAll_checkBox.isChecked():
				self.selectAllROIs()  # update ROIs selected

		except Exception as e:
			print(e)
		return onacid_mask

#%% Display functions:
	def random_color(self):
		r = random.randrange(0, 255)
		g = random.randrange(0, 255)
		b = random.randrange(0, 255)
		return QColor(r, g, b)

	def updateTable(self):
		# re-create the whole table in the table to those in the ROIlist
		self.Thresh_tableWidget.clear()  # empty the table
		self.Thresh_tableWidget.setHorizontalHeaderLabels(['X','Y','Thresh','W','STA'])  # column names disappear with clearing
		self.Thresh_tableWidget.setColumnCount(5)

		if self.thisROIIdx > self.MaxNumROIs:  # do not allow to miss detected cells
			self.MaxNumROIs = self.thisROIIdx + 10
			print('Extending max num ROIs to ' + str(self.MaxNumROIs))
			self.MaxNumROIs_spinBox.setValue(self.MaxNumROIs)

		NumROIs = self.MaxNumROIs
		self.Thresh_tableWidget.setRowCount(NumROIs)
		self.Thresh_labels = [QTableWidgetItem() for i in range(NumROIs)]
		self.STA_labels = [QTableWidgetItem() for i in range(NumROIs)]
		self.X_labels = [QTableWidgetItem() for i in range(NumROIs)]
		self.Y_labels = [QTableWidgetItem() for i in range(NumROIs)]
		self.W_labels = [QTableWidgetItem() for i in range(NumROIs)]


		for ROIIdx in range(self.thisROIIdx):

			thisROIx = self.ROIlist[ROIIdx]["ROIx"]
			thisROIy = self.ROIlist[ROIIdx]["ROIy"]
			thisW = self.ROIlist[ROIIdx]["weight"]
			qcolor = self.ROIlist[ROIIdx]["ROIcolor"]

			tempbrush = QBrush(qcolor)
			self.X_labels[ROIIdx].setForeground(tempbrush)
			self.Y_labels[ROIIdx].setForeground(tempbrush)
			self.Thresh_labels[ROIIdx].setForeground(tempbrush)
			self.W_labels[ROIIdx].setForeground(tempbrush)

			self.X_labels[ROIIdx].setText(str('{:.0f}'.format(thisROIx)))
			self.Y_labels[ROIIdx].setText(str('{:.0f}'.format(thisROIy)))
			self.W_labels[ROIIdx].setText(str('{:.3f}'.format(thisW)))

			# make items read-only
			self.X_labels[ROIIdx].setFlags(Qt.ItemIsEditable)
			self.Y_labels[ROIIdx].setFlags(Qt.ItemIsEditable)

			# put item in table
			self.Thresh_tableWidget.setItem(ROIIdx,0,self.X_labels[ROIIdx])
			self.Thresh_tableWidget.setItem(ROIIdx,1,self.Y_labels[ROIIdx])
			self.Thresh_tableWidget.setItem(ROIIdx,2,self.Thresh_labels[ROIIdx])
			self.Thresh_tableWidget.setItem(ROIIdx,3,self.W_labels[ROIIdx])

		self.Thresh_tableWidget.resizeColumnsToContents



	def tableItemChanged(self,item):
		row = item.row()
		col = item.column()
		if col==3: #only deal with change in weights
			self.ROIlist[row]["weight"] = item.text()
		else:
			return

	def displayImg(self,event):
		if self.opsinMaskOn:
			self.updateImage(cv2.resize(np.squeeze(self.opsin_img), (512, 512), interpolation=cv2.INTER_CUBIC))

#        elif self.calciumMaskOn:
#             self.updateImage(cv2.resize(self.calcium_img, (512, 512), interpolation=cv2.INTER_CUBIC))

	def displayMask(self,event):
		if self.opsinMaskOn:
			self.updateImage(cv2.resize(np.squeeze(self.opsin_mask).astype('u1'), (512, 512), interpolation=cv2.INTER_CUBIC))

#        elif self.calciumMaskOn:
#            self.updateImage(cv2.resize(self.reject_mask.astype('uint8'), (512, 512), interpolation=cv2.INTER_CUBIC))

	def showIdx(self):
		self.deleteTextItems()
		self.showROIIdx()

	def showFOV(self):
		print('show FOV clicked')
		self.opsinMaskOn = False
		self.showROIsOn = False
#        self.calciumMaskOn = False
		self.imageItem.setImage(cv2.resize(self.c['img_norm'],(512,512),interpolation=cv2.INTER_CUBIC)) # display FOV

	def plotSTAonMasks(self,sta_amp):
		# show cells detected by caiman
		# make a image with sta levels
		print('show STA on mask clicked')
		self.opsinMaskOn = False
		cnm = self.c['cnm2']
		print(cnm.gnb)

		A = cnm.Ab[:, cnm.gnb:cnm.M]

		if issparse(A):
			A = np.array(A.todense())
		else:
			A = np.array(A)

		print('shape A',A.shape)

		# do not show rejected cells
		nr = self.c['cnm2'].N
		print('N = ',nr)
		accepted = self.c['cnm2'].accepted
#		print('accepted index')
#		print(accepted)

		# use sta value, otherwise use one
		try:
			if sta_amp is None: # will show scaled amplitude of A
				print('input sta is None')
				sta_amp = np.ones((nr,))*255
			elif len(sta_amp) == 0:
				print('making sta trace')
				frames_to_average = range(p['staAvgStartFrame']+p['staPreFrame'],p['staAvgStopFrame']+p['staPreFrame'])
				sta_amp = STATraceMaker.plotSTAtraces(self.sta_traces,frames_to_average = frames_to_average, IF_PLOT = False) # load from temp file
				sta_amp = sta_amp(accepted)
			else: # normalise within component before multiply with sta
				print('normalise A')
				for i in accepted: # range(nr):
					A[:,i] = A[:,i]/sum(A[:,i])

#			print('using sta:')
#			print(sta_amp)
#			print(len(sta_amp))
		except Exception as e:
			print(e)
		# put value into mask
		try:
			cellMask = np.zeros((cnm.dims[1]*cnm.dims[0],))

			j = 0
			for i in accepted: # range(np.minimum(len(sta_amp),nr)):
				if not np.isnan(sta_amp[j]):
					cellMask+=A[:,i].flatten()*sta_amp[j]
					j +=1

			cellMask2D = np.reshape(cellMask,cnm.dims,order='F')
			cellMask2D = cellMask2D/max(cellMask)*255

			# show on a separate plot
	#		norm = plt.Normalize(0,1)
	#		im = plt.imshow(norm(cellMask2D),aspect = 'equal',cmap = 'Greys')
	#		plt.colorbar(im, orientation='horizontal')
	#		plt.show()
			cellMask2D = np.repeat(cellMask2D[:,:,None],3,axis=-1)
			self.imageItem.setImage(cv2.resize(cellMask2D, (512, 512), interpolation=cv2.INTER_CUBIC))
			self.showROIsOn = True
		except Exception as e:
			print('adding value to mask error')
			print(e)
	def plotSTA(self):
		frames_to_average = range(p['staAvgStartFrame']+p['staPreFrame'],p['staAvgStopFrame']+p['staPreFrame'])
		self.sta_amp=STATraceMaker.plotSTAtraces(self.sta_traces,frames_to_average = frames_to_average) # load from temp file

		#  add figure to layout
		self.staTracesFig = fig
		self.staTracesCanvas = FigureCanvas(self.staTracesFig)
		self.plot_gridLayout.addWidget(self.staTracesCanvas)
		self.updateStatusBar('STA traces plotted on plot tab')

	def resetFigure(self):
		try:
	#        plt.close('all')
			self.ax1.cla()
			self.ax2.cla()
			self.ax3.cla()
			self.axcomp.cla()

			self.figCanvas.draw()
	#            self.figCanvas.update()
			self.figCanvas.flush_events()

			g = gc.collect()
			print('gc: ', g)
		except: pass

	def plotOnacidTraces(self, t, online=False):
		try:
			try: self.resetFigure()
			except: pass

			cnm_struct = self.c['cnm2']

			if online:
				C,f = self.c['online_C'][cnm_struct.gnb:cnm_struct.M], self.c['online_C'][:cnm_struct.gnb]
			else:
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

			cell_mask = 1
			if cell_mask:
				img = np.reshape(np.array(A.mean(axis=1)), (d1, d2), order='F')
			else:
				img = self.c['Cn_init']

			all_ = list(range(0, cnm_struct.N))

			try:
				accepted = cnm_struct.accepted
				rejected = [item for item in all_ if item not in accepted]
				print('accepted idx:')
				print(accepted)
			except:
				accepted = all_
				rejected = []

			self.axcomp = self.fig.add_axes([0.05, 0.05, 0.9, 0.03])
			self.ax1 = self.fig.add_axes([0.05, 0.55, 0.4, 0.4])
			self.ax3 = self.fig.add_axes([0.55, 0.55, 0.4, 0.4])
			self.ax2 = self.fig.add_axes([0.05, 0.1, 0.9, 0.4])
			self.s_comp = Slider(self.axcomp, 'Component', 0, nr + nb - 1, valinit=0, facecolor = [.5, .5, .5])

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
						self.ax1.set_title('Spatial component ' + str(i+1) + '(glb idx ' + str(j) + ')')
					else:
						self.ax1.set_title('Rejected spatial component ' + str(i-len(accepted)+1)+ 'glb idx' + str(j))
					self.ax1.axis('off')


					self.ax2.cla()
					self.ax2.plot(np.arange(T), Y_r[j], 'c', linewidth=2)
					self.ax2.plot(np.arange(T), C[j], 'r', linewidth=2)
					if not rej:
						self.ax2.set_title('Temporal component ' + str(i+1))
					else:
						self.ax2.set_title('Rejected temporal component ' + str(i-len(accepted)+1))

					if online:
						self.ax2.legend(labels=['Filtered raw data', 'Online trace'])
					else:
						self.ax2.legend(labels=['Filtered raw data', 'Inferred trace'])


					self.ax3.cla()
					if vmax > 0:
						self.ax3.imshow(img, interpolation='None', cmap=pl.cm.gray, vmax=vmax)
					else:
						self.ax3.imshow(img, interpolation='None', cmap=pl.cm.gray)
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
		except Exception as e:
			print('plot exception:')
			print(e)

	def arrow_key_image_control(self, event):
#		print(event.key())
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

	def showCellsOnMask(self):
		if self.A_opsin.size and self.c != {}:
			self.checkOpsin()  # for overlap update
			opsin = self.c['cnm2'].opsin
			use_mask = self.useMask_checkBox.isChecked()

			# visualise all comp masks
			if use_mask:
				summed_A = np.hstack((self.A_opsin, self.onacid_mask))
				summed_mask = np.reshape(np.array(summed_A.sum(axis=1)), self.dims, order='F')
				pl.figure();pl.imshow(summed_mask)
				pl.colorbar()
				pl.axis('off')

			# visualise coms and contours on opsin mask
			else:
				cnm_struct = self.c['cnm2']
				A = cnm_struct.Ab[:, cnm_struct.gnb:cnm_struct.M]
				coords = get_contours(A, self.dims)

				pl.figure(); pl.imshow(self.opsin_mask, cmap='gray')

				try:
					coms = self.c['coms']   # coms received from worker post onacid
				except:
					coms = self.c['coms_init']

				try:
					cell = 0
					for idx in self.c['cnm2'].accepted: # range(cnm_struct.N):
						c = coords[idx]
						contours = c['coordinates']
						if opsin[cell]:
							col = 'g'
						else:
							col = 'r'

						pl.plot(*contours.T, c=col)
						pl.plot(round(coms[cell,1]), round(coms[cell,0]),
								c=col, marker = '.')
						cell += 1

					pl.axis('off')
				except Exception as e:
					print(e)

		elif self.A_opsin.size and self.A_loaded.size:
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
		try:
			print('show num rois = '+str(self.thisROIIdx))
			font = QFont("Arial", 12, QFont.Bold)
			for idx in range(self.thisROIIdx): # self.c['cnm2'].accepted:
				thisText = pg.TextItem(text=str(idx+1), color=self.ROIlist[idx]["ROIcolor"].getRgb()[:3],
									   html=None, anchor=(0,0), border=None, fill=None, angle=0, rotateAxis=None)
				thisText.setPos(self.ROIlist[idx]["ROIx"],self.ROIlist[idx]["ROIy"])
				thisText.setFont(font)
				self.graphicsView.addItem(thisText)
		except Exception as e:
			print(e)
#            i += 1

	def updateROI(self,img):
		# import pdb
		# pdb.set_trace()
		self.updateImage(img)

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
			self.STA_labels[i].setText(str('{0:.1f}'.format(sta_arr[i])))
			self.Thresh_tableWidget.setItem(i,3,self.STA_labels[i])
		print('updated sta value')

	def updateFrameInfo(self,FrameIdx):
		self.CurrentFrameIdx_label.setText(str(FrameIdx))

	def updateNumPhotostimLabel(self, totalNumFrames):
		self.totalNumPhotostim_label.setText(str(totalNumFrames))

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
#%% close and cleanup
	def closeEvent(self,event):
		# override closeEvent method in Qt
		print('form closing, PV connection = ' + str(self.FLAG_PV_CONNECTED))
		if self.FlAG_stimSOCK_READY:
			self.stimSOCK.close()
			print('sensory stim socket closed')

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
			self.niTrialStartWriter.write_one_sample(0,10.0)
			self.niTrialStartTask.stop()
			self.niTrialStartTask.close()
			del self.niTrialStartTask
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
			print('closed photostim full trigger task')

		except Exception as e:
			print('photostim full task clean up error')
			print(e)

		try:
			self.niPhotoStimWriter.write_one_sample(0,10.0)
			self.niPhotoStimTask.stop()
			self.niPhotoStimTask.close()
			del self.niPhotoStimTask
			print('closed photostim single trigger task')
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


#%% load/save configuration
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

	def saveResults(self):
		try:
			save_dict = dict()
			save_dict['cnm2'] = self.c['cnm2'] #self.cnm2
			save_dict['com_count_init'] = self.c['K_init']
			save_dict['coms'] = self.c['coms']
			save_dict['accepted'] = self.accepted
			save_dict['t_cnm'] = self.t_cnm
			save_dict['params'] = p
			save_object(save_dict, p['saveResultPath'])
		except Exception as e:
			self.updateStatusBar('Save results error: '+str(e))

	def browseResultPath(self):
		saveResultsPath = QFileDialog.getSaveFileName (self, 'Select save path', p['saveResultPath'], '.pkl')
		p['saveResultPath'] = saveResultsPath[0]+saveResultsPath[1]
		self.saveResultPath_lineEdit.setText(p['saveResultPath'])

	def loadTrialOrder(self):
		# load trial order from pybehav
		filepath = str(QFileDialog.getOpenFileName(self, 'Load trial sequence', p['saveResultPath'], '*.txt')[0])
		if filepath:
			arr = np.genfromtxt(filepath, delimiter=',')
			p['trialOrder'] = arr[0].astype(np.int)

			self.updateStatusBar('Loaded trial type from trialOrder file')
			p['useExternalTrialOrder'] = True
			self.useExternalTrialOrder_checkBox.setCheckState(Qt.Checked)
			self.updateTrialType()

	def loadStimOrder(self):
		# load target idx list and trial order from matlab
		filepath = str(QFileDialog.getOpenFileName(self, 'Load stimOrder mat', p['saveResultPath'], '*.mat')[0])
		print(filepath)
		if filepath:
			try:
				[self.TargetIdxList,p['trialOrder'],p['fixedTargetSeqX'],p['fixedTargetSeqY']] = get_stimOrder(filepath)
				self.updateStatusBar('Loaded trial type from stimOrder file')
				p['useExternalTrialOrder'] = True
				p['ConditionOnStimType'] = False
				self.useExternalTrialOrder_checkBox.setCheckState(Qt.Checked)
				self.ConditionOnStimType_checkBox.setCheckState(Qt.Unchecked)
				self.updateTrialType()
				self.updateStatusBar('Loaded StimOrder file')

			except Exception as e:
				print(e)
				return
	def previewStimOrder(self):
		try:
#			for i in range(len(p['trialOrder'])):
				i = random.randint(1,len(p['trialOrder']))
				p['currentTargetX'] = p['fixedTargetSeqX'][i].copy()
				p['currentTargetY'] = p['fixedTargetSeqY'][i].copy()
				if len(p['currentTargetX'])>0:
					self.CurrentFrameIdx_label.setText(str(p['sen_stim_frames'][i]))
					self.updateTargetROIs()
				else:
					print('this trial as no photostim')


		except Exception as e:
			print(e)
			return
	def loadTriggerConfig(self):
		trigger_config_path = str(QFileDialog.getOpenFileName(self, 'Load trigger configeration', self.movie_folder, 'mat (*.mat);;All files (*.*)')[0])
		if trigger_config_path:
			self.triggerConfigPath_lineEdit.setText(trigger_config_path)
			[self.TriggerIdx,self.TriggerWeights,self.TriggerFrames, self.TriggerThresh,self.TargetIdx] = get_triggertargets_params(trigger_config_path)
#			preview trigger cell positions
			print('loaded trigger idx:')
			print(self.TriggerIdx)
			print('loaded target idx:')
			print(self.TargetIdx)

			# show trigger cells
			self.Triggercontour_item.clear()
			self.Triggercontour_item.addPoints(x = [self.ROIlist[i]["ROIx"] for i in self.TriggerIdx], y = [self.ROIlist[i]["ROIy"] for i in self.TriggerIdx], pen = self.TriggerPen, size = self.RoiRadius*2+3)

			# config target cells
#			self.clearTargets()
#			p['currentTargetX'].extend([self.ROIlist[i]["ROIx"] for i in self.TargetIdx])
#			p['currentTargetY'].extend([self.ROIlist[i]["ROIy"] for i in self.TargetIdx])
			self.updateTargets()

	def loadTargetCentroid(self):
		self.target_img_path = str(QFileDialog.getOpenFileName(self, 'Load target centroids', self.movie_folder, 'MPTIFF (*.tif);;All files (*.*)')[0])
		if self.target_img_path:
			self.targetCentroidsPath_lineEdit.setText(self.target_img_path)
			target_image = cm.load(self.target_img_path)
			cols,rows = np.where(target_image>0)

			# reset targets
			self.clearTargets()

			# delet existing extra targets if they are too close to the new loaded ones
			det_dist = 10
			idx_remove = []

			for roi_idx in range(len(rows)):
				for extra_idx in range(len(p['ExtraTargetX'])):
						x = p['ExtraTargetX'][extra_idx]
						y = p['ExtraTargetY'][extra_idx]

						detected = abs(x - rows[roi_idx]) <= det_dist and abs(y - cols[roi_idx]) <= det_dist
						if detected:
							print('An ROI overlaps with an extra target. Removing the extra duplicate.')
							idx_remove.append(extra_idx)

			for idx in sorted(idx_remove, reverse=True):
				del p['ExtraTargetX'][idx]
				del p['ExtraTargetY'][idx]

			p['ExtraTargetX'].extend(rows)
			p['ExtraTargetY'].extend(cols)

			self.updateTargets()

	def saveCurrentTargetCentroids(self):
		save_img = np.zeros([512, 512], dtype=np.float)
		for idx in range(self.numTargets):
			save_img[int(p['currentTargetY'][idx])][int(p['currentTargetX'][idx])] = 255
		filepath = str(QFileDialog.getSaveFileName(self, 'Save as preset...', 'Tiffs', 'Tiff file (*.tif)')[0])
		cv2.imwrite(filepath,save_img)

#%% calibration check (burn spots)
	def calCheck(self):
		# burn patterns of spots in a sequence; even patterns are blanked out (0 power on sample)
		ERROR_FLAG = False
		current_dir = os.getcwd()
		print(current_dir)
		pattern_dir = current_dir+'/Patterns'
#        pattern_dir =r"C:\Users\User\Desktop\pyRTAOI-rig\Patterns"
		file_list = glob.glob(os.path.join(pattern_dir,'*.tif'))
		num_patterns = len(file_list)
		xcoord_list = []
		ycoord_list = []

		for file_name in file_list:
			this_pattern = tifffile.imread(file_name)
			this_coords = np.where(this_pattern==255)
			xcoord_list.append(this_coords[1])
			ycoord_list.append(this_coords[0])

		print(str(num_patterns)+' patterns loaded')

		power_array = np.full(num_patterns,p['photoPowerPerCell'])
		power_array[1::2] =0

		# show expected burnt locations
		burntPen = pg.mkPen(color = (255,112,75), width = 2)
		skipPen = pg.mkPen(color = (64,64,64), width = 2)
		skipx = [item for sublist in xcoord_list[1::2] for item in sublist ]
		skipy = [item for sublist in ycoord_list[1::2] for item in sublist ]
		burntx =  [item for sublist in xcoord_list[::2] for item in sublist ]
		burnty =  [item for sublist in ycoord_list[::2] for item in sublist ]
		self.Targetcontour_item.clear()
		self.Targetcontour_item.addPoints(x = burntx, y = burnty, pen = burntPen, size = self.RoiRadius*2+5)
		self.Targetcontour_item.addPoints(x = skipx, y = skipy, pen = skipPen, size = self.RoiRadius*2)

		# confirm if start buring
		msg = QMessageBox()
		msg.setText("Start burning spots?")
		msg.setWindowTitle('calCheck Message')
		msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
		retval = msg.exec_()
		if(retval==QMessageBox.Cancel):
			ERROR_FLAG = True
			self.Targetcontour_item.clear()
			return
		# send to blink and trigger photostim
		for i in range (num_patterns):
			p['currentTargetX'] = xcoord_list[i]
			p['currentTargetY'] = ycoord_list[i]
			num_stim_targets = len(xcoord_list[i])

			send_start = time.time()
			if not p['stimFromBlink']:
				# send from mainwindow
				if p['FLAG_BLINK_CONNECTED']:
					self.sendCoords()
					print('coords sent')
				else:
					ERROR_FLAG = True
				print('send coords time'+ str("%.4f"%(time.time()-send_start)))


				if (not ERROR_FLAG) and self.POWER_CONTROL_READY:
					print('total power:'+str(power_array[i]*num_stim_targets))
					try:
						if power_array[i]==0:
							p['NI_2D_ARRAY'][1,:] = p['NI_UNIT_POWER_ARRAY']*0
						else:
							p['NI_2D_ARRAY'][1,:] = p['NI_UNIT_POWER_ARRAY']*np.polyval(p['power_polyfit_p'],power_array[i]*num_stim_targets)
						self.sendPhotoStimTrigger()
					except Exception as e:
						print(e)
						ERROR_FLAG = True
					print('stim trigger sent')
				else:
					ERROR_FLAG = True
			else:
				# send from photostimer and blink
				qtarget.put([p['currentTargetX'].copy(),p['currentTargetY'].copy()])
				if power_array[i]==0:
					ERROR_FLAG = self.photostimObject.updateNewTargets()
				else:
					ERROR_FLAG = self.photostimObject.photostimNewTargets()

			if not ERROR_FLAG:
				print('pattern'+str(i)+' stimulated')
			else:
				print('calcheck error, check daq and/or blink connection!')
#%% photoexcitability check
	def configSeqStim(self):
		# setup configuration for stimulating current targets one by one in randomised sequence
		self.photoProto_comboBox.setCurrentIndex(CONSTANTS.PHOTO_SEQUENCE)
				# update trial parameters
		self.interStimInterval_spinBox.setValue(30)
		self.staPreFrame_spinBox.setValue(30)
		self.staPostFrame_spinBox.setValue(60)
		self.stimStartFrame_spinBox.setValue(150)
		self.enableStimTrigger_checkBox.setCheckState(Qt.Unchecked)
		self.getValues()

		p['photo_sequence_idx'] = []
		self.updateFixTargets() # save current targets to fix target list
		num_targets = len(p['fixedTargetX'])
		print('seq stim num targets:')
		print(num_targets)
		num_repeats = p['numberRepeats']  # using stim condif
		target_idx = [i for i in range(num_targets)]
		tot_num_trials = num_targets*num_repeats
		p['photo_stim_frames'] =  p['stimStartFrame']+p['interStimInterval']*np.arange(0,tot_num_trials)
		self.numberStims_spinBox.setValue(tot_num_trials)
		for r in range(num_repeats):
			p['photo_sequence_idx'].extend(random.sample(target_idx,num_targets))
		print(p['photo_sequence_idx'])
		self.getValues()


#%% drag drop
	def dragEnterEvent( self, event ):
		data = event.mimeData()
		urls = data.urls()
		if ( urls and urls[0].scheme() == 'file' ):
			event.acceptProposedAction()

	def dragMoveEvent( self, event ):
		data = event.mimeData()
		urls = data.urls()
		if ( urls and urls[0].scheme() == 'file' ):
			event.acceptProposedAction()

	def dropEvent( self, event ):
		data = event.mimeData()
		urls = data.urls()
		if ( urls and urls[0].scheme() == 'file' ):
			# for some reason, this doubles up the intro slash
			if os.name == 'nt':
				filepath = str(urls[0].path()[1:])
			else:
				filepath = str(urls[0].path())
			droppedIMG = cv2.imread(filepath,cv2.IMREAD_ANYDEPTH)
			self.imageItem.setImage(cv2.resize(droppedIMG,(512,512),interpolation=cv2.INTER_CUBIC))
			self.opsinMaskOn = False
			print('current image:'+str(filepath))



#%%  Main entry to program.  Sets up the main app and create a new window.
if __name__ == '__main__':
	# create Qt application
	app = QCoreApplication.instance()
	if app is None:  # stops kernel crashing
		app = QApplication([])

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

	# for debug use
	global test_cnm
	test_cnm = {}

   # global data buffer
	global qbuffer
	qbuffer = Queue(maxsize = 0)

	# global photostim target que
	global qtarget
	qtarget = Queue(maxsize = 0)

	# global stop flag
	global STOP_JOB_FLAG
	STOP_JOB_FLAG = False


	# create main window
	window = MainWindow()

	# show it and bring to front
	window.show()
	window.raise_()

	# start the app
	sys.exit(app.exec_())
	app.deleteLater()


#	# launch program
#	try:
#		main(sys.argv)
#	except Exception as e:
#		logger.exception(e)
#		print(str(e))
