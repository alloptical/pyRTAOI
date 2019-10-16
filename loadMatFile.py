
import scipy.io as sio
import os
import numpy as np
import time



def get_power_params():

# load PowerFile (saved out by matlab SavePowerFile)
# map power to voltage
# control voltage = (outMax/displayMax)*PV value;
# outMax and displayMax are defined in prairie view configuration.xml
# add bruker1 dropbox path to system variable, as BRUKER1_DROPBOX_PATH


	# parameters
	outMax = 1
	displayMax = 1000
	powerfile_name = os.environ['BRUKER1_DROPBOX_PATH'] + '/PowerUtilities/Bruker1_PowerFile.mat'
	mat_file = sio.loadmat(powerfile_name)
	power_file = mat_file['power_file']
	power_polyfit_p = np.ravel(power_file['p'])[0][0]
	slop = outMax/displayMax
	power_polyfit_p *=slop
	return power_polyfit_p


def get_triggertargets_params(file_full_name):
	'''
	read file saved by generate_cell_idx_file for texture discrimination task
	get cell idx, weight and threshold used for population trigger

	'''
	mat_file = sio.loadmat(file_full_name)
	trigger_file = mat_file['output']
	trigger_idx = trigger_file['trigger_idx'][0][0].flatten()-1
	target_idx = trigger_file['target_idx'][0][0].flatten()-1
	trigger_weights = trigger_file['trigger_weights'][0][0].flatten()
	trigger_frames = trigger_file['trigger_frames'][0][0].flatten()
	trigger_thresh = trigger_file['trigger_thresh'][0][0].flatten()

	return trigger_idx,trigger_weights,trigger_frames,trigger_thresh, target_idx


def get_stimOrder(file_full_name):
		# this is for loading PHOTO_FIX_SEQUENCE
		mat_file = sio.loadmat(file_full_name)
		stimOrder = mat_file['pyrtaoi_stimOrder']
		input_target_idx_list = stimOrder['target_idx_list'][0][0]
		input_target_centroid_x = stimOrder['target_centroid_x'][0][0]
		input_target_centroid_y = stimOrder['target_centroid_y'][0][0]
		trialOrder = stimOrder['trialOrder'][0][0].flatten()
		target_idx_list = []
		target_centroid_x = []
		target_centroid_y = []

		for i in range(len(trialOrder)):
			temp_target_idx = input_target_idx_list[i]
			target_idx = [i for i in temp_target_idx if i!=0] # exclude dummy targets
			if all(idx ==0 for idx in target_idx): # exclude trials without photostim
				target_idx_list.append([]) # set lists of zeros to empty
				target_centroid_x.append([])
				target_centroid_y.append([])
			else:
				target_idx_list.append([idx-1 for idx in target_idx]) # -1 for python indexing
				target_centroid_x.append([x for x in input_target_centroid_x[i].tolist() if x>=0])
				target_centroid_y.append([y for y in input_target_centroid_y[i].tolist() if y>=0])

		return target_idx_list,trialOrder,target_centroid_x,target_centroid_y


if __name__ == '__main__':

	# test loadpowerfile:
	# power = 0
	# power_polyfit_p = get_power_params()
	# print(power_polyfit_p)

	# timer_start = time.time()
	# volt = np.polyval(power_polyfit_p,power)
	# print('convert time: '+str("%.4f"%(time.time()-timer_start)))
	# print(volt)

	# test loadtriggerfile
	file_name = r'D:\TextureData\data\cb217\20190805\pyrtaoi_results\\20190805_cb217_t_0006_rtaoi_DS_2.0_OnlineProc_153400proc_OutputParams_20191010_1138.mat'
	[trigger_idx,trigger_weights,trigger_frames,trigger_thresh, target_idx] = get_triggertargets_params(file_name)




