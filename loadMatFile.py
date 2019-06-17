
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


def get_trigger_params(file_full_name):
	'''
	read file saved by generate_cell_idx_file for texture discrimination task
	get cell idx, weight and threshold used for population trigger

	'''
	mat_file = sio.loadmat(file_full_name)
	trigger_file = mat_file['output']
	trigger_idx = trigger_file['trigger_idx'][0][0].flatten()-1
	trigger_weights = trigger_file['trigger_weights'][0][0].flatten()
	trigger_frames = trigger_file['trigger_frames'][0][0].flatten()
	trigger_thresh = trigger_file['trigger_thresh'][0][0].flatten()

	return trigger_idx,trigger_weights,trigger_frames,trigger_thresh



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
	file_name = r'D:\TextureData\pyrtaoi_proc_data\20190429_CB183\output_files\20190429-CB183_OutputParams_20190605_1405.mat'

	[trigger_idx,trigger_weights,trigger_frames,trigger_thresh] = get_trigger_params(file_name)
	print(trigger_idx)
	print(trigger_weights)
	print(trigger_frames)

