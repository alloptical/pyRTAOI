'''
# load PowerFile (saved out by matlab SavePowerFile) 
# map power to voltage
	 control voltage = (outMax/displayMax)*PV value; 
	 outMax and displayMax are defined in prairie view configuration.xml
	 add bruker1 dropbox path to system variable, as BRUKER1_DROPBOX_PATH
'''

import scipy.io as sio
import os
import numpy as np
import time



def get_power_params():
	# parameters 
	outMax = 1.9
	displayMax = 1000
	powerfile_name = os.environ['BRUKER1_DROPBOX_PATH'] + '/PowerUtilities/Bruker1_PowerFile.mat'
	mat_file = sio.loadmat(powerfile_name)
	power_file = mat_file['power_file']
	power_polyfit_p = np.ravel(power_file['p'])[0][0]
	slop = outMax/displayMax
	power_polyfit_p *=slop
	return power_polyfit_p




if __name__ == '__main__':
	power = 100
	power_polyfit_p = get_power_params()
	print(power_polyfit_p)

	timer_start = time.time()
	volt = np.polyval(power_polyfit_p,power)
	print('convert time: '+str("%.4f"%(time.time()-timer_start)))
	print(volt)