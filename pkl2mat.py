# convert caiman results saved by pickle to .mat
# will save in the same dir as the same name
import os
import _pickle
import scipy.io
from scipy.sparse import issparse
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
import numpy as np
root = tk.Tk()
root.withdraw()

def pkl2mat(file_full_name = '',save_full_name = '', IF_CNMF = True):
	if file_full_name == '':
		file_full_name = filedialog.askopenfilename()
	print(file_full_name)
	with open(file_full_name, 'rb') as f:
		# file opened
		file_data = _pickle.load(f)
	p = Path(file_full_name)
	file_dir = str(p.parents[0])
	file_name = p.name
	save_name = file_name.replace('.pkl','.mat')
	if save_full_name == '':
		save_full_name = os.path.join(file_dir,save_name)

	save_dict = dict()

	if IF_CNMF:
	# parse cnmf object
		cnm = file_data['cnm2']
		A, b = cnm.Ab[:, cnm.gnb:], cnm.Ab[:, :cnm.gnb].toarray()
		C, f = cnm.C_on[cnm.gnb:cnm.M], cnm.C_on[:cnm.gnb]
		S = cnm.S
		t = cnm.t



		if issparse(A):
			A = np.array(A.todense())
		save_dict['cnm_A'] = A
		save_dict['cnm_b'] = b
		save_dict['cnm_C'] = C #cnm.C
		save_dict['cnm_f'] = f #cnm.f
		save_dict['cnm_t'] = f #cnm.t
		save_dict['opsin_positive'] = cnm.opsin


		# cell traces
		save_dict['noisyC'] = cnm.noisyC
		save_dict['cnm_dims'] = cnm.dims

		# other
		save_dict['cnm_N'] = cnm.N
		save_dict['cnm_gnb'] = cnm.gnb
		save_dict['num_frames_init'] = cnm.initbatch
	

	# parse params
	params = file_data['p']

	# copy parameters
	param_names = ['ds_factor', 'photo_stim_frames_caiman','K','min_SNR','gSig',
				   'rval_thr','thresh_overlap','merge_thresh','expected_comps',
				   'frame_added','online_photo_frames','online_photo_targets','repeated_idx','accepted_idx','rejected_idx',
				   't_cnm','coms','opsin_mask','overlap','stim_frames_caiman','online_C','filt_C','online_thresh','tottime','t_init',
				   'frames_skipped','sensory_stim_frames','frame_detected','init_com_count','trialOrder','trialVar','photo_sequence_idx','photoDuration','keep_prev',
				   'online_traj','bs_level','sd_level','ROIsumThresh','ROIw','offsetFrames','CNN_thresh','CNN_predictions','online_oppo_photo_frames']   # added record of new cells in pyrtaoi

	for param in param_names:
		try:
			save_dict[param] = file_data[param]
		except:
			try:
				save_dict[param] = params[param]
			except:
				print('Parameter not found:'+param)



	print('saving...')
	scipy.io.savemat(save_full_name, mdict = save_dict)
	print('saved as '+save_full_name)


if __name__ == '__main__':
	pkl2mat()