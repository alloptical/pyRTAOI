import os
import glob
import _pickle
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

root = tk.Tk()
root.withdraw()

def make_sta_file(file_full_name = '',save_full_name = '', stim_frames = [], pre_samples = 30, post_samples = 90):
	# get file names
	if file_full_name == '':
		file_full_name = filedialog.askopenfilename()
	print(file_full_name)

	with open(file_full_name, 'rb') as f:
		# file opened
		file_data = _pickle.load(f)

	p = Path(file_full_name)
	file_dir = str(p.parents[0])
	file_name = p.name
	save_name = file_name.replace('.pkl','_sta.mat')
	if save_full_name == '':
		save_full_name = os.path.join(file_dir,save_name)


	# get data
	cnm = file_data['cnm2']
	C, f = cnm.C_on[cnm.gnb:cnm.M], cnm.C_on[:cnm.gnb]
	cnm_C = C
	print(C.shape)

	if stim_frames ==[]:
		try:
			stim_frames = file_data['sensory_stim_frames']
		except Exception as e:
			print('load sensory stim frames error')
			print(e)

	# # plot trace and stim time
	# plot_range = np.arange(0,2700)
	# stim_trace = np.zeros(plot_range.shape)
	# stim_trace[stim_frames] = np.max(C)
	# plt.figure
	# plt.plot(np.transpose(C[:,plot_range]))
	# plt.plot(stim_trace,color='black')
	# plt.show('hold')


	# make STA template
	sta_template = np.arange(-pre_samples, post_samples)
	num_trials = len(stim_frames)
	num_frames = C.shape[1]
	num_rois = C.shape[0]

	all_trials_sta_frames = []
	for stim_frame_idx in stim_frames:
		all_trials_sta_frames.append(sta_template + stim_frame_idx)

	trials = np.zeros([num_rois, num_trials, len(sta_template)], dtype=np.float32)

	for i in range(num_rois):
		for j, trial_sta_frames in enumerate(all_trials_sta_frames):
			for k, frame_idx in enumerate(trial_sta_frames):
				trials[i,j,k] = C[i,frame_idx]
	print(trials.shape)

	# save to file
	np.save(save_name,trials)

	return trials

def plotSTA(sta_traces=[]):
	if sta_traces == []:
		# load file in Temp folder
		try:
			filename = os.getcwd()+'/Temp/sta_traces.npy'
			sta_traces = np.load(filename)
		except:
			try:
				filename = str(QFileDialog.getOpenFileName(self, 'Select STA traces file', '', 'npy (*.npy);;All files (*.*)')[0])
			except:
				return

	sta_trial_avg = np.nanmean(sta_traces,1)

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
		axs[i].set_aspect('auto')
		axs[i].set_frame_on(False)
		axs[i].set_adjustable('box')

	plt.show('hold')


if __name__ == '__main__':
	stim_frames = np.array([238,	548,859, 1169,1480,1790,2101,2411])
	stim_frames = stim_frames + 200
	sta_traces = make_sta_file(r'C:\Users\Zihui\Documents\Test movies\20171229_OG245_t-052_Cycle00001_Ch2_substack1-2700_DS_2.0_OnlineProc.pkl', stim_frames = stim_frames)
	plotSTA(sta_traces)