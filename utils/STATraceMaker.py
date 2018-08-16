# copied scripts from STAMovieMaker 

import os
import glob
import _pickle
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import colorsys

root = tk.Tk()
root.withdraw()

def make_sta_file(file_full_name = '',save_full_name = '', stim_frames = [], 
    pre_samples = 30, post_samples = 90,stim_frames_field = 'photo_stim_frames_caiman',
    num_init_frames = 500):
    # get file names
    if file_full_name == '':
        file_full_name = filedialog.askopenfilename()

    print(file_full_name)

    with open(file_full_name, 'rb') as f:
        # file opened
        file_data = _pickle.load(f)

    if save_full_name == '':
        save_full_name = file_full_name.replace('.pkl','_sta.npy')


    # get data
    cnm = file_data['cnm2']
    init_com_count = file_data['init_com_count']

    cnm_C, f = cnm.C_on[cnm.gnb:cnm.M], cnm.C_on[:cnm.gnb]
#    cnm_C = C

    # get frame detected
    try:
        frame_detected = file_data['frame_added'] # as detected in gui
    except:
        frame_detected_roi = [x[0] for x in cnm.time_neuron_added]
        frame_detected_ = [x[1] for x in cnm.time_neuron_added]
        first_new_roi_idx = [x for x in range(len(frame_detected_roi)-1) if np.diff(frame_detected_roi)[x] <0]
        first_new_roi_idx = first_new_roi_idx[0]+1
        frame_detected =[num_init_frames]*init_com_count +frame_detected_[first_new_roi_idx:]
    print(frame_detected)


    if stim_frames ==[]:
        try:
            stim_frames = file_data[stim_frames_field]
            stim_frames = [x+num_init_frames for x in stim_frames] # add number of initialisation frames, 500 is default
            print('stim frames :')
            print(stim_frames)
        except Exception as e:
            print('load stim frames error')
            print(e)

    # make STA template
    sta_template = np.arange(-pre_samples, post_samples)
    num_trials = len(stim_frames)
    num_frames = C.shape[1]
    num_rois = C.shape[0]

    # set F before detection to nan
    for i in range(num_rois):
        C[i,range(frame_detected[i])] = np.nan


    # plot trace and stim time
    plot_range = np.arange(stim_frames[0]-pre_samples,stim_frames[-1]+post_samples)
    stim_trace = np.zeros(num_frames)
    stim_trace[stim_frames] = np.nanmax(C)
    plt.figure
    for i in range(num_rois):
        plt.plot(C[i,plot_range])
    plt.plot(stim_trace[plot_range],color='black')
    plt.show('hold')

    # get raw sta traces
    all_trials_sta_frames = []
    for stim_frame_idx in stim_frames:
        all_trials_sta_frames.append(sta_template + stim_frame_idx)

    trials = np.empty([num_rois, num_trials, len(sta_template)])
    trials[:] = np.nan
    norm_trials = np.copy(trials)
    print(trials.shape)

    for i in range(num_rois):
        for j, trial_sta_frames in enumerate(all_trials_sta_frames):
            for k, frame_idx in enumerate(trial_sta_frames):
                trials[i,j,k] = C[i,frame_idx]

    # normalise by mean value during pre_sample
    for i in range(num_rois):
        for j in range(num_trials):
            this_baseline = np.nanmean(trials[i,j,range(pre_samples)])
            norm_trials[i,j,:] = [ value - this_baseline for value in trials[i,j,:]]

    # save to file
    np.save(save_full_name,norm_trials)

    return norm_trials

def plotSTAtraces(sta_traces=[],frames_to_avgerage = range(30,60)):
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
    print(sta_traces.shape)
    sta_trial_avg = np.nanmean(sta_traces,1)

    sta_trial_avg_amp = np.nanmean(sta_trial_avg[:,frames_to_avgerage],1)

    num_rois = sta_traces.shape[0]
    num_trials = sta_traces.shape[1]
    plot_rows = np.ceil(np.sqrt(num_rois))
    plot_cols = np.ceil((num_rois+1)/plot_rows)

    trial_colors = get_trial_color(num_trials)


    fig,axs = plt.subplots(int(plot_rows),int(plot_cols),
                           facecolor='w', edgecolor='k', frameon=False)
    axs = axs.ravel()
    for i in range(num_rois):
        for j in range(num_trials):
            axs[i].plot(sta_traces[i,j,:],color = trial_colors[j])
        axs[i].plot(sta_trial_avg[i,:], color=(.3,.3,.3))

        axs[i].set_title('ROI'+str(i+1))
    axs[num_rois-1].get_xaxis().set_visible(True)

    for i in range(num_rois+1):
        axs[i].set_aspect('auto')
        axs[i].set_frame_on(False)
        axs[i].set_adjustable('box')
        axs[i].get_xaxis().set_visible(False)

    # show color map in the last plot
    for i in range(num_trials):
        axs[num_rois].bar(i+1,1,1,color= trial_colors[i])
        axs[num_rois].get_yaxis().set_visible(False)
        axs[num_rois].get_xaxis().set_visible(True)
        axs[num_rois].set_xlabel('Trials')

    for i in range(num_rois+1,int(plot_rows*plot_cols)):
        fig.delaxes(axs[i])

    plt.show('hold')
    print(sta_trial_avg_amp)
    return sta_trial_avg_amp

def get_trial_color(num_stims):
    colors = []
    for i in range(num_stims):
        H = np.float32(i / (num_stims+1))
        S = .3
        V = .9
        colors.append(colorsys.hsv_to_rgb(H,S,V))
    return colors


if __name__ == '__main__':
    # stim_frames = np.array([238,    548,859, 1169,1480,1790,2101,2411])
    # stim_frames = stim_frames + 200
    # sta_traces = make_sta_file(r'C:\Users\Zihui\Documents\Test movies\20171229_OG245_t-052_Cycle00001_Ch2_substack1-2700_DS_2.0_OnlineProc.pkl', stim_frames = stim_frames)
    sta_traces = make_sta_file()
    
    plotSTAtraces(sta_traces)