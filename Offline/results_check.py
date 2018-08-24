# Offline check of saved pyrtoai data

import os
import glob
import tkinter as tk
import copy
import numpy as np
import matplotlib.pyplot as pl
import cv2
from scipy.sparse import issparse

import caiman as cm
from caiman.utils.utils import load_object
from caiman.utils.visualization import view_patches_bar, plot_contours

#%% Load pkl-ed cnm object

# online recordings
#folder = r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\tests on rig\20180811\pyrtaoi_results'
#
#files = sorted(glob.glob(os.path.join(folder,'*.pkl')))
#
#file = files[1]  # choose file

file = r'F:\Data\Zoe\20180822\pyrtaoi_results\20180822_OG299_t_0008_rtaoi_DS_2.0.tirtaoi_OnlineProc_DS_2.0_180806.pkl'

unsaved = 1
save_figs = 0

# example of 'correct' file i.e. with frame_extra_added data
#file = r'//live.rd.ucl.ac.uk/ritd-ag-project-rd00g6-mhaus91/forPat/samples/example1\pyrtaoi_results\20171229_OG245_t-052_Cycle00001_Ch2_substack1-2700_DS_1.5_rtaoi_OnlineProc_20180808_144712.pkl'
#unsaved = 0

cnm_object = load_object(file)
locals().update(cnm_object)

A, b = cnm2.Ab[:, cnm2.gnb:cnm2.M], cnm2.Ab[:, :cnm2.gnb]
C, f = cnm2.C_on[cnm2.gnb:cnm2.M,:t_cnm], cnm2.C_on[:cnm2.gnb,:t_cnm]
dims = cnm2.dims

# Convert A to numpy array
if issparse(A):
    A = np.array(A.todense())
else:
    A = np.array(A)

#YrA=cnm2.YrA # residual signal
noisyC = cnm2.noisyC[:,:t_cnm]
YrA = noisyC[cnm2.gnb:cnm2.M] - C
Y_r=cnm2.YrA + cnm2.C # raw signal

try:
    frames_skipped = framesSkipped
except: pass

#%% Visualise OnACID output
view_patches_bar([], A, C, b, f, dims[0], dims[1], YrA=YrA, img=None) #img=Cn for background frame, None for cell mask

#%% Check cell mask and init coms
cell_mask = np.reshape(np.array(A.mean(axis=1)), dims, order='F')
#cell_mask = np.reshape(A[:,1], (d1,d2), order='F') # for single cells

cell_mask_ = cv2.resize(cell_mask, (512, 512), interpolation=cv2.INTER_CUBIC)
pl.figure();
pl.imshow(cell_mask_)

coms_ = coms*ds_factor
for cell in range(init_com_count):
    x = coms_[cell,1]
    y = coms_[cell,0]
    pl.plot(x, y, '.r')
    pl.text(x-10, y-10, str(cell), fontsize=12)

#pl.plot(coms_[:init_com_count,1], coms_[:init_com_count,0], '.r')

#vmax = np.percentile(cell_mask_, 95)
#pl.imshow(cell_mask_, interpolation='None', cmap=pl.cm.gray, vmax=vmax)  # for more 'binarised' look

#%% Find original indeces of cells deleted
coms_init_orig = c['coms_init']

orig_keep_idx = []

for pair in coms: # coms_init
    idx = np.where(pair == coms_init_orig)[0]
    if idx.size:
        if idx[0] == idx[1]:
            orig_keep_idx.append(idx[0])
        else:
            print('Different indeces - check')

orig_K = coms_init_orig.shape[0]
orig_removed_idx = list(set(range(orig_K)) - set(orig_keep_idx))

#%% Compare frame_extra_added with time neuron added
if not unsaved:
    time_added = [time[1]-cnm2.initbatch for time in cnm2.time_neuron_added]

    diff = np.array(time_added)[2+K:] - np.array(frame_extra_added)  # for some reason time added for 3 and 4 repeated for that example
    print(diff)  # all 0s - good

#%% Extract info on cells to get frame-by-frame time
try:
    K = init_com_count
except: pass

#    if file == files[2]:  # last file: incorrect K value saved for 0807 session
#        K = 9

tot_onacid_cells = cnm2.N - K
onacid_cells_ix = list(range(K,cnm2.N))

accepted_ix = accepted
rejected = cnm2.N - len(accepted)
rejected_ix = list(set(range(cnm2.N)) ^ set(accepted))

onacid_accepted_ix = copy.copy(onacid_cells_ix)
for ix in rejected_ix:
    if ix in onacid_accepted_ix:
        onacid_accepted_ix.remove(ix)

# unsaved - used for first online files
if unsaved:
    offset = 0
    if cnm2.N != len(cnm2.time_neuron_added):
        print('Some repeats in time_neuron_added: manually check and choose which values to keep!')
        offset = 8   # set this to get the right values

    accepted_spotted = []
    for cell_ix in onacid_accepted_ix:
#        trace = C[cell_ix,:t_cnm]
#        t_spotted = np.where(trace>0)[0][0] - cnm2.initbatch
        t_spotted = cnm2.time_neuron_added[cell_ix+offset][1] - cnm2.initbatch
        accepted_spotted.append(t_spotted)

    rejected_spotted = []
    for cell_ix in rejected_ix:
#        trace = C[cell_ix,:t_cnm]
#        t_spotted = np.where(trace>0)[0][0] - cnm2.initbatch
        t_spotted = cnm2.time_neuron_added[cell_ix+offset][1] - cnm2.initbatch
        rejected_spotted.append(t_spotted)

    all_spotted = np.sort(accepted_spotted+rejected_spotted)

# if frame_extra_added saved
else:
    try:
        frame_extra_added = frame_added[K:]
    except: pass

    accepted_spotted = []
    rejected_spotted = []
    all_spotted = np.array(frame_extra_added)

    i = 0;
    for cell_ix in onacid_cells_ix:
        if cell_ix in accepted_ix:
            accepted_spotted.append(frame_extra_added[i])
        elif cell_ix in rejected_ix:
            rejected_spotted.append(frame_extra_added[i])
        i += 1

    accepted_spotted = np.array(accepted_spotted)
    rejected_spotted = np.array(rejected_spotted)

#%% Plot time of each frame to see the delay introduced by various functionalities
tottime = np.array(tottime)
tottime_ms = tottime*1000
duration = tottime.shape[0]

mean_rate = np.mean(tottime_ms)# sum(tottime)/duration
std_time = np.std(tottime_ms)
print('\nMean rate was ' + str(mean_rate)[:5] + ' +- ' + str(std_time)[:4] + ' ms per frame')

shape_refresh = cnm2.minibatch_shape
shape_update = np.arange(shape_refresh,duration+1,shape_refresh)-1

gui_refresh = 30
gui_update = np.arange(gui_refresh,duration+1,gui_refresh)-1

root = tk.Tk()
width = root.winfo_screenwidth()
height = root.winfo_screenheight()

display_unsaved = 1
if not unsaved: display_unsaved = 0;

font = 18
mark = 10

fig = pl.figure(figsize=(width/100., height/100.), dpi=100)
pl.plot(tottime_ms)
pl.plot([0,duration],[33,33],'--m',label = 'Image acquisition time')
pl.plot(shape_update, np.ones([int((t_cnm-cnm2.initbatch)/shape_refresh),1])*9.5, 'b.', markersize=mark, label='Cell shapes updated')
pl.plot(gui_update, np.ones([int((t_cnm-cnm2.initbatch)/gui_refresh),1])*9.5, 'y.', markersize=5, label='GUI display updated')

distinguish = 0 # flag for whether to distinguish between accepted and rejected cells

#if display_unsaved:
if distinguish:
    pl.plot(accepted_spotted,np.ones([len(accepted_spotted),1])*11,'g*', markersize=mark, label='New accepted cell detected')
    pl.plot(rejected_spotted,np.ones([len(rejected_spotted),1])*11,'r*', markersize=mark, label='New rejected cell detected')
else:
    pl.plot(all_spotted,np.ones([len(all_spotted),1])*80,'g*', markersize=6, label='New cell detected')


pl.title('Time per frame for the online pipeline\n' + r'$mean = ' + str(mean_rate)[:4] + '\pm' + str(std_time)[:3] + ' ms$',fontsize=font)

pl.xlabel('Frame',fontsize=font)
pl.ylabel('Processing time (ms)',fontsize=font)
pl.legend(fontsize=font)

try:
    figManager = pl.get_current_fig_manager()
    figManager.window.showMaximized()
except:
    mng = pl.get_current_fig_manager()
    mng.frame.Maximize(True)

if save_figs:
    if display_unsaved:
        fig.savefig(file[:-4] + '_time_plot_unsaved.png')
    else:
        fig.savefig(file[:-4] + '_time_plot.png')


#%% Time-per-frame histogram
fig = pl.figure();
#fig = pl.figure(figsize=(width/100., height/100.), dpi=100)
weights = np.ones_like(tottime_ms)/float(len(tottime_ms))
detailed = 1

# time for all frames
hist, bins_, patches = pl.hist(tottime_ms, weights=weights,bins=15)
pl.axvline(x=33,color='m',label='Acquisition rate')

if detailed:
    # time for shape update
    shape_update_time = []
    shape_update_weights = np.ones(len(shape_update))*weights[0] # all weights same
    for frame in shape_update:
        shape_update_time.append(tottime_ms[frame])

    pl.hist(shape_update_time,weights=shape_update_weights,bins=bins_,color='b',label='Shape update')


    # time for gui update
    gui_update_time = []
    gui_update_weights = np.ones(len(gui_update))*weights[0] # all weights same
    for frame in gui_update:
        gui_update_time.append(tottime_ms[frame])

    pl.hist(gui_update_time,weights=gui_update_weights,bins=bins_,color='y',label='GUI update')


    # time for new cell detection
    new_cell_time = []
    new_cell_weights = np.ones(len(all_spotted))*weights[0] # all weights same
    for frame in all_spotted:
        new_cell_time.append(tottime_ms[frame])

    pl.hist(new_cell_time,weights=new_cell_weights,bins=bins_,color='g',label='New cell detected')


pl.title('Time per frame')
pl.xlabel('Time (ms)')
pl.ylabel('Fraction of frames')
pl.legend()


if save_figs:
    if detailed:
        fig.savefig(file[:-4] + '_time_hist_det.png')
    else:
        fig.savefig(file[:-4] + '_time_hist.png')

#%% Visualise noise
cell = 5
pl.figure();
pl.suptitle('Signal and noise for cell ' + str(cell+1))

pl.subplot(121)
pl.plot(YrA[cell,:],'r',label='noise')
pl.plot(C[cell,:],'g',label='C')
y_max = int(pl.ylim()[1] + 2)
y_min = int(pl.ylim()[0] - 0.5)
pl.ylim(y_min, y_max)
pl.legend()

pl.subplot(122)
pl.plot(noisyC[cell+1,:],'b',label='noisyC')
pl.legend()
pl.ylim(y_min, y_max)

#%% Check photostim frames and effects on traces
pl.close('all')

for cell in range(init_com_count,cnm2.N):
#cell = 10
    pl.figure(); pl.plot(C[cell,:])
    pl.title('Cell ' + str(cell))

    try:
        frames_skipped = framesSkipped
    except: pass

    for frame in photo_stim_frames_caiman:
        pl.axvline(x=frame+cnm2.initbatch,color='r')

# also check neuropil signal
pl.figure(); pl.plot(cnm2.C_on[0,:t_cnm])
pl.title('Neuropil signal')
for frame in photo_stim_frames_caiman:
    pl.axvline(x=frame+cnm2.initbatch,color='r')

#%% Check online traces
pl.close('all')

save_plots = 0

spotted = list(accepted_spotted.copy() + cnm2.initbatch)  # accepted for now because online traces are only of those accepted
spotted.reverse()

i = 0; # separate counter for online_C if online_C saves only accepted traces


for cell in accepted: # range(cnm2.N):
    pl.figure()

    pl.subplot(311); pl.plot(cnm2.noisyC[cell+1,:t_cnm]); pl.ylabel('noisyC')
    pl.title('cell ' + str(cell));
    ylim_ = pl.ylim()
    pl.subplot(312); pl.plot(cnm2.C_on[cell+1,:t_cnm]); pl.ylabel('C_on')
    pl.ylim(ylim_)

    try:
        if store_all_online:
            pl.subplot(313); pl.plot(online_C[cell,:t_cnm]); pl.ylabel('online_C')
        else:
            pl.subplot(313); pl.plot(online_C[i,:t_cnm]); pl.ylabel('online_C')
            i += 1;
        pl.ylim(ylim_)
    except: pass
    if cell >= K:
        fr = spotted.pop()
        print(fr)
        pl.plot([fr, fr], pl.ylim(), 'r')   # something off with plotting this time now

    if save_plots:
#        subfolder = 'on_cell_plots' # 'off_cell_plots' #
#        save_folder = os.path.join(folder, subfolder)
        save_folder = folder
        if not os.path.exists(save_folder):
            os.makedirs(save_folder)

        pl.savefig(os.path.join(save_folder, 'cell_' + str(cell)))
        pl.close()