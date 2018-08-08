# Offline check of saved pyrtoai data

import os
import glob
import tkinter as tk
import numpy as np
import matplotlib.pyplot as pl

import caiman as cm
from caiman.utils.utils import load_object
from caiman.utils.visualization import view_patches_bar, plot_contours

#%% Load pkl-ed cnm object
folder = r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\tests on rig\20180807\pyrtaoi_results'
files = glob.glob(os.path.join(folder,'*.pkl'))

file = files[2]

cnm_object = load_object(file)
locals().update(cnm_object)

A, b = cnm2.Ab[:, cnm2.gnb:cnm2.M], cnm2.Ab[:, :cnm2.gnb]
C, f = cnm2.C_on[cnm2.gnb:cnm2.M,:t_cnm], cnm2.C_on[:cnm2.gnb,:t_cnm]
dims = cnm2.dims

#YrA=cnm2.YrA # residual signal
noisyC = cnm2.noisyC[:,:t_cnm]
YrA = noisyC[cnm2.gnb:cnm2.M] - C
Y_r=cnm2.YrA + cnm2.C # raw signal

#%% Visualise OnACID output
view_patches_bar([], A, C, b, f, dims[0], dims[1], YrA=YrA, img=None) #img=Cn for background frame, None for cell mask

#%% Frame-by-frame time
try:
    K = init_com_count
except: pass

if file == files[2]:  # last file: incorrect K value saved
    K = 9
        
tot_onacid_cells = cnm2.N - K
onacid_cells_ix = list(range(K,cnm2.N))

rejected = cnm2.N - len(accepted)
rejected_ix = list(set(range(cnm2.N)) ^ set(accepted))

onacid_accepted_ix = onacid_cells_ix
for ix in rejected_ix:
    if ix in onacid_accepted_ix:
        onacid_accepted_ix.remove(ix)

new_spotted = []
for cell_ix in onacid_accepted_ix:
    trace = C[cell_ix,:t_cnm]
    t_spotted = np.where(trace>0)[0][0] - cnm2.initbatch
    new_spotted.append(t_spotted)
    
new_spotted_rej = []
for cell_ix in rejected_ix:
    trace = C[cell_ix,:t_cnm]
    t_spotted = np.where(trace>0)[0][0] - cnm2.initbatch
    new_spotted_rej.append(t_spotted)

#%% Plot time of each frame to see the delay introduced by various functionalities
tottime = np.array(tottime)
duration = tottime.shape[0] 
shape_ref = cnm2.minibatch_shape

root = tk.Tk()
width = root.winfo_screenwidth()
height = root.winfo_screenheight()

fig = pl.figure(figsize=(width/100., height/100.), dpi=100)
pl.plot(tottime*1000)
pl.plot([0,duration],[33,33],'--m',label = 'Image acquisition time')
pl.plot(np.arange(shape_ref,duration+1,shape_ref)-1, np.ones([int((t_cnm-cnm2.initbatch)/shape_ref),1])*9.5, 'b.', label='Cell shapes refreshed')
pl.plot(new_spotted,np.ones([len(new_spotted),1])*10,'g*', label='New accepted cell detected')
pl.plot(new_spotted_rej,np.ones([len(new_spotted_rej),1])*10,'r*', label='New rejected cell detected')

    
pl.title('Time per frame for the online pipeline  (incorrect accepted cell times)')
pl.xlabel('Frame')
pl.ylabel('Processing time (ms)')
pl.legend()

try:
    figManager = pl.get_current_fig_manager()
    figManager.window.showMaximized()
except:
    mng = pl.get_current_fig_manager()
    mng.frame.Maximize(True)
    
fig.savefig(file[:-4] + '_time_plot.png')
