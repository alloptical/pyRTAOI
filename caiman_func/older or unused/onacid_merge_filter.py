# -*- coding: utf-8 -*-
"""
Created on Sun Dec  3 12:40:25 2017

@author: Patrycja
"""

import numpy as np
try:
    if __IPYTHON__:
        print('Debugging!')
        # this is used for debugging purposes only. allows to reload classes when changed
        get_ipython().magic('load_ext autoreload')
        get_ipython().magic('autoreload 2')
except NameError:
    print('Not IPYTHON')
    pass

from time import time
import caiman as cm
from caiman.utils.visualization import view_patches_bar
from caiman.utils.utils import download_demo, load_object, save_object
import pylab as pl
import scipy
from caiman.motion_correction import motion_correct_iteration_fast
import cv2
from caiman.utils.visualization import plot_contours
import glob
from caiman.source_extraction.cnmf.online_cnmf import bare_initialization
from copy import deepcopy

#%%
folder_name = 'C:/Users/Patrycja/CaImAn/sample_vid/spontaneous_cropped'                       # folder where files are located
extension = '.tif'                                  # extension of files
fls = glob.glob(folder_name+'/*' + extension)       # read all files to be processed 
print(fls)

#%%
import os
os.chdir(folder_name)
os.getcwd()

#%%   Set up some parameters
fr = 30                                                             # frame rate (Hz)
decay_time = 0.5                                                    # approximate length of transient event in seconds
gSig = (10,10)                                                      # expected half size of neurons
p = 1                                                               # order of AR indicator dynamics
min_SNR = 2.5                                                       # minimum SNR for accepting new components
rval_thr = 0.7                                                      # correlation threshold for new component inclusion
ds_factor = 1                                                       # spatial downsampling factor (increases speed but may lose some fine structure)
gnb = 1                                                             # number of background components
gSig = tuple(np.ceil(np.array(gSig)/ds_factor).astype('int'))       # recompute gSig if downsampling is involved
mot_corr = True                                                     # flag for online motion correction 
max_shift = np.ceil(10./ds_factor).astype('int')                    # maximum allowed shift during motion correction

# set up some additional supporting parameters needed for the algorithm (these are default values but change according to dataset characteristics)

max_comp_update_shape = np.inf                                      # number of shapes to be updated each time (put this to a finite small value to increase speed)
init_files = 1                                                      # number of files used for initialization
online_files = len(fls) - 1                                         # number of files used for online
initbatch = 200                                                     # number of frames for initialization (presumably from the first file)
expected_comps = 50                                                # maximum number of expected components used for memory pre-allocation (exaggerate here)
K = 2                                                               # initial number of components
N_samples = np.ceil(fr*decay_time)                                  # number of timesteps to consider when testing new neuron candidates
thresh_fitness_raw = scipy.special.log_ndtr(-min_SNR)*N_samples     # exceptionality threshold
epochs = 1                                                          # number of passes over the data
len_file = 3000                                                     # upper bound for number of frames in each file (used right below)
T1 = len(fls)*len_file*epochs                                       # total length of all files (if not known use a large number, then truncate at the end)
merge_thresh = 0.7

#%% try cnmf/other on c1v1?



#%%    Initialize movie
if ds_factor > 1:                                   # load only the first initbatch frames and possibly downsample them
    Y = cm.load(fls[0], subindices = slice(0,initbatch,None)).astype(np.float32).resize(1. / ds_factor, 1. / ds_factor)
else:
    Y =  cm.load(fls[0], subindices = slice(0,initbatch,None)).astype(np.float32)
    
if mot_corr:                                        # perform motion correction on the first initbatch frames
    mc = Y.motion_correct(max_shift, max_shift)
    Y = mc[0].astype(np.float32)
    borders = np.max(mc[1])
else:
    Y = Y.astype(np.float32)
      
img_min = Y.min()                                   # minimum value of movie. Subtract it to make the data non-negative
Y -= img_min
img_norm = np.std(Y, axis=0)                        
img_norm += np.median(img_norm)                     # normalizing factor to equalize the FOV
Y = Y / img_norm[None, :, :]                        # normalize data

_, d1, d2 = Y.shape
dims = (d1, d2)                                     # dimensions of FOV
Yr = Y.to_2D().T                                    # convert data into 2D array                                    

Cn_init = Y.local_correlations(swap_dim = False)    # compute correlation image
pl.imshow(Cn_init); pl.title('Correlation Image on initial batch'); pl.colorbar()

#%% initialize OnACID with bare initialization
cnm_init = bare_initialization(Y[:initbatch].transpose(1, 2, 0), init_batch=initbatch, k=K, gnb=gnb,
                                 gSig=gSig, p=p, minibatch_shape=100, minibatch_suff_stat=5,
                                 method_init = 'greedy_roi',  #doesn't work with sparse_nmf? 
                                 #options: greedy_roi, sparse_nmf, local_nmf, pca_ica, corr_pnr 
                                 merge_thresh = merge_thresh, remove_very_bad_comps=True,  #added/changed
                                 update_num_comps = True, rval_thr=rval_thr,
                                 thresh_fitness_raw = thresh_fitness_raw,
                                 batch_update_suff_stat=True, max_comp_update_shape = max_comp_update_shape, 
                                 deconv_flag = False, use_dense = True,
                                 simultaneously=False, n_refit=0)

#%% Plot initialization results
pl.figure()
crd = plot_contours(cnm_init.A.tocsc(), Cn_init, thr=0.9)
A, C, b, f, YrA, sn = cnm_init.A, cnm_init.C, cnm_init.b, cnm_init.f, cnm_init.YrA, cnm_init.sn
view_patches_bar(Yr, scipy.sparse.coo_matrix(A.tocsc()[:, :]), C[:, :], b, f, dims[0], dims[1], YrA=YrA[:, :], img=Cn_init)

#%% Prepare object for OnACID
save_init = False     # flag for saving initialization object. Useful if you want to check OnACID with different parameters but same initialization
if save_init:   
    cnm_init.dview = None
    save_object(cnm_init, fls[0][:-4] + '_DS_' + str(ds_factor) + '.pkl')
    cnm_init = load_object(fls[0][:-4] + '_DS_' + str(ds_factor) + '.pkl')
    
cnm_init._prepare_object(np.asarray(Yr), T1, expected_comps, idx_components=None, 
                         min_num_trial = 2, N_samples_exceptionality = int(N_samples))


#%% create a function for plotting results in real time if needed
def create_frame(cnm2,img_norm,captions):
    A, b = cnm2.Ab[:, cnm2.gnb:], cnm2.Ab[:, :cnm2.gnb].toarray()
    C, f = cnm2.C_on[cnm2.gnb:cnm2.M, :], cnm2.C_on[:cnm2.gnb, :]
    comps_frame = A.dot(C[:,t-1]).reshape(cnm2.dims, order = 'F')*img_norm/np.max(img_norm)   # inferred activity due to components (no background)
    bgkrnd_frame = b.dot(f[:,t-1]).reshape(cnm2.dims, order = 'F')*img_norm/np.max(img_norm)  # denoised frame (components + background)
    if show_residuals:
        all_comps = np.reshape(cnm2.Yres_buf.mean(0),cnm2.dims, order='F')*img_norm/np.max(img_norm)
        all_comps = np.minimum(np.maximum(all_comps*10,0),255)
    else:
        all_comps = (np.array(A.sum(-1)).reshape(cnm2.dims, order = 'F'))                         # spatial shapes
    frame_comp_1 = cv2.resize(np.concatenate([frame_/np.max(img_norm),all_comps*3.],axis = -1),(2*np.int(cnm2.dims[1]*resize_fact),np.int(cnm2.dims[0]*resize_fact) ))
    frame_comp_2 = cv2.resize(np.concatenate([comps_frame*10.,comps_frame+bgkrnd_frame],axis = -1),(2*np.int(cnm2.dims[1]*resize_fact),np.int(cnm2.dims[0]*resize_fact) ))
    frame_pn = np.concatenate([frame_comp_1,frame_comp_2],axis=0).T
    vid_frame = np.repeat(frame_pn[:,:,None],3,axis=-1)
    vid_frame = np.minimum((vid_frame*255.),255).astype('u1')
    cv2.putText(vid_frame,captions[0],(5,20),fontFace = 5, fontScale = 1.2, color = (0,255,0), thickness = 1)
    cv2.putText(vid_frame,captions[1],(np.int(cnm2.dims[0]*resize_fact) + 5,20),fontFace = 5, fontScale = 1.2, color = (0,255,0), thickness = 1)
    cv2.putText(vid_frame,captions[2],(5,np.int(cnm2.dims[1]*resize_fact)  + 20),fontFace = 5, fontScale = 1.2, color = (0,255,0), thickness = 1)
    cv2.putText(vid_frame,captions[3],(np.int(cnm2.dims[0]*resize_fact) + 5 ,np.int(cnm2.dims[1]*resize_fact)  + 20),fontFace = 5, fontScale = 1.2, color = (0,255,0), thickness = 1)
    cv2.putText(vid_frame,'Frame = '+str(t),(vid_frame.shape[1]//2-vid_frame.shape[1]//10,vid_frame.shape[0]-20),fontFace = 5, fontScale = 1.2, color = (0,255,255), thickness = 1)      
    return vid_frame

#%%
YYY = np.zeros([3000,384,432])#dims[0],dims[1]])   #memory error. seemed to work before though?

#YYY[:initbatch,:,:] = Y
#%% Run OnACID and optionally plot results in real time
cnm2 = deepcopy(cnm_init)
cnm2.Ab_epoch = []                       # save the shapes at the end of each epoch
t = cnm2.initbatch                       # current timestep
tottime = []
Cn = Cn_init.copy()

#*************** SOMETHING OFF RN SO LEAVE PLOT_CONTOURS_FLAG FALSE *********************
plot_contours_flag = False               # flag for plotting contours of detected components at the end of each file
play_reconstr = False                    # flag for showing video with results online (turn off flags for improving speed)
save_movie = False                       # flag for saving movie (file could be quite large..)
movie_name = folder_name + '/output.avi' # name of movie to be saved
resize_fact = 0.8                        # image resizing factor

if online_files == 0:                    # check whether there are any additional files
    process_files = fls[:init_files]     # end processing at this file
    init_batc_iter = [initbatch]         # place where to start
    end_batch = T1              
else:
    process_files = fls[:init_files + online_files]     # additional files
    init_batc_iter = [initbatch] + [0]*online_files     # where to start reading at each file


shifts = []
show_residuals = False
if show_residuals:
    caption = 'Mean Residual Bufer'
else:
    caption = 'Identified Components'
captions = ['Raw Data','Inferred Activity',caption,'Denoised Data']
if save_movie and play_reconstr:
    #fourcc = cv2.VideoWriter_fourcc('8', 'B', 'P', 'S') 
    fourcc = cv2.VideoWriter_fourcc(*'XVID')
    out = cv2.VideoWriter(movie_name,fourcc, 30.0, tuple([int(2*x*resize_fact) for x in cnm2.dims]))

for iter in range(epochs):    
    if iter > 0:
        process_files = fls[:init_files + online_files]     # if not on first epoch process all files from scratch
        init_batc_iter = [0]*(online_files+init_files)      
        
    for file_count, ffll in enumerate(process_files):  # np.array(fls)[np.array([1,2,3,4,5,-5,-4,-3,-2,-1])]:
        print('Now processing file ' + ffll)
        Y_ = cm.load(ffll, subindices=slice(init_batc_iter[file_count],T1,None))
        
        if plot_contours_flag:   # update max-correlation (and perform offline motion correction) just for illustration purposes
            if ds_factor > 1:
                Y_1 = Y_.resize(1. / ds_factor, 1. / ds_factor, 1)
            else:
                Y_1 = Y_.copy()                    
            if mot_corr:
                templ = (cnm2.Ab.data[:cnm2.Ab.indptr[1]] * cnm2.C_on[0, t - 1]).reshape(cnm2.dims, order='F') * img_norm        
                newcn = (Y_1 - img_min).motion_correct(max_shift, max_shift, template=templ)[0].local_correlations(swap_dim=False)                
                Cn = np.maximum(Cn, newcn)
            else:
                Cn = np.maximum(Cn, Y_1.local_correlations(swap_dim=False))
    
        old_comps = cnm2.N                              # number of existing components
        i = 0
        for frame_count, frame in enumerate(Y_):        # now process each file
            if np.isnan(np.sum(frame)):
                raise Exception('Frame ' + str(frame_count) + ' contains nan')
            if t % 100 == 0:
                print('Epoch: ' + str(iter+1) + '. ' + str(t)+' frames have beeen processed in total. '+str(cnm2.N - old_comps)+' new components were added. Total number of components is '+str(cnm2.Ab.shape[-1]-gnb))
                old_comps = cnm2.N
    
            t1 = time()                                 # count time only for the processing part
            frame_ = frame.copy().astype(np.float32)    # 
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
            tottime.append(time() - t1)                             # store time
    
            #print(frame_cor.shape) -- (512,512)
            ### works but super slow:
            curr = np.array(frame_cor)
            YYY[initbatch+i,:,:] = curr[np.newaxis,:,:]
            i += 1
            #YYY = np.concatenate([YYY,curr[np.newaxis,:,:]])
    
            t += 1
            
            if t % 1000 == 0 and plot_contours_flag:
                pl.cla()
                A = cnm2.Ab[:, cnm2.gnb:]
                crd = cm.utils.visualization.plot_contours(A, Cn, thr=0.9)  # update the contour plot every 1000 frames
                pl.pause(1)
                
            if play_reconstr:                                               # generate movie with the results
                vid_frame = create_frame(cnm2,img_norm,captions)             
                if save_movie:
                    out.write(vid_frame)
                cv2.imshow('frame',vid_frame)                
                if cv2.waitKey(1) & 0xFF == ord('q'):
                    break                                
    
        print('Cumulative processing speed is ' + str((t - initbatch) / np.sum(tottime))[:5] + ' frames per second.')
    cnm2.Ab_epoch.append(cnm2.Ab.copy())                        # save the shapes at the end of each epoch
        
if save_movie:
    out.release()
cv2.destroyAllWindows()
#%%  save results (optional)
save_results = False

if save_results:
    np.savez('results_analysis_online_MOT_CORR.npz',
             Cn=Cn, Ab=cnm2.Ab, Cf=cnm2.C_on, b=cnm2.b, f=cnm2.f,
             dims=cnm2.dims, tottime=tottime, noisyC=cnm2.noisyC, shifts=shifts)

#%% extract results from the objects and do some plotting
A, b = cnm2.Ab[:, cnm2.gnb:], cnm2.Ab[:, :cnm2.gnb].toarray()
C, f = cnm2.C_on[cnm2.gnb:cnm2.M, t-t//epochs:t], cnm2.C_on[:cnm2.gnb, t-t//epochs:t]
noisyC = cnm2.noisyC[:,t-t//epochs:t]
b_trace = [osi.b for osi in cnm2.OASISinstances] if hasattr(cnm2, 'OASISinstances') else [0]*C.shape[0]

pl.figure()
crd = cm.utils.visualization.plot_contours(A, Cn, thr=0.9)
view_patches_bar(Yr, scipy.sparse.coo_matrix(A.tocsc()[:, :]), C[:, :], b, f,
                 dims[0], dims[1], YrA=noisyC[cnm2.gnb:cnm2.M] - C, img=Cn)

#%% ******************** 52 COMPONENTS DETECTED **************************

#%% ****************** MERGING CELLS AND FILTERING NOISY COMPONENTS**********************
import os
import psutil
import caiman.source_extraction.cnmf as cnmf
from past.utils import old_div
from caiman.components_evaluation import evaluate_components

#%%
backend='local'
if backend == 'SLURM':
    n_processes = np.int(os.environ.get('SLURM_NPROCS'))
else:
    n_processes = np.maximum(np.int(psutil.cpu_count()),1) # roughly number of cores on your machine minus 1
print(('using ' + str(n_processes) + ' processes'))
#%% set parameters for full field of view analysis

#if type(Y) is tuple:
#    dims, T = Y[:-1], Y[-1]
merge_thresh = 0.75  # vary this

dim1 = dims[0]
dim2 = dims[1]

fY = (dim1, dim2, T1) # fake Y

options = cnmf.utilities.CNMFSetParms(fY,n_processes,p=1,gSig=gSig,K=A.shape[-1],thr=merge_thresh)
T = T1 #length of all files

pix_proc=np.minimum(np.int((d1*d2)/n_processes/(old_div(T,2000.))),np.int(old_div((d1*d2),n_processes))) # regulates the amount of memory used
options['spatial_params']['n_pixels_per_process']=pix_proc
options['temporal_params']['n_pixels_per_process']=pix_proc

#%% NEED YR OF THE WHOLE VIDEO I THINK FOR THIS. HMMMMMM
#Y_tot = cm.load(fls[0]).astype(np.float32)       #something wrong with Y in denoising section

Y_tot = cm.load(ref_movie_path, subindices = slice(0,500,None)).astype(np.float32).resize(1. / ds_factor, 1. / ds_factor)

#TODO: probably don't want it over the whole video: but then it is not included in merging at all

img_min = Y_tot.min()                                   # minimum value of movie. Subtract it to make the data non-negative
Y_tot -= img_min
img_norm = np.std(Y_tot, axis=0)                        
img_norm += np.median(img_norm)                     # normalizing factor to equalize the FOV
Y_tot = Y_tot / img_norm[None, :, :]                        # normalize data

_, d1, d2 = Y_tot.shape
dims = (d1, d2)                                     # dimensions of FOV

Yr_tot = Y_tot.to_2D().T
#Y_tot = np.swapaxes(Y_tot,1,2).T.shape
Y_tot = Y_tot.transpose(1,2,0) # apparently not the same as Y from caiman_patches -- but why??
#TODO: maybe figure this out 
#%% merge spatially overlaping and temporally correlated components   
dview = None #depends on clusters enabled  
t = 500

#Yr_tot = np.zeros([dims[0]*dims[1],t]) # [d, t] = np.shape(Y)

A_m,C_m,nr_m,merged_ROIs,S_m,bl_m,c1_m,sn_m,g_m=cnmf.merging.merge_components(
        Yr_tot,A,[],np.array(C[:,:t]),[],np.array(C[:,:t]),[],options['temporal_params'],options['spatial_params'],
        dview=dview,thr=options['merging']['thr'],mx=np.Inf,fast_merge=False)

#%% ******************* COMPONENTS REDUCED TO 37 ***********************
#%%#%% update temporal to get Y_r
options['temporal_params']['p']=1#0
options['temporal_params']['fudge_factor']=0.96 #change if denoised traces time constant is wrong
options['temporal_params']['backend']='ipyparallel'
options['temporal_params']['block_size'] = 5000 #ADDED
C_m,A_m,b,f_m,S_m,bl_m,c1_m,neurons_sn_m,g2_m,YrA_m,lam = cnmf.temporal.update_temporal_components(Yr_tot,A_m,b,C_m,f[:,:t],dview=dview,bl=None,c1=None,sn=None,g=None,**options['temporal_params'])
 
#%% plotting
pl.figure()
crd = cm.utils.visualization.plot_contours(A_m, Cn, thr=0.9)
view_patches_bar(Yr_tot, scipy.sparse.coo_matrix(A_m.tocsc()[:, :]), C_m[:, :], b, f_m,
                 dims[0], dims[1], YrA=YrA_m, img=Cn)

#%% get rid of evenrually noisy components. 
# But check by visual inspection to have a feeling fot the threshold. Try to be loose, you will be able to get rid of more of them later!
final_frate = 30
Npeaks=10
traces=C_m+YrA_m
#        traces_a=traces-scipy.ndimage.percentile_filter(traces,8,size=[1,np.shape(traces)[-1]/5])
#        traces_b=np.diff(traces,axis=1)

fitness_raw, fitness_delta, erfc_raw, erfc_delta, r_values, significant_samples =\
             evaluate_components(Y_tot, traces, A_m, C_m, b, f_m, \
             remove_baseline=True, N=5, robust_std=False, Athresh = 0.1, Npeaks = Npeaks, thresh_C = 0.3,
             final_frate = final_frate) # tB=tB, tA = tA,

idx_components_r=np.where(r_values>=.5)[0]
idx_components_raw=np.where(fitness_raw<-20)[0]        
idx_components_delta=np.where(fitness_delta<-10)[0]   


idx_components=np.union1d(idx_components_r,idx_components_raw)
idx_components=np.union1d(idx_components,idx_components_delta)  
idx_components_bad=np.setdiff1d(list(range(len(traces))),idx_components)

print(' ***** ')
print((len(traces)))
print((len(idx_components)))

#%%
A_m=A_m.tocsc()[:,idx_components]
C_m=C_m[idx_components,:]
YrA_m = YrA_m[idx_components,:]

#%% display components  DO NOT RUN IF YOU HAVE TOO MANY COMPONENTS
pl.figure()
crd = plot_contours(A_m,Cn,thr=0.9)
print(('Number of components:' + str(A_m.shape[-1]))) 

#%% plotting
pl.figure()
crd = cm.utils.visualization.plot_contours(A_m, Cn, thr=0.9)
view_patches_bar(Yr, scipy.sparse.coo_matrix(A_m.tocsc()[:, :]), C_m[:, :], b, f_m,
                 dims[0], dims[1], YrA=YrA_m, img=Cn)


#%% TRY UPDATING SPATIAL COMPONENTS FUNCTION TO REMOVE EMPTY COMPONENTS -- HOW ARE THESE DEFINED?

