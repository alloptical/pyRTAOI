### OFFLINE CNMF ANALYSIS PIPELINE ### 

"""
    Offline analysis can also be implemented using caiman batch i.e. running
    CNMF algorithm on the whole recording. It uses memory mapping to
    efficiently store and access big video files. Motion correction is done in
    advance and is fairly slow.
    
    This generally finds many components, based on K per patch and number of
    patches. If high number of expected components are procided, many of them
    overlap or are not of great quality. Filtering the algorithm output using
    CNN can be adopted; this should significantly reduce the number of 
    components to a reasonable number.
    
    Motion correction and memmaping are slow processes but the rest of the 
    pipeline runs smoothly on large files. The algorithm find similar number
    of cells compared to OnACID run multiple times.
    
"""

from __future__ import division
from __future__ import print_function
from builtins import range

import os
import sys
import cv2
import glob

try:
    cv2.setNumThreads(0)
except:
    pass

try:
    if __IPYTHON__:
        print("Running under iPython")
        # this is used for debugging purposes only. allows to reload classes
        # when changed
        get_ipython().magic('reload_ext autoreload')
        get_ipython().magic('autoreload 2')
except NameError:
    pass

import caiman as cm
import numpy as np
import time
import matplotlib.pyplot as plt

from caiman.utils.visualization import plot_contours, view_patches_bar
from caiman.source_extraction.cnmf import cnmf as cnmf
from caiman.motion_correction import MotionCorrect
from caiman.source_extraction.cnmf.utilities import detrend_df_f
from caiman.components_evaluation import estimate_components_quality_auto

#%% First setup some parameters
#fname = [r'C:\Users\intrinsic\caiman_data\sample_vid\stimulated_test\20170329_OG151_t-008_Substack (1-3000)--(1-500).tif'] # filename to be processed

# example movies
example = 1

sample_folder = r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\samples'
example_folder = os.path.join(sample_folder, 'example' + str(example))
files = glob.glob(os.path.join(example_folder,'*Ch2.tif'))

if len(files)==1:
    fname = files
print('Ref movie path: ', fname)

# dataset dependent parameters
display_images = False              # Set this to true to show movies and plots

#fname = ['Sue_2x_3000_40_-46.tif'] 

fr = 30                             # imaging rate in frames per second
decay_time = 0.2                    # length of a typical transient in seconds


######## motion correction parameters
niter_rig = 1               # number of iterations for rigid motion correction
max_shifts = (10, 10)   # (6, 6)      # maximum allow rigid shift
# for parallelization split the movies in  num_splits chuncks across time
splits_rig = 100 # 56
# start a new patch for pw-rigid motion correction every x pixels
strides = (144, 144) #(146, 146) # (48, 48)
# overlap between pathes (size of patch strides+overlaps)
overlaps = (24, 24) # (73, 73) # (24, 24)
# for parallelization split the movies in  num_splits chuncks across time
splits_els = 100 # 56
upsample_factor_grid = 4    # upsample factor to avoid smearing when merging patches
# maximum deviation allowed for patch with respect to rigid shifts
max_deviation_rigid = 3


######## parameters for source extraction and deconvolution
p = 1                       # order of the autoregressive system
gnb = 1                     # number of global background components
merge_thresh = 0.8          # merging threshold, max correlation allowed
# half-size of the patches in pixels. e.g., if rf=25, patches are 50x50
rf = 50 #15
stride_cnmf = 5 #6             # amount of overlap between the patches in pixels
K = 10 #3                       # number of components per patch
gSig = [10, 10]               # expected half size of neurons
# initialization method (if analyzing dendritic data using 'sparse_nmf')
init_method = 'greedy_roi'
is_dendrites = False        # flag for analyzing dendritic data
# sparsity penalty for dendritic data analysis through sparse NMF
alpha_snmf = None

# parameters for component evaluation
min_SNR = 2.5              # signal to noise ratio for accepting a component
rval_thr = 0.75              # space correlation threshold for accepting a component
cnn_thr = 0.8               # threshold for CNN based classifier
use_cnn = True

#%%
#import tifffile
#M = tifffile.TiffFile(fname[0])
#print(M.is_imagej)

#%% play the movie
# playing the movie using opencv. It requires loading the movie in memory. To
# close the video press q

if display_images:
    m_orig = cm.load_movie_chain(fname[:1])
    downsample_ratio = 0.2
    offset_mov = -np.min(m_orig[:100])
    moviehandle = m_orig.resize(1, 1, downsample_ratio)
    moviehandle.play(gain=10, offset=offset_mov, fr=30, magnification=2)

#%% start a cluster for parallel processing
try:
    dview.terminate()
#     cm.stop_server(dview=dview)
except:
    print('OK')

c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=None, single_thread=False)


#%%% MOTION CORRECTION
# first we create a motion correction object with the parameters specified
min_mov = cm.load(fname[0], subindices=range(200)).min()
# this will be subtracted from the movie to make it non-negative

mc = MotionCorrect(fname[0], min_mov,
                   dview=dview,  # None for now - then it runs
                   max_shifts=max_shifts, niter_rig=niter_rig,
                   splits_rig=splits_rig,
                   strides=strides, overlaps=overlaps, splits_els=splits_els,
                   upsample_factor_grid=upsample_factor_grid,
                   max_deviation_rigid=max_deviation_rigid,
                   shifts_opencv=True, nonneg_movie=True, use_cuda=False,
                   border_nan='copy')

# note that the file is not loaded in memory

# currently error below when use_cuda=True:
#  File "C:\ProgramData\Anaconda3\envs\caiman\lib\multiprocessing\pool.py", line 644, in get
#    raise self._value
#
#CompileError: nvcc preprocessing of C:\Users\INTRIN~1\AppData\Local\Temp\tmpwan_5_ki.cu failed
#%% Run piecewise-rigid motion correction using NoRMCorre
#mc.motion_correct_rigid(save_movie=True)
mc.motion_correct_pwrigid(save_movie=True)

m_els = cm.load(mc.fname_tot_els)
bord_px_els = np.ceil(np.maximum(np.max(np.abs(mc.x_shifts_els)),
                                 np.max(np.abs(mc.y_shifts_els)))).astype(np.int)
# maximum shift to be used for trimming against NaNs

# save bord_px_els
folder = os.path.dirname(fname[0])
save = os.path.join(folder, 'bord_px_els.npz')
np.savez(save, bord_px_els=bord_px_els)

#%% compare with original movie
display_images = False
if display_images:
    moviehandle = cm.concatenate([m_orig.resize(1, 1, downsample_ratio) + offset_mov,
                m_els.resize(1, 1, downsample_ratio)],
               axis=2)
    moviehandle.play(fr=60, q_max=99.5, magnification=2, offset=0)  # press q to exit

#%% memory map the file in order 'C'
fnames = mc.fname_tot_els   # name of the pw-rigidly corrected file.
border_to_0 = bord_px_els     # number of pixels to exclude
fname_new = cm.save_memmap(fnames, base_name='memmap_', order='C',
                           border_to_0=bord_px_els)  # exclude borders

#%% now load the file
#fname_new = r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\samples\example8\memmap__d1_512_d2_512_d3_1_order_C_frames_4956_.mmap'

fname_new = r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\samples\example1\memmap__d1_512_d2_512_d3_1_order_C_frames_5400_.mmap'
b = np.load(r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\samples\example1\bord_px_els.npz')
bord_px_els = b['bord_px_els']

Yr, dims, T = cm.load_memmap(fname_new)
d1, d2 = dims
images = np.reshape(Yr.T, [T] + list(dims), order='F')
# load frames in python format (T x X x Y)

# load board_px_els
#d = np.load(r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\samples\example8\bord_px_els.npz')
#try:
#    bord_px_els = bord_px_els
#except:
#    bord_px_els = d['bord_px_els']
    
# restart cluster to clean up memory
try:
    dview.terminate()
except: pass
c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=None, single_thread=False)

#%% Normalise data
#img_min = Y.min()                                   # minimum value of movie. Subtract it to make the data non-negative
#Y -= img_min
#img_norm = np.std(Y, axis=0)                        
#img_norm += np.median(img_norm)                     # normalizing factor to equalize the FOV
#Y = Y / img_norm[None, :, :]                        # normalize data
#
#_, d1, d2 = Y.shape
#dims = (d1, d2)                                     # dimensions of FOV
#Yr = Y.to_2D().T                                    # convert data into 2D array                                    
#    
#images = np.reshape(Yr.T, [T] + list(dims), order='F')

#%% RUN CNMF ON PATCHES

# First extract spatial and temporal components on patches and combine them
# for this step deconvolution is turned off (p=0)
t1 = time.time()

cnm = cnmf.CNMF(n_processes=1,
                k=20, 
                gSig=gSig, merge_thresh=merge_thresh,
                p=0, 
                dview=dview, 
                rf=rf, 
                stride=stride_cnmf,
                memory_fact=0.5,
                method_init=init_method, alpha_snmf=alpha_snmf,
                only_init_patch=False, gnb=gnb, border_pix=bord_px_els)
cnm = cnm.fit(images)

#%% plot contours of found components
#Cn = cm.local_correlations(images.transpose(1, 2, 0))
#Cn[np.isnan(Cn)] = 0
#plt.figure()
#crd = plot_contours(cnm.A, Cn, thr=0.9)
#plt.title('Contour plots of found components')


#%% COMPONENT EVALUATION
# the components are evaluated in three ways:
#   a) the shape of each component must be correlated with the data
#   b) a minimum peak SNR is required over the length of a transient
#   c) each shape passes a CNN based classifier

idx_components, idx_components_bad, SNR_comp, r_values, cnn_preds = \
    estimate_components_quality_auto(images, cnm.A, cnm.C, cnm.b, cnm.f,
                                     cnm.YrA, fr, decay_time, gSig, dims,
                                     dview=dview, min_SNR=1.5,#min_SNR,
                                     r_values_min=rval_thr, use_cnn=use_cnn,
                                     thresh_cnn_min=cnn_thr)

print('Number of components accepted: ', len(idx_components))

#%% PLOT COMPONENTS
try:
    Cn
except:
    Cn = cm.local_correlations(images[:100,:,:].transpose(1, 2, 0))   # approx img of first frames for quicker display
    Cn[np.isnan(Cn)] = 0

if display_images:
    plt.figure()
    plt.subplot(121)
    crd_good = cm.utils.visualization.plot_contours(
        cnm.A[:, idx_components], Cn, thr=.8, vmax=0.75)
    plt.title('Contour plots of accepted components')
    plt.subplot(122)
    crd_bad = cm.utils.visualization.plot_contours(
        cnm.A[:, idx_components_bad], Cn, thr=.8, vmax=0.75)
    plt.title('Contour plots of rejected components')

#%% VIEW TRACES (accepted and rejected)

if display_images:
    view_patches_bar(Yr, cnm.A.tocsc()[:, idx_components], cnm.C[idx_components],
                     cnm.b, cnm.f, dims[0], dims[1], YrA=cnm.YrA[idx_components],
                     img=Cn)

    view_patches_bar(Yr, cnm.A.tocsc()[:, idx_components_bad], cnm.C[idx_components_bad],
                     cnm.b, cnm.f, dims[0], dims[1], YrA=cnm.YrA[idx_components_bad],
                     img=Cn)
    
#%% View only accepted
view_patches_bar(Yr, cnm.A.tocsc()[:, idx_components], cnm.C[idx_components],
             cnm.b, cnm.f, dims[0], dims[1], YrA=cnm.YrA[idx_components],
             img=None)

#%% RE-RUN seeded CNMF on accepted patches to refine and perform deconvolution
A_in, C_in, b_in, f_in = cnm.A[:,
                               idx_components], cnm.C[idx_components], cnm.b, cnm.f
cnm2 = cnmf.CNMF(n_processes=1, k=A_in.shape[-1], gSig=gSig, p=p, dview=dview,
                 merge_thresh=merge_thresh, Ain=A_in, Cin=C_in, b_in=b_in,
                 f_in=f_in, rf=None, stride=None, gnb=gnb,
                 method_deconvolution='oasis', check_nan=True)

cnm2 = cnm2.fit(images)

#%% SIMPLE CNMF
simple = 0

if simple:
    thresh_overlap = 0.2
    K = 100 # tot K expected 
    
    cnm2 = cnmf.CNMF(n_processes=1, k=K, gSig=gSig,
                     merge_thresh=merge_thresh, p=p,
                     #rf=patch_size//2,
                     stride=None, rf=None,
                     simultaneously=False,
                     del_duplicates=True, 
                     use_dense=True,
                     thresh_overlap=thresh_overlap,
                     remove_very_bad_comps=True,
                     skip_refinement=False,
                     normalize_init=False, options_local_NMF=None,
    #                 minibatch_shape=minibatch_shape, minibatch_suff_stat=5,
                     update_num_comps=True, rval_thr=rval_thr,
                     thresh_fitness_delta=-50, gnb=gnb,
    #                 thresh_fitness_raw=thresh_fitness_raw,
                     batch_update_suff_stat=True)#, max_comp_update_shape=max_comp_update_shape)
    
    cnm2 = cnm2.fit(images)

#%% Detrend DF/F values without loading the movie 

#    F_df:
#        the computed Calcium activity to the derivative of f

F_dff = detrend_df_f(cnm2.A, cnm2.b, cnm2.C, cnm2.f, YrA=cnm2.YrA,
                     quantileMin=8, frames_window=250)

F_dff_no_noise = detrend_df_f(cnm2.A, cnm2.b, cnm2.C, cnm2.f,
                     quantileMin=8, frames_window=250)

#%% Display DF/F -- very noisy if YrA provided
cell = 0
plt.figure()
plt.plot(F_dff[cell,:])

#%% Extract DF/F values -- memory inefficient (a bit) but traces a lot less noisy

#    Cdf:
#        the computed Calcium activity to the derivative of f
        
from caiman.source_extraction.cnmf.utilities import extract_DF_F

bl = cnm2.bl # baseline for each component

C_df = extract_DF_F(Yr, cnm2.A, cnm2.C, bl)

#%% Display C_df
cell = 0
plt.figure()
plt.plot(C_df[cell,:])


#%% Show final traces
cnm2.view_patches(Yr, dims=dims, img=None) #=Cn

#%% Extract coms of cells found
from caiman.base.rois import com

coms = com(cnm2.A, *dims)

#%% Save results -- NotImplementedError: pool objects cannot be passed between processes or pickled
#from caiman.utils.utils import save_object
#
#save_dict = dict()
#
#save_dict['cnm2'] = cnm2
##save_dict['F_dff'] = F_dff
#save_dict['coms'] = coms
#
#folder = os.path.join(os.path.dirname(fname[0]), 'offline_results')
#saveResultPath = os.path.join(folder, 'ex1_cnmf_results.pkl')
#
#save_object(save_dict, saveResultPath)
#
#save_mat = False # additionally save the mat file
#if save_mat:
#    from pkl2mat import convert2mat
#    convert2mat(file_full_name = saveResultPath)

#%% Load npz file
file = r'C:\Users\intrinsic\Desktop\pyRTAOI-rig\samples\example1\offline_results\ex3_cnmf_results.npz'
file_data = np.load(file)
locals().update(file_data)
#
#cnm2 = cnmf.CNMF
#cnm2.A = A
#cnm2.b = b
#cnm2.C = C
#cnm2.f = f
#cnm2.YrA = YrA

#%% Saving as npz works
folder = os.path.dirname(fname[0])
results_folder = os.path.join(folder, 'offline_results')

if not os.path.exists(results_folder):
    os.makedirs(results_folder)
    
saved_data = os.path.join(results_folder, 'ex' + str(example) + '_cnmf_results_minSRN1.5.npz')

np.savez(saved_data, A=cnm2.A, b=cnm2.b, C=cnm2.C, f=cnm2.f, YrA=cnm2.YrA, # YrA is residual signal
         Y_r=cnm2.YrA+cnm2.C, Cn=Cn, dims=cnm2.dims, 
#         F_dff=F_dff, 
#         F_dff_no_noise=F_dff_no_noise,
         C_df=C_df, coms=coms, gnb=gnb,
         thresh_overlap=cnm2.thresh_overlap, cnm_N=cnm2.A.shape[1],
         merge_thresh=merge_thresh, gSig=gSig, min_SNR=min_SNR,
         rval_thr=rval_thr, cnn_thr=cnn_thr, use_cnn=use_cnn, rf=rf, K=K)

# Save npz to mat
from npz2mat import npz2mat

save_mat = 1
if save_mat:
    npz2mat(saved_data)

#%% STOP CLUSTER and clean up log files
cm.stop_server(dview=dview)
log_files = glob.glob('*_LOG_*')
for log_file in log_files:
    os.remove(log_file)
    
#%% reconstruct denoised movie
denoised = cm.movie(cnm2.A.dot(cnm2.C) +
                    cnm2.b.dot(cnm2.f)).reshape(dims + (-1,), order='F').transpose([2, 0, 1])

#%% play along side original data
moviehandle = cm.concatenate([m_els.resize(1, 1, downsample_ratio),
                denoised.resize(1, 1, downsample_ratio)],
               axis=2)
if display_images:
        moviehandle.play(fr=60, gain=15, magnification=2, offset=0)  # press q to exit
