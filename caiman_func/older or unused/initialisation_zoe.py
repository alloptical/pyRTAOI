import numpy as np
from numpy import linalg as LA

try:
    if __IPYTHON__:
        print('Debugging!')
        # this is used for debugging purposes only. allows to reload classes when changed
        get_ipython().magic('reload_ext autoreload')
        get_ipython().magic('autoreload 2')
except NameError:
    print('Not IPYTHON')
    pass

import time
import sys

if sys.platform == 'win32':
    time_ = time.clock
else:
    time_ = time.time

import caiman as cm
from caiman.source_extraction import cnmf as cnmf
from caiman.utils.visualization import view_patches_bar
from caiman.utils.utils import load_object, save_object
import pylab as pl
import psutil
import scipy
#from caiman.motion_correction import motion_correct_iteration_fast
import cv2
from caiman.utils.visualization import plot_contours
from caiman.source_extraction.cnmf.online_cnmf import bare_initialization
from caiman.base.rois import com
from caiman.utils.visualization import get_contours
#from copy import deepcopy

def initialise(ref_movie, init_method='cnmf', K=3, ds_factor=1, initbatch=300,  # 500?
               T1 = 10000, mot_corr = True,  save_init = False, expected_comps = 100,
               rval_thr = 0.85, NumROIs = None, thresh_overlap = 0.1,
               del_duplicates = False):
    """
    Inputs:
    -------
    
    ref_movie               string
                            reference movie file location
    
    init_method             string
                            choose initialisation method: 'cnmf' (recommended) or 'bare'
    
    K                       int
                            number of initial components to be found
    
    ds_factor               float
                            downsampling factor
    
    T1                      int
                            expected length of the whole recording, in frames (overestimate)
                            
    mot_corr                Boolean
                            flag for online motion correction
                            
    save_init               Boolean
                            flag for saving initialization object
                            
    """
    
    t1 = time_()
    fr = 30                                                             # frame rate (Hz)
    decay_time = 0.5                                                    # approximate length of transient event in seconds
    gSig = (10,10)                                                        # expected half size of neurons
    p = 1                                                               # order of AR indicator dynamics
    min_SNR = 2.5                                                       # minimum SNR for accepting new components
    #rval_thr = 0.85                                                     # correlation threshold for new component inclusion
    #ds_factor = 1                                                       # spatial downsampling factor (increases speed but may lose some fine structure)
    gnb = 1                                                             # number of background components
    gSig = tuple(np.ceil(np.array(gSig)/ds_factor).astype('int'))       # recompute gSig if downsampling is involved
    max_shift = np.ceil(10./ds_factor).astype('int')                    # maximum allowed shift during motion correction
    
    # set up some additional supporting parameters needed for the algorithm (these are default values but change according to dataset characteristics)
    
    # not necessary but useful
    if NumROIs == None:
        max_comp_update_shape = expected_comps
    else:
        max_comp_update_shape = NumROIs
        
    max_comp_update_shape = expected_comps #np.inf                                      # number of shapes to be updated each time (put this to a finite small value to increase speed)

    # init_files = 1                                                      # number of files used for initialization
   # online_files = len(fls) - 1                                         # number of files used for online
    #initbatch = 100                                                     # number of frames for initialization (presumably from the first file)
    #expected_comps = 50                                                # maximum number of expected components used for memory pre-allocation (exaggerate here)
   # K = 2                                                               # initial number of components
    N_samples = np.ceil(fr*decay_time)                                  # number of timesteps to consider when testing new neuron candidates
    thresh_fitness_raw = scipy.special.log_ndtr(-min_SNR)*N_samples     # exceptionality threshold
    #len_file = 5000                                                     # upper bound for number of frames in each file (used right below)
    #T1 = 10000#len(fls)*len_file*epochs                                       # total length of all files (if not known use a large number, then truncate at the end)
    
    
    # Initialize movie
    if ds_factor > 1:                                   # load only the first initbatch frames and possibly downsample them
        Y = cm.load(ref_movie, subindices = slice(0,initbatch,None)).astype(np.float32).resize(1. / ds_factor, 1. / ds_factor)
    else:
        Y = cm.load(ref_movie, subindices = slice(0,initbatch,None)).astype(np.float32)
        
#    if initbatch is None:
#        initbatch = Y.shape[0]
        
    if mot_corr:                                        # perform motion correction on the first initbatch frames
        mc = Y.motion_correct(max_shift, max_shift)
        Y = mc[0].astype(np.float32)
        #borders = np.max(mc[1])
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
    
    
    if init_method == 'bare':
        cnm_init = bare_initialization(Y.transpose(1, 2, 0), init_batch=initbatch, k=K, gnb=gnb,
                                     gSig=gSig, p=p, minibatch_shape=100, minibatch_suff_stat=5,
                                     update_num_comps = True, rval_thr=rval_thr,
                                     thresh_fitness_raw = thresh_fitness_raw,
                                     thresh_overlap = thresh_overlap,  # TODO: check if it works
                                     batch_update_suff_stat=True, max_comp_update_shape = max_comp_update_shape, 
                                     deconv_flag = False, use_dense = True,
                                     simultaneously=False, n_refit=0)
    
    elif init_method == 'cnmf':
        merge_thresh = 0.65
        stride = None   # amount of overlap between patches
        T = initbatch
        n_processes = np.maximum(np.int(psutil.cpu_count()),1)
        images = np.reshape(Yr.T, [T] + list(dims), order='F')
    
        cnm_init = cnmf.CNMF(n_processes, k=K, gSig=gSig, 
                         merge_thresh=merge_thresh,
                         p=p, 
                         #rf=patch_size//2, 
                         #simultaneously = True, #slower but more accurate
                         del_duplicates = True,
                         use_dense = True,
                         thresh_overlap = thresh_overlap,
                         remove_very_bad_comps = True, # apparently can be less precise if true?
                         stride=stride, skip_refinement=False,
                         normalize_init=False, options_local_NMF=None,
                         minibatch_shape=100, minibatch_suff_stat=5,
                         update_num_comps=True, rval_thr=rval_thr,
                         thresh_fitness_delta=-50, gnb = gnb,
                         thresh_fitness_raw=thresh_fitness_raw,
                         batch_update_suff_stat=True, max_comp_update_shape=max_comp_update_shape)
    
        cnm_init = cnm_init.fit(images)
    
    print(('Number of components:' + str(cnm_init.A.shape[-1])))

    # Find centre and radius for cells detected
    coms_init = com(cnm_init.A, d1, d2)
    #coords_init = get_contours(cnm_init.A, dims, thr=0.9) # currently duplicated calculation of coms - NOT NEEDED
    
    # calculate radius as max/mean/median? distance from com to contour --- NOT NEEDED
#    comp_radius = np.zeros((cnm_init.A.shape[-1],1))
#    for index, cell in enumerate(coords_init):
#        cell['coordinates'] = cell['coordinates'][~np.isnan(cell['coordinates']).all(1)]  # remove nan values from contour coordinates
#        cell_coms = np.array([cell['CoM'][1], cell['CoM'][0].copy()])  # reverse x and y for coms
#        diff = cell['coordinates'] - cell_coms
#        dist = LA.norm(diff, axis = 1)
#        comp_radius[index] = int(np.median(dist)) # radius as a max -- kind of big values. may want to use something else, like mean or median?
#        # convert to int?
        



        
    
    if ds_factor == 1:
        last_frame = Y[-1]
    else:
        #h, w = dims
        #f1 = int(512/h)
        #f2 = int(512/h)
        last_frame = cv2.resize(Y[-1], (512, 512), interpolation=cv2.INTER_CUBIC)
        #last_frame = cm.load(ref_movie, subindices = slice(initbatch-1,initbatch,None)).astype(np.float32)


    # Store all relevant intialisation values
    init_values = {}
    
    init_values['mot_corr'] = mot_corr
    init_values['ds_factor'] = ds_factor
    init_values['T1'] = T1
    init_values['N_samples'] = N_samples
    init_values['expected_comps'] = expected_comps
    init_values['K'] = K
    init_values['img_min'] = img_min
    init_values['img_norm'] = img_norm
    init_values['Cn_init'] = Cn_init
    init_values['cnm_init'] = cnm_init
    init_values['coms_init'] = coms_init
    #init_values['coords_init'] = coords_init
#    init_values['comp_radius_init'] = comp_radius
    
    if save_init: #
        cnm_init.dview = None
        init_values['Yr'] = Yr
        save_object(init_values, ref_movie[:-4] + '_DS_' + str(ds_factor) + '.pkl')
#        init_values = load_object(ref_movie[:-4] + '_DS_' + str(ds_factor) + '.pkl')
        
    # Prepare object for OnACID
    cnm_init._prepare_object(np.asarray(Yr), T1, expected_comps, idx_components=None,
                             min_num_trial = 2, N_samples_exceptionality = int(N_samples))
    t2 = time_()
    print('Initialisation time: ' + str(t2-t1)[:6] + ' s')
    
    return last_frame, init_values



if __name__ == '__main__':
    ref_movie = 'C:\\Users\\Patrycja\\CaImAn\\sample_vid\\stimulated_trimmed\\20170329_OG151_t-008_Substack (1-3000)--(1-500).tif'
    lframe, initv = initialise(ref_movie, ds_factor=1.5, initbatch=300)
