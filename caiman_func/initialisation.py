import numpy as np

try:
    if __IPYTHON__:
        print('Debugging!')
        # this is used for debugging purposes only. allows to reload classes when changed
        get_ipython().magic('reload_ext autoreload')
        get_ipython().magic('autoreload 2')
except NameError:
    print('Not IPYTHON')
    pass

import os
import time
import sys
from copy import deepcopy

if sys.platform == 'win32':
    time_ = time.clock
else:
    time_ = time.time

import psutil
import scipy
import cv2

import caiman as cm
from caiman.source_extraction import cnmf as cnmf
from caiman.utils.utils import load_object, save_object
from caiman.source_extraction.cnmf.online_cnmf import bare_initialization, seeded_initialization
from caiman.base.rois import com
#from caiman.paths import caiman_datadir
from caiman.components_evaluation import evaluate_components_CNN



def initialise(ref_movie, init_method='cnmf', Ain=None, K=3, ds_factor=1, initbatch=500,
               T1=20000, mot_corr=True,  save_init=False, expected_comps=150,
               rval_thr=0.85, NumROIs=None, thresh_overlap=0.1, decay_time=0.2,
               min_SNR=2.5, merge_thresh=0.85, minibatch_shape=100,
               CNN_filter=False):
               #del_duplicates=False):
    """
    Inputs:
    -------
    
    ref_movie               string
                            reference movie file location
    
    init_method             string
                            choose initialisation method: 'cnmf' (recommended) or 'bare'
                                
    Ain                     2D numpy array
                            binary cell mask for seeded initialisation method
    
    K                       int
                            number of initial components to be found
    
    ds_factor               float
                            spatial downsampling factor (increases speed but may lose some fine structure)
                            
    initbatch               int
                            number of frames for initialization
    
    T1                      int
                            expected length of the whole recording, in frames (overestimate)
                            
    mot_corr                Boolean
                            flag for online motion correction
                            
    save_init               Boolean
                            flag for saving initialization object
                            
    expected_compts         int
                            maximum number of expected components used for memory pre-allocation (exaggerate here)
                            
    rval_thr                float
                            correlation threshold for new component inclusion; preferrable range: 0.8-0.9;
                            default for CNMF: 0.9
                            
    NumROIs                 int
                            maximum number of ROIs to be tracked    
    
    thresh_overlap          int
                            allowed overlap between detected cells; set to 0 for online analysis (0.5 default)
                            
    decay time              float
                            approximate length of transient event in seconds (0.5 default)
                            
    min_SNR                 float
                            minimum SNR for accepting new components
    
    merge_thresh            float
                            merging correlation threshold for initialisation
                            default: 0.8-0.85, strict: 0.65 (may merge seperate cells if close)
                            
    """
    
    
    # set up some additional supporting parameters needed for the algorithm (these are default values but change according to dataset characteristics)
    t1 = time_()
    fr = 30                                                             # frame rate (Hz)
    gSig = (10,10)                                                      # expected half size of neurons
    p = 1                                                               # order of AR indicator dynamics
    gnb = 1                                                             # number of background components
    gSig = tuple(np.ceil(np.array(gSig)/ds_factor).astype('int'))       # recompute gSig if downsampling is involved
    max_shift = np.ceil(10./ds_factor).astype('int')                    # maximum allowed shift during motion correction
    thresh_cnn = 0.1                                                    # threshold for CNN classifier
        
    
    # not necessary but useful
    if NumROIs == None:
        max_comp_update_shape = expected_comps
    else:
        max_comp_update_shape = NumROIs
        
    max_comp_update_shape = expected_comps                              # number of shapes to be updated each time (put this to a finite small value to increase speed - otherwise to np.inf)
    N_samples = np.ceil(fr*decay_time)                                  # number of timesteps to consider when testing new neuron candidates
    thresh_fitness_raw = scipy.special.log_ndtr(-min_SNR)*N_samples     # exceptionality threshold
    
    # Initialize movie
    if ds_factor > 1:                                   # load only the first initbatch frames and possibly downsample them
        Y = cm.load(ref_movie, subindices = slice(0,initbatch,None)).astype(np.float32).resize(1. / ds_factor, 1. / ds_factor)
    else:
        Y = cm.load(ref_movie, subindices = slice(0,initbatch,None)).astype(np.float32)
        
    initbatch = Y.shape[0] # ensure the correct initbatch info is stored given the file used

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
                                     gSig=gSig, p=p, 
                                     minibatch_shape=100,
                                     minibatch_suff_stat=5,
                                     update_num_comps=True, rval_thr=rval_thr,
                                     thresh_fitness_raw=thresh_fitness_raw,
                                     thresh_overlap=thresh_overlap,
                                     batch_update_suff_stat=True, max_comp_update_shape=max_comp_update_shape, 
                                     deconv_flag=False, use_dense=True,
                                     simultaneously=False, n_refit=0)
    
    elif init_method == 'cnmf':
        stride = None
        T = Y.shape[0]  # initbatch
        n_processes = np.maximum(np.int(psutil.cpu_count()),1)
        images = np.reshape(Yr.T, [T] + list(dims), order='F')
    
        cnm_init = cnmf.CNMF(n_processes, k=K, gSig=gSig,
                         merge_thresh=merge_thresh, p=p,
                         #rf=patch_size//2,
                         stride=stride,
                         simultaneously=False, # True: slower but more accurate (doesn't seem slower actually)
                         del_duplicates=True, 
                         use_dense=True,
                         thresh_overlap=thresh_overlap,
                         remove_very_bad_comps=True, # apparently can be less precise if true?
                         skip_refinement=False,
                         normalize_init=False, options_local_NMF=None,
                         minibatch_shape=minibatch_shape, minibatch_suff_stat=5,
                         update_num_comps=True, rval_thr=rval_thr,
                         thresh_fitness_delta=-50, gnb=gnb,
                         thresh_fitness_raw=thresh_fitness_raw,
                         batch_update_suff_stat=True, max_comp_update_shape=max_comp_update_shape)
    

        cnm_init = cnm_init.fit(images)
        
    elif init_method == 'seeded':
        stride = None
        T = Y.shape[0]  # initbatch
        n_processes = np.maximum(np.int(psutil.cpu_count()),1)
        images = np.reshape(Yr.T, [T] + list(dims), order='F')
        
        cnm_init = seeded_initialization(Y.transpose(1, 2, 0), Ain=Ain, init_batch=initbatch, gnb=gnb,
                                        # n_processes - default of 2 inside seeded_init
                                         #k=K, 
                                         gSig=gSig,
                                         merge_thresh=merge_thresh, p=p,
                                         stride=stride,
                                         simultaneously=True, # True: slower but more accurate (doesn't seem slower actually)
                                         del_duplicates=True, 
                                         use_dense=True,
                                         thresh_overlap=thresh_overlap,
                                         remove_very_bad_comps=True, # apparently can be less precise if true?
                                         skip_refinement=False,
                                         normalize_init=False, options_local_NMF=None,
                                         minibatch_shape=100, minibatch_suff_stat=5,
                                         update_num_comps=True, rval_thr=rval_thr,
                                         thresh_fitness_delta=-50, #gnb=gnb,
                                         thresh_fitness_raw=thresh_fitness_raw,
                                         batch_update_suff_stat=True, max_comp_update_shape=max_comp_update_shape)
    
        
    print(('Number of components:' + str(cnm_init.A.shape[-1])))

    # Find centre and radius for cells detected
    coms_init = com(cnm_init.A, d1, d2)
    
    if CNN_filter:
        print('Filtering the components')
        A = cnm_init.A
        
        predictions, final_crops = evaluate_components_CNN(
            A, dims, gSig, model_name=os.path.join(caiman_datadir(), 'model', 'cnn_model'))

        idx_keep = np.where(predictions[:, 1] >= thresh_cnn)[0].tolist()
#        coms_init = coms_init[idx_keep]


    # Setting idx_components to select which cells are kept
    if CNN_filter:
        idx_components=idx_keep
    else:
        idx_components=None

    if ds_factor == 1:
        last_frame = Y[-1]
    else:
        last_frame = cv2.resize(Y[-1], (512, 512), interpolation=cv2.INTER_CUBIC)
        #last_frame = cm.load(ref_movie, subindices = slice(initbatch-1,initbatch,None)).astype(np.float32)
        
        
#    cnm2 = deepcopy(cnm_init)
#    cnm2.opsin = None
    
    
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
    init_values['Yr'] = Yr
    init_values['idx_components'] = idx_components
    #    init_values['cnm2'] = cnm2

    if CNN_filter:
        init_values['thresh_cnn'] = thresh_cnn
    
    if save_init:
        cnm_init.dview = None
        filename = ref_movie[:-4] + '_init_' + init_method + '_DS_' + str(ds_factor) + '.pkl'
        save_object(init_values, filename)
        print('Init file saved as ' + filename)
#        init_values = load_object(ref_movie[:-4] + '_DS_' + str(ds_factor) + '.pkl')
    

    t2 = time_()
    print('Initialisation time: ' + str(t2-t1)[:6] + ' s')
    
    return last_frame, init_values


if __name__ == '__main__':
    ref_movie = r'C:\Users\Patrycja\caiman_data\sample_vid\20170728_A329_t001Substack (195-540).tif'
    lframe, initv = initialise(ref_movie, ds_factor=1.5, initbatch=500)