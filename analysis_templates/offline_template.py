### OFFLINE ONACID ANALYSIS PIPELINE ### 

"""
    For offline purposes, OnACID can be run multiple times (epochs) or 
    initialised with a mask (seeded initialisatin).
    
    Running the algorithm two times (epochs) can find quite many cells which
    were missed during the first run. Third epoch finds a few extra cells but
    not as many.
    
    The OnACID output cell mask (A) can be used to initialise the algorithm. 
    This allows for tracking the activity of cells in the mask throughout the 
    whole recording.
    
    Strict offline analysis can be done by running OnACID for two epochs
    followed by seeded init + another OnACID run.
    
    C1V1 mask can also be used to initialise the algorithm with cells 
    expressing the opsin. This should be followed by filtering of the 
    initialisation output with the CNN filter as the mask does not separate
    cells well.
    
"""


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
import glob
import numpy as np
from copy import deepcopy
import time
import sys
sys.path.append(r'C:\Users\intrinsic\Desktop\pyRTAOI-rig') # for caiman_func access

if sys.platform == 'win32':
    time_ = time.clock
else:
    time_ = time.time

import cv2
import matplotlib.pyplot as pl
import scipy
from scipy.sparse import issparse

import caiman as cm
from caiman.utils.utils import load_object, save_object
from caiman.base.rois import com
from caiman.utils.visualization import view_patches_bar, plot_contours
from caiman_func.initialisation import initialise
from caiman.motion_correction import motion_correct_iteration_fast
from caiman.paths import caiman_datadir
from caiman.components_evaluation import evaluate_components_CNN
from caiman.base.rois import extract_binary_masks_from_structural_channel

#%% Select reference movie, init parameters and initialise the algorithm or load an init .pkl file

#ref_movie_path = r'C:\Users\intrinsic\Desktop\samples\example5\init_results\20171228_OG241_t-024_Cycle00001_Ch2_init_cnmf_DS_1.5.pkl'
#ref_movie_path = r'C:\Users\intrinsic\Desktop\samples\example1\init_results\20171229_OG245_t-052_Cycle00001_Ch2_substack1-200_init_seeded_DS_1.5.pkl'

# example movies
example = 8

sample_folder = r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\samples'
example_folder = os.path.join(sample_folder, 'example' + str(example))
files = glob.glob(os.path.join(example_folder,'*Ch2.tif'))

if len(files)==1:
    ref_movie_path = files[0]
print('Ref movie path: ', ref_movie_path)

#ref_movie_path = r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\samples\example8\init_results\20180107_OG242_t-001_Cycle00001_Ch2_init_cnmf_DS_1.5.pkl'

movie_ext = ref_movie_path[-4:]

if movie_ext  == '.tif':
    K = 100
    ds_factor = 1.5
    initbatch = 500
    minibatch_shape = 100
    gSig = (4,4)
    expected_comps = 700
    
    lenient = 0
    
    if lenient:
        rval_thr = 0.7
        min_SNR = 1
        thresh_overlap = 0.2 # 0.5
    else:
        rval_thr = 0.8
        min_SNR = 2.5
        thresh_overlap = 0.2
    
    merge_thresh = 0.95     # 0.85
    
    lframe, init_values = initialise(ref_movie_path, init_method='cnmf', K=K,
                                     ds_factor=ds_factor, gSig=gSig,
                                     initbatch=initbatch, rval_thr=rval_thr, 
                                     expected_comps=expected_comps,
                                     thresh_overlap=thresh_overlap, 
                                     save_init=False, mot_corr=True,
                                     merge_thresh=merge_thresh,
                                     minibatch_shape=minibatch_shape,
                                     min_SNR=min_SNR, T1=10000)
    
elif movie_ext == '.pkl':
    init_values = load_object(ref_movie_path)
    

c = init_values

# Extract initialisation parameters
mot_corr = c['mot_corr']
img_norm = c['img_norm']
img_min = c['img_min']
T1 =  c['T1']
ds_factor =  c['ds_factor']
Cn_init =  c['Cn_init']
cnm_init =  c['cnm_init']
gnb =  c['cnm_init'].gnb
gSig = c['cnm_init'].gSig
dims =  c['cnm_init'].dims
initbatch = cnm_init.initbatch
NumROIs = c['expected_comps'] # as set in the interface OR c['expected_comps']
N_samples = c['N_samples']
K = c['K']
coms_init = c['coms_init']
Yr = c['Yr']
coms_init = c['coms_init']


#%% Visualise the results of initialisation
visualise_init = True

if visualise_init:
    pl.figure()
    crd = plot_contours(cnm_init.A.tocsc(), Cn_init, thr=0.9)
    pl.plot(coms_init[:,1], coms_init[:,0], '.r')
    A, C, b, f, YrA, sn = cnm_init.A, cnm_init.C, cnm_init.b, cnm_init.f, cnm_init.YrA, cnm_init.sn
    
    view_patches_bar([], scipy.sparse.coo_matrix(
    A.tocsc()[:, :]), C[:, :], b, f, dims[0], dims[1], YrA=YrA[:, :], img=None)#Cn_init)

#%% Load c1v1 image and create a binary cell mask
mask_file = r'C:\Users\intrinsic\Desktop\pyRTAOI2018\samples\example1\20171229_OG245_s-026_Cycle00001_Ch1_000001.ome.tif'
ds_factor = 1.5

opsin_img = cm.load(mask_file, subindices = slice(0,1,None))
opsin_img = opsin_img[np.newaxis,:,:].resize(1. / ds_factor, 1. / ds_factor)
dims = opsin_img.shape[-2:]

pl.figure(); pl.imshow(np.squeeze(opsin_img), cmap='gray')


# specify cell extraction parameters
min_area_size = 100
gSig_ = 9 #gSig[0]

A_opsin, mr_opsin = extract_binary_masks_from_structural_channel(opsin_img, gSig=gSig_, min_area_size=min_area_size, min_hole_size=15) 
A_opsin = A_opsin.astype('int')
opsin_mask = np.reshape(np.array(A_opsin.max(axis=1)), dims, order='F').astype('int')  # changed mean to max because it's binary

pl.figure(); pl.imshow(opsin_mask)

#%% Compare detected cells with the c1v1 mask
A = cnm_init.A

# Convert A to numpy array
if issparse(A):
    A = np.array(A.todense())
else:
    A = np.array(A)


onacid_mask = (deepcopy(A)>0).astype('int')  # binarise the onacid output mask

opsin = []
overlap_ratio = []

for cell in range(onacid_mask.shape[-1]):
    cell_mask = (np.reshape(onacid_mask[:,cell], dims, order='F'))
    mask_inx = np.where(cell_mask==1)
    cell_pix = sum(sum(cell_mask == 1))
    
    inter = cv2.bitwise_and(opsin_mask, cell_mask)
    inter_pix = sum(sum(inter))
    overlap = inter_pix/cell_pix
    overlap_ratio.append(overlap)
    
    thresh = 0.5
    if overlap <= thresh:
        onacid_mask[:,cell][onacid_mask[:,cell] == 1] = -3
    else:
        onacid_mask[:,cell][onacid_mask[:,cell] == 1] = 3
        
    opsin.append(overlap > thresh)
    print('cell ' + str(cell) + ': ' + str(overlap > thresh))

no_opsin = A.shape[-1] - sum(opsin)

# visualise all comps
summed_A = np.hstack((A_opsin, onacid_mask))
summed_mask = np.reshape(np.array(summed_A.sum(axis=1)), dims, order='F')
pl.figure();pl.imshow(summed_mask)
pl.colorbar()

#%% Seeded initialisation
#ref_movie_path = r'T:\ForPatrycja\pyRTAOI\samples\example2\20171229_OG245_t-053\20171229_OG245_t-053_Cycle00001_Ch2.tif'
ref_movie_path = r'C:\Users\intrinsic\Desktop\pyRTAOI2018\samples\example1\20171229_OG245_t-052_Cycle00001_Ch2_substack1-2700.tif'

opsin_seeded = True
if opsin_seeded:
    Ain = A_opsin # or A from OnACID output (as np array) for seeded init
else:
    Ain = A

ds_factor = 1.5
initbatch = 500
minibatch_shape = 100
    
lframe, seeded_init_values = initialise(ref_movie_path, init_method='seeded', Ain=Ain,
                                 ds_factor=ds_factor, initbatch=initbatch, rval_thr=0.85,
                                 thresh_overlap=0.2, save_init=False, mot_corr=True,
                                 merge_thresh=0.85) # , min_SNR=15)

c = seeded_init_values

# Extract initialisation parameters
mot_corr = c['mot_corr']
img_norm = c['img_norm']
img_min = c['img_min']
T1 =  c['T1']
ds_factor =  c['ds_factor']
Cn_init =  c['Cn_init']
cnm_init =  c['cnm_init']
#cnm2 = c['cnm2']
gnb =  c['cnm_init'].gnb
gSig = c['cnm_init'].gSig # changes to (7,7) on its own?
dims =  c['cnm_init'].dims
initbatch = cnm_init.initbatch
NumROIs = c['expected_comps'] # as set in the interface OR c['expected_comps']
N_samples = c['N_samples']
K = c['K']
coms_init = c['coms_init']
Yr = c['Yr']
coms_init = c['coms_init']

#%% Filter seeded
A, C, b, f, YrA, sn = cnm_init.A, cnm_init.C, cnm_init.b, cnm_init.f, cnm_init.YrA, cnm_init.sn

# threshold for CNN classifier
thresh_cnn = 0.1 # 0.1 default


predictions, final_crops = evaluate_components_CNN(
    A, dims, gSig, model_name=os.path.join(caiman_datadir(), 'model', 'cnn_model'), isGPU=1)

A_exclude, C_exclude = A[:, predictions[:, 1] <
                         thresh_cnn], C[predictions[:, 1] < thresh_cnn]  # elements removed from analysis

ind_rem = np.where(predictions[:,1] < thresh_cnn)[0].tolist()
ind_keep = np.where(predictions[:, 1] >= thresh_cnn)[0].tolist()

A, C = A[:, predictions[:, 1] >=
         thresh_cnn], C[predictions[:, 1] >= thresh_cnn]

noisyC = cnm_init.Cin # [cnm_init.gnb:cnm_init.A.shape[-1]+1] #cnm2.noisyC[cnm2.gnb:cnm2.M]
YrA = noisyC[predictions[:, 1] >= thresh_cnn] - C

#%% Prepare object for OnACID

# Setting idx_components to ind_keep to select only cells that passed filtering (None means all kept)
try:
    idx_components = ind_keep
    coms_init = coms_init[idx_components]
except:
    idx_components = None
    
try: opsin_seeded
except: check_opsin = False

cnm2 = deepcopy(cnm_init)
path_to_cnn_residual = os.path.join(caiman_datadir(), 'model', 'cnn_model_online.h5')

cnm2._prepare_object(np.asarray(c['Yr']), c['T1'], c['expected_comps'], idx_components=idx_components,
                         min_num_trial=2, N_samples_exceptionality=int(c['N_samples']),
                         path_to_model=path_to_cnn_residual)#,
#                         sniper_mode=False, use_peak_max=False, q=0.5) # default


#%% Store info on opsin as object property
if opsin_seeded:
    cnm2.opsin = [True]*len(ind_keep)
else:
    cnm2.opsin = opsin


#%% Visualise filtering of seeded
C, f = cnm2.C_on[cnm2.gnb:cnm2.M], cnm2.C_on[:cnm2.gnb]
A, b = cnm2.Ab[:, cnm2.gnb:cnm2.M], cnm2.Ab[:, :cnm2.gnb]

noisyC = cnm2.noisyC[cnm2.gnb:cnm2.M]
YrA = noisyC - C

visualise_init = True

if visualise_init:
    pl.figure()
    crd = plot_contours(A, Cn_init, thr=0.9)
    coms = c['coms_init']
    pl.plot(coms[:,1], coms[:,0], '.r')
    #b, f, YrA, sn = cnm_init.b, cnm_init.f, cnm_init.YrA, cnm_init.sn
    
    view_patches_bar([], scipy.sparse.coo_matrix(
    A.tocsc()[:, :]), C[:, :500], b, f, dims[0], dims[1], #YrA=np.zeros([A.shape[-1],cnm_init.initbatch]),
    YrA=YrA[:, :500], 
    img=Cn_init)


#%% COMPONENT EVALUATION - doesn't work now. ValueError: cannot reshape array of size 15652 into shape (100,5,28)
# the components are evaluated in three ways:
#   a) the shape of each component must be correlated with the data
#   b) a minimum peak SNR is required over the length of a transient
#   c) each shape passes a CNN based classifier
    
#from caiman.components_evaluation import estimate_components_quality_auto
#
#T = cnm_init.initbatch
#images = np.reshape(Yr.T, [T] + list(dims), order='F')
#
#idx_components, idx_components_bad, SNR_comp, r_values, cnn_preds = \
#    estimate_components_quality_auto(images, cnm_init.A, cnm_init.C, cnm_init.b, cnm_init.f,
#                                     cnm_init.YrA, frate=30, decay_time=0.5, gSig=gSig, dims=dims,
#                                     min_SNR=2.5,
#                                     r_values_min=0.85, use_cnn=False,
#                                     thresh_cnn_min=0.8)

#%% Remove any undesirable cells from initialisation
inx = [0,1] # from visualisation
cnm2.remove_components(inx) 


#%% Define OnACID parameters
NumROIs = cnm2.expected_comps # temp

max_shift = np.ceil(10./ds_factor).astype('int')  # max shift allowed
t = cnm2.initbatch
tottime = []
Cn = Cn_init.copy()
shifts = []
coms = coms_init.copy()
com_count = cnm2.N

BufferLength = 200
BufferPointer = 0
RoiBuffer = np.zeros([NumROIs,BufferLength]).astype('float32')
ROIlist_threshold = np.zeros(NumROIs)

rejected = 0
accepted = list(range(1, cnm2.N+1)) # 0th comp is background
expect_components = True
dist_thresh = 21
display_shift = False   # option to display motion correction shift inside the interface
check_overlap = False    # removed in pyRTAOI
new_spotted = []        # for timing purposes
removed = []


# cnn classifier for removal of unfit cells
remove_flag = False    # flag for removing components with bad shapes during onacid
T_rm = 650    # remove bad components every T_rm frames
rm_thr = 0.1  # CNN classifier removal threshold
cnm2.thresh_CNN_noisy = 0.5

epochs = 1
cnm2.Ab_epoch = []                       # save the shapes at the end of each epoch

#%% Run OnACID multiple times (epochs)
#movie_path = r'C:\Users\intrinsic\caiman_data\sample_vid\stimulated_test\20170329_OG151_t-008_Substack (1-3000)--(500-3000).tif'
movie_path = ref_movie_path

print('Loading video')
try:
    print(Y_.shape)
except:
    Y_ = cm.load(movie_path, subindices=slice(initbatch,T1,None)) # 0,T1,None
    
print('Video loaded')

for iter in range(epochs):  
    for frame_count, frame in enumerate(Y_):        # now process each file
        t1 = time_()
        
        if np.isnan(np.sum(frame)):
            raise Exception('Frame ' + str(frame_count) + ' contains nan')
        
        frame_ = frame.copy().astype(np.float32)
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
            
        if expect_components:
#            if cnm2.C_on[com_count+1,t] != 0:
#                print(cnm2.C_on[com_count+1,t])
#                print(cnm2.N)
             if cnm2.N - (com_count + rejected) == 1:  # cnm2.N - (com_count) == 1: if removing too close cells
                new_coms = com(cnm2.Ab[:, -1], dims[0], dims[1])[0]
                new_spotted.append(t-initbatch)
                
                # Check for repeated components
                if check_overlap:
                    close = abs(coms - new_coms) < dist_thresh/ds_factor
                    repeated = any(np.all(close,axis=1))
                else:
                    repeated = False
                
                if repeated == False:
                    coms = np.vstack((coms, new_coms))
                    y, x = new_coms   # reversed
                    
                    com_count += 1
                    accepted.append(cnm2.N)
                    print('New cell detected (' + str(cnm2.N-rejected) + ')')

#                    tt = time_()
                    if check_opsin:
                        cell_A = np.array(cnm2.Ab[:,-1].todense())
                        cell_mask = (np.reshape(cell_A, dims, order='F') > 0).astype('int')
                        mask_inx = np.where(cell_mask==1)
                        cell_pix = sum(sum(cell_mask == 1))
                        
                        inter = cv2.bitwise_and(opsin_mask, cell_mask)
                        inter_pix = sum(sum(inter))
                        overlap = inter_pix/cell_pix
                                   
                        cnm2.opsin.append(overlap > thresh)
                        print(overlap>thresh)
#                    print(time_()-tt)
                    
                    if com_count == NumROIs:
                        expect_components = False
                        print('Not accepting more components')
                else:
                    print('Repeated component found')
                   # cnm2.remove_components([cnm2.N-1])
                    removed.append(t-initbatch)
                    rejected += 1
    
     
        if mot_corr:
            shift_ = [round(shift[0]), round(shift[1])]
        
        try:
            RoiBuffer[:com_count, BufferPointer] = cnm2.C_on[accepted,t] # -1
#            time.sleep(1)
            
            ROIlist_threshold[:com_count] = np.nanmean(RoiBuffer[:com_count,:], axis=1) + 3*np.nanstd(RoiBuffer[:com_count,:], axis=1)
        except Exception as e:
            print(e)
            print(RoiBuffer[:com_count,:])
        
        t += 1
        
        if BufferPointer==BufferLength-1:
            BufferPointer = 0
        else:
            BufferPointer +=1
            
        
        if remove_flag and t % T_rm == 0:
            prd, _ = evaluate_components_CNN(cnm2.Ab[:, gnb:], dims, gSig)
            ind_rem = np.where(prd[:, 1] < rm_thr)[0].tolist()
            cnm2.remove_components(ind_rem)
            print('Removing '+str(len(ind_rem))+' components')
        
        tottime.append(time_() - t1)                             # store time

    
    print('Cumulative processing speed is ' + str(t/np.sum(tottime))[:5] + ' frames per second.')
    
    # save the shapes at the end of each epoch
    cnm2.Ab_epoch.append(cnm2.Ab.copy())

 
#%% Extract results from object
# Can pass through the CNN classifier with a low threshold (keeps clearer neuron shapes and excludes processes)

use_CNN = 0
if use_CNN:
    A, b = cnm2.Ab[:, cnm2.gnb:cnm2.M], cnm2.Ab[:, :cnm2.gnb]
    C, f = cnm2.C_on[cnm2.gnb:cnm2.M], cnm2.C_on[:cnm2.gnb]

    # threshold for CNN classifier
    thresh_cnn = 0.1 # 0.1
    from caiman.components_evaluation import evaluate_components_CNN
    predictions, final_crops = evaluate_components_CNN(
        A, dims, gSig, model_name=os.path.join(caiman_datadir(), 'model', 'cnn_model'))
    A_exclude, C_exclude = A[:, predictions[:, 1] <
                             thresh_cnn], C[predictions[:, 1] < thresh_cnn]  # elements removed from analysis

    A, C = A[:, predictions[:, 1] >=
             thresh_cnn], C[predictions[:, 1] >= thresh_cnn]
    noisyC = cnm2.noisyC[cnm2.gnb:cnm2.M]
    YrA = noisyC[predictions[:, 1] >= thresh_cnn] - C
    
    ind_rem_CNN = np.where(predictions[:,1] < thresh_cnn)[0].tolist()
    ind_keep_CNN = np.where(predictions[:, 1] >= thresh_cnn)[0].tolist()
    
    # trim long files
    C = C[:,:t]
    f = f[:,:t]
    noisyC = noisyC[:,:t]
    YrA = YrA[:,:t]
else:
    A, b = cnm2.Ab[:, cnm2.gnb:], cnm2.Ab[:, :cnm2.gnb].toarray()
    C, f = cnm2.C_on[cnm2.gnb:cnm2.M, t - t //
                 epochs:t], cnm2.C_on[:cnm2.gnb, t - t // epochs:t] # this does not show init frames for epochs > 1 and may not show all frames of last epochs

# could use t - (t - initbatch) // epochs:t for just the last epoch
# can also stack init (hstack)
    
    noisyC = cnm2.noisyC[:, t - t // epochs:t]
    YrA = noisyC[cnm2.gnb:cnm2.M] - C
    
print(('Number of components:' + str(A.shape[-1])))

# Convert A to numpy array
if issparse(A):
    A = np.array(A.todense())
else:
    A = np.array(A)


b_trace = [osi.b for osi in cnm2.OASISinstances] if hasattr(
    cnm2, 'OASISinstances') else [0] * C.shape[0]
deconvolved = [osi.s for osi in cnm2.OASISinstances]

coms_post = com(A, *dims)

#%% Visualise OnACID output
pl.figure()
crd = cm.utils.visualization.plot_contours(A, Cn, thr=0.9)
show_coms = False
if show_coms:
    pl.plot(coms_post[:,1], coms_post[:,0], '.r') # reverse if swap_dim = False above

view_patches_bar([], A, # scipy.sparse.coo_matrix(A.tocsc()[:, :]),
                 C[:, :t], b, f[:,:t],
                 dims[0], dims[1], YrA=YrA[:,:t], img=None) #img=Cn for background frame, None for cell mask


#%% Plot deconvolved signal for checking
pl.figure();pl.plot(deconvolved[cell][initbatch:])

#%% Save results individually, as npz
save_results = False
folder = os.path.dirname(ref_movie_path)

# can't save Ab or cnm2 (can't access the inside of it) directly

if save_results:        
    saved_data = folder + '\_offline_DS_ ' + str(ds_factor) + '.npz'
    
    np.savez(saved_data,
             Cn=Cn,
             Cf=cnm2.C_on,
             A=A, # can be fully saved if converted to dense numpy array; for sparse would have to save 
             b=b, C=C, f=f, noisyC=noisyC, b_trace=b_trace,
             OASISinstances = cnm2.OASISinstances, # can be saved; no need for b_trace or deconvolved but may be easier
             deconvolved = deconvolved,
             coms=coms, opsin=cnm2.opsin,
             dims=cnm2.dims, ds_factor=ds_factor,
             min_SNR=min_SNR, thresh_overlap=cnm2.thresh_overlap,
             
             tottime=tottime,
             shifts=shifts)

#%% Reload and unpack saved data
file = r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\samples\example1\offline_results\ex1_cnmf_results.npz'
npz_datafile = np.load(file)
locals().update(npz_datafile)

#%%
#A = npz_datafile['cnm_A']

#%% Save results (whole cnm object) as pkl

save_dict = dict()
save_dict['accepted'] = accepted
save_dict['t_cnm'] = t
save_dict['Cn'] = Cn
save_dict['coms'] = coms

# parameters
save_dict['ds_factor'] = ds_factor
save_dict['K'] = K
save_dict['min_SNR'] = min_SNR
save_dict['gSig'] = gSig
save_dict['rval_thr'] = rval_thr
save_dict['thresh_overlap'] = thresh_overlap
save_dict['merge_thresh'] = merge_thresh
save_dict['expected_comps'] = expected_comps
    
    
folder = os.path.join(os.path.dirname(ref_movie_path), 'offline_results')
if not os.path.exists(folder):
    os.makedirs(folder)
    
if use_CNN:
    saveResultPath = os.path.join(folder, 'ex' + str(example) + '_onacid_results_ds_' + str(ds_factor) + '_CNN_' + str(thresh_cnn) + '.pkl')
    
    cnm3 = deepcopy(cnm2)
    cnm3.remove_components(ind_rem_CNN)  # remove cells rejected by CNN
    
    save_dict['cnm2'] = cnm3
    save_dict['coms'] = coms[ind_keep_CNN]

else:
    if lenient:
        saveResultPath = os.path.join(folder, 'lenient_ex' + str(example) + '_onacid_results_ds_' + str(ds_factor) + '.pkl')
    else: 
        saveResultPath = os.path.join(folder, 'ex' + str(example) + '_onacid_results_ds_' + str(ds_factor) + '.pkl')

    save_dict['cnm2'] = cnm2
    save_dict['coms'] = coms

save_object(save_dict, saveResultPath)

save_mat = True # additionally save the mat file
if save_mat:
    from pkl2mat import pkl2mat

    pkl2mat(file_full_name = saveResultPath)

#%% Load pkl-ed cnm object
file1 = 'substack1-200_DS_1.5_OnlineProc.pkl'
file2 = 'substack1-1000_DS_1.5_OnlineProc.pkl'
saveResultPath = r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\samples\example1\pyrtaoi_results\20171229_OG245_t-052_Cycle00001_Ch2_' + file1
cnm_object = load_object(saveResultPath)

cnm2 = cnm_object['cnm2']
t = cnm_object['t_cnm']
epochs = 1

C, f = cnm2.C_on[cnm2.gnb:cnm2.M], cnm2.C_on[:cnm2.gnb]

#%% Load pkl init object
saved_pkl = r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\samples\example1\init_results\20171229_OG245_t-052_Cycle00001_Ch2_init_cnmf_DS_1.5.pkl'
pkl_datafile = load_object(saved_pkl)

#mask = pkl_datafile['cnm_init'].A # access mask

#%%
init_values = pkl_datafile
idx_components = init_values['idx_components']
path_to_cnn_residual = os.path.join(caiman_datadir(), 'model', 'cnn_model_online.h5')

cnm2 = deepcopy(init_values['cnm_init'])
cnm2._prepare_object(np.asarray(init_values['Yr']), init_values['T1'], 
                     init_values['expected_comps'], idx_components=idx_components,
                     min_num_trial=2, N_samples_exceptionality=int(init_values['N_samples']),
                     path_to_model=path_to_cnn_residual)

#%% Load npz object
saved_npz = r'T:\ForPatrycja\pyRTAOI\samples\example1\_offline_DS_ 1.5.npz'
npz_datafile = np.load(saved_npz)

mask = npz_datafile['A']

#%% Check and plot opsin overlap offline
onacid_mask = (deepcopy(A)>0).astype('int')  # binarise the onacid output mask

opsin = []
overlap_ratio = []

for cell in range(onacid_mask.shape[-1]):
    cell_mask = (np.reshape(onacid_mask[:,cell], dims, order='F'))
    mask_inx = np.where(cell_mask==1)
    cell_pix = sum(sum(cell_mask == 1))
    
    inter = cv2.bitwise_and(opsin_mask, cell_mask)
    inter_pix = sum(sum(inter))
    overlap = inter_pix/cell_pix
    overlap_ratio.append(overlap)
    
    thresh = 0.5
    if overlap <= thresh:
        onacid_mask[:,cell][onacid_mask[:,cell] == 1] = -3
    else:
        onacid_mask[:,cell][onacid_mask[:,cell] == 1] = 3
        
    opsin.append(overlap > thresh)
    print('cell ' + str(cell) + ' : ' + str(overlap > thresh))

no_opsin = cnm2.N - sum(opsin)

# visualise all comps
summed_A = np.hstack((A_opsin, onacid_mask))
summed_mask = np.reshape(np.array(summed_A.sum(axis=1)), dims, order='F')
pl.figure();pl.imshow(summed_mask)
pl.colorbar()


# Check offline == online
print(opsin == cnm2.opsin)

#%% Testing RoiBuffer
folder = r'C:\Users\intrinsic\Desktop\pyRTAOI20180530\analysis templates\signal plotting online testing\plots_f0-2500_mbs100_suffstat1_simultaneouslyFalse'

save_plots = 0

spotted = deepcopy(new_spotted)
spotted.reverse()

for cell in range(cnm2.N):
#cell = 3
    pl.figure()
    
    pl.subplot(311); pl.plot(cnm2.noisyC[cell+1,500:t]); pl.ylabel('noisyC') # t-1
    pl.title('cell ' + str(cell));
    pl.subplot(312); pl.plot(cnm2.C_on[cell+1,500:t]); pl.ylabel('C_on')
    pl.subplot(313); pl.plot(RoiBuffer[cell,:]); pl.ylabel('RoiBuffer')
    if cell >= K:
        fr = spotted.pop()
        pl.plot([fr, fr], pl.ylim(), 'r')
        
    if save_plots:
        pl.savefig(folder + '\cell_' + str(cell))
        pl.close()

        
#%% Plot time of each frame to see the delay introduced by various functionalities
tottime = np.array(tottime)

duration = tottime.shape[0] 
shape_ref = minibatch_shape # this is minibatch_shape i think
cnn_check = []

for i in range(1,t):
    if i % T_rm == 0:
        cnn_check.append(i - initbatch)

pl.figure(); pl.plot(tottime*1000)
pl.plot([0,duration],[33,33],'--m',label = 'Image acquisition time')
pl.plot(np.arange(shape_ref,duration+1,shape_ref)-1, np.ones([int((t-cnm2.initbatch)/shape_ref),1])*9.5, 'g.', label='Cell shapes refreshed')
pl.plot(new_spotted,np.ones([len(new_spotted),1])*10,'r.', label='New cell detected')

removed_cells = False
if removed_cells:
    pl.plot(removed,np.ones([len(removed),1])*12, 'm*', label='Cell removed')

if remove_flag:
    pl.plot(cnn_check, np.ones([len(cnn_check),1])*8.5, 'b^', label='cnn check')
    
    
pl.title('Time per frame for the online pipeline')
pl.xlabel('Frame')
pl.ylabel('Processing time (ms)')
#pl.ylim([8,100])
pl.legend()

#%% Visualise noise
cell = 5
pl.figure();

pl.subplot(131)
pl.plot(YrA[cell,:],'r',label='noise')
pl.plot(C[cell,:],'g',label='C')
y_max = int(pl.ylim()[1] + 2)
y_min = int(pl.ylim()[0] - 0.5)
pl.ylim(y_min, y_max)
pl.legend()

pl.subplot(132)
pl.plot(noisyC[cell+1,:],'b',label='noisyC')
pl.legend()
pl.ylim(y_min, y_max)

pl.subplot(133)
pl.plot(RoiBuffer[cell,:t],'m',label='buffer')
pl.legend()
pl.ylim(y_min, y_max)
