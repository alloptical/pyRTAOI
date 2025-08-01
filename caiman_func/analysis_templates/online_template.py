### ONLINE ONACID ANALYSIS PIPELINE ### 
"""
    Basic OnACID pipeline for online purposes.
    
    Additional functionalities such as seeded initialisation or c1v1 masking
    shown in offline_template.py.
    
"""


import os
import numpy as np
from copy import deepcopy
import time
import sys
sys.path.append(r'C:\Users\Patrycja\Desktop\pyRTAOI20180530') # for caiman_func access

if sys.platform == 'win32':
    time_ = time.clock
else:
    time_ = time.time

import cv2
import matplotlib.pyplot as pl
import scipy
from scipy.sparse import issparse

import caiman as cm
from caiman.base.rois import com
from caiman.utils.visualization import view_patches_bar, plot_contours
from caiman_func.initialisation import initialise
from caiman.motion_correction import motion_correct_iteration_fast
from caiman.paths import caiman_datadir
from caiman.components_evaluation import evaluate_components_CNN

#%% select reference movie, init parameters and initialise the algorithm
ref_movie_path = r'C:\Users\intrinsic\caiman_data\sample_vid\stimulated_test\20170329_OG151_t-008_Substack (1-3000)--(1-500).tif'

K = 20
ds_factor = 1.5
initbatch = 500

lframe, init_values = initialise(ref_movie_path, init_method='cnmf', K=K, 
                                 ds_factor=ds_factor, initbatch=initbatch, 
                                 rval_thr=0.85, thresh_overlap=0, save_init=True, 
                                 mot_corr=True, merge_thresh=0.85)
#                                 NumROIs = 50)

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
gSig = c['cnm_init'].gSig # changed to (7,7)?
dims =  c['cnm_init'].dims
NumROIs = 40 # as set in the interface (and input into initialise) OR c['expected_comps']
N_samples = c['N_samples']
K = c['K']
coms_init = c['coms_init']

#%% Visualise the results of initialisation
visualise_init = True

if visualise_init:
    pl.figure()
    crd = plot_contours(cnm_init.A.tocsc(), Cn_init, thr=0.9)
    A, C, b, f, YrA, sn = cnm_init.A, cnm_init.C, cnm_init.b, cnm_init.f, cnm_init.YrA, cnm_init.sn
    view_patches_bar([], scipy.sparse.coo_matrix(
    A.tocsc()[:, :]), C[:, :], b, f, dims[0], dims[1], YrA=YrA[:, :], img=Cn_init)

#%% Define OnACID parameters
pl.close('all')

max_shift = np.ceil(10./ds_factor).astype('int')  # max shift allowed
cnm2 = deepcopy(cnm_init)
t = cnm2.initbatch
tottime = []
Cn = Cn_init.copy()
shifts = []
coms = coms_init.copy()
com_count = cnm2.N

BufferLength = 200
BufferPointer = 0
RoiMeanBuffer = np.zeros([NumROIs,BufferLength])
ROIlist_threshold = np.zeros(NumROIs)

rejected = 0
accepted = list(range(1,cnm2.N+1))
expect_components = True
dist_thresh = 21
display_shift = False   # option to display motion correction shift inside the interface
check_overlap = False    # option removed in pyRTAOI (always on)
new_spotted = []        # for timing purposes
removed = []


# cnn classifier for removal of unfit cells
remove_flag = False    # flag for removing components with bad shapes
T_rm = 650    # remove bad components every T_rm frames
rm_thr = 0.1  # CNN classifier removal threshold
path_to_cnn_residual = os.path.join(caiman_datadir(), 'model', 'cnn_model_online.h5')
cnm2.thresh_CNN_noisy = 0.5

#%% Run OnACID once (online)
movie_path = r'C:\Users\Patrycja\caiman_data\sample_vid\stimulated_test\20170329_OG151_t-008_Substack (1-3000)--(500-3000).tif'
#movie_path = r'C:\Users\Patrycja\CaImAn\sample_vid\spontaneous_test\20170329_OG151_t-006_2_Substack_(1000-3000).tif'

print('Loading video')
Y_ = cm.load(movie_path, subindices=slice(0,1400,None)) # T1
print('Video loaded')
            
for frame_count, frame in enumerate(Y_):        # now process each file
    t1 = time.clock()
    
    if BufferPointer==BufferLength-1:
        BufferPointer = 0
    else:
        BufferPointer +=1
    
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
        RoiMeanBuffer[:com_count, BufferPointer] = cnm2.C_on[accepted,t]
        ROIlist_threshold[:com_count] = np.nanmean(RoiMeanBuffer[:com_count,:], axis=1) + 3*np.nanstd(RoiMeanBuffer[:com_count,:], axis=1)
    except Exception as e:
        print(e)
        print(RoiMeanBuffer[:com_count,:])
    
    t += 1
    
    if remove_flag and t % T_rm == 0:
        prd, _ = evaluate_components_CNN(cnm2.Ab[:, gnb:], dims, gSig)
        ind_rem = np.where(prd[:, 1] < rm_thr)[0].tolist()
        cnm2.remove_components(ind_rem)
        print('Removing '+str(len(ind_rem))+' components')
    
    tottime.append(time.clock() - t1)                             # store time

    
print('Cumulative processing speed is ' + str(t/np.sum(tottime))[:5] + ' frames per second.')


#%% Extract results from the objects and visualise them
epochs = 1
A, b = cnm2.Ab[:, cnm2.gnb:], cnm2.Ab[:, :cnm2.gnb].toarray()

if issparse(A):
    A_ = np.array(A.todense())
else:
    A_ = np.array(A)

C, f = cnm2.C_on[cnm2.gnb:cnm2.M, :], cnm2.C_on[:cnm2.gnb, :]
noisyC = cnm2.noisyC
b_trace = [osi.b for osi in cnm2.OASISinstances] if hasattr(
    cnm2, 'OASISinstances') else [0] * C.shape[0]

deconvolved = [osi.s for osi in cnm2.OASISinstances]

#%% Save results
save_results = False

if save_results:        
    saved_data = movie_path[:-4] + '_offline_DS_ ' + str(ds_factor) + '.npz'
    
    np.savez(movie_path[:-4] + '_offline_DS_ ' + str(ds_factor) + '.npz',
             Cn=Cn,
#             Ab=cnm2.Ab, # can't access it either now
             Cf=cnm2.C_on, 
             A=A_, # seems to work now
             b=b, C=C, f=f, noisyC=noisyC, b_trace=b_trace,
             OASISinstances = cnm2.OASISinstances, # can be saved like this too! no need for b_trace or deconvolved, but may be easier
             deconvolved = deconvolved,
             coms=coms,
             dims=cnm2.dims, ds_factor=ds_factor,
             tottime=tottime,
             shifts=shifts)


#%% Reload and unpack saved data
npz_datafile = np.load(saved_data)
locals().update(npz_datafile)

#%% Visualise OnACID output
pl.figure()
crd = cm.utils.visualization.plot_contours(A, Cn, thr=0.9)
pl.plot(coms[:,1], coms[:,0], '.r') # reverse if swap_dim = False above

view_patches_bar([], 
#                 A,
                 scipy.sparse.coo_matrix(A.tocsc()[:, :]),
                 C[:, :], b, f,
                 dims[0], dims[1], YrA=noisyC[cnm2.gnb:cnm2.M] - C, img=Cn)


#%% Plot time of each frame to see the delay introduced by various functionalities
tottime = np.array(tottime)

duration = tottime.shape[0] 
shape_ref = 100 # this is minibatch_shape i think
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