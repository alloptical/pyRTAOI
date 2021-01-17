#%% ########### UNUSED // NOT WORKING IDEAS FOR C1V1 MASKING #############


#%% all components from A1 and A2 detected as empty -- not useful
from caiman.base.rois import extract_binary_masks_blob
import scipy 

masks, pos_examples, neg_examples = extract_binary_masks_blob(scipy.sparse.csr_matrix(A1), neuron_radius=7, dims=dims) # mask has to be sparse

#%% visualise masks - downsampling actually helped with masks (fewer 'cells' found)
mask = A1
#d1, d2 = Y.shape[-2:]
img = np.reshape(np.array(mask.mean(axis=1)), (d1, d2), order='F') # 341, 341
pl.figure();pl.imshow(img)

#%% doesnt work here cause A1 and A2 are already in 2D
#from caiman.base.rois import mask_to_2d
#A1_2D = mask_to_2d(A1)


#%% 3d mask for opsin
comps = A1.shape[-1]

for cell in range(comps):
    single_mask = np.reshape(np.array(A1[:,cell]),(d1,d2),order='F')[np.newaxis,:,:]
    
    if cell == 0: 
        mask1_3d = single_mask[np.newaxis,:,:]
    elif cell > 0:
        mask1_3d = np.concatenate([mask1_3d, single_mask[np.newaxis,:,:]])
    
mask1_3d = np.squeeze(mask1_3d)

#%% 3d mask for calcium
comps = A2.shape[-1]

for cell in range(comps):
    single_mask = np.reshape(np.array(A2[:,cell]),(d1,d2),order='F')[np.newaxis,:,:]
    
    if cell == 0: 
        mask2_3d = single_mask[np.newaxis,:,:]
    elif cell > 0:
        mask2_3d = np.concatenate([mask2_3d, single_mask[np.newaxis,:,:]])
    
mask2_3d = np.squeeze(mask2_3d)

#%% doesnt seem useful for extracting mask either... 
# masks to input: components x d1 x d2
idx_tp_gt, idx_tp_comp, idx_fn_gt, idx_fp_comp, performance = nf_match_neurons_in_binary_masks(mask1_3d, mask2_3d, labels = ['opsin mask', 'calcium mask'])

# output: 
#FOV: 0, shape: 282,232 total cost: 231.236668
#0.35302019119262695
#{'recall': 0.0, 'precision': 0.0, 'accuracy': 0.0, 'f1_score': 0.0}

#%% Correlation between images
#from scipy.signal import correlate2d

# suuuper slow and wrong: corr_img shape = (1023, 1023)
#corr_img = correlate2d(calcium_img, opsin_img)

#%%
from caiman.base.rois import extract_binary_masks_from_structural_channel

img = opsin_img

Y = np.reshape(img,(1,512,512)) # for cm function
A, mr = extract_binary_masks_from_structural_channel(Y)
# A.shape = (262144, 581) for c1v1 image

d1, d2 = Y.shape[-2:]
mask1 = np.reshape(np.array(A.mean(axis=1)), (d1, d2), order='F') # 341, 341
pl.figure();pl.imshow(mask1)



#%%
from scipy.signal import fftconvolve
# Correlation is convolution with one input reversed
corr = fftconvolve(opsin_img, calcium_img[::-1, ::-1])  # basically same ouput as corr_img but fast
pl.figure();pl.imshow(corr)

#%% masks should be of dimensions (dim1*dim2) x num of cells
# masks A: ndarray or csc_matrix  # pixels x # of components
dims = d1, d2
from caiman.base.rois import register_ROIs


matched_ROIs1, matched_ROIs2, non_matched1, non_matched2, performance, A2 =  register_ROIs(A1, A2, dims, plot_results = True)

# error:
#No loop matching the specified signature and casting
#was found for ufunc true_divide


#%% load data  -- error appears but loaded
#C:\Users\Patrycja\Anaconda3\envs\caiman\lib\site-packages\tifffile\tifffile.py:6794: RuntimeWarning: invalid value encountered in true_divide
#  'ZDistance': values[:, 0] / values[:, 1],
#fname = r'C:\Users\Patrycja\Desktop\Project M\Zoe\orig-forPatrycja\sample\spontaneous\20171116_OG229\20171116_OG229_s-017_Cycle00001_Ch1_000001.ome.tif'

fname = r'T:\ForPatrycja\pyRTAOI\samples\example1\20171229_OG245_s-026_Cycle00001_Ch1_000001.ome.tif'
Y = cm.load(fname).astype(np.float32)
#Y = Y[np.newaxis,:,:]

#%% visualise
pl.figure();pl.imshow(Y[0,:,:], cmap='gray') # this works

#%% try smoothing // edge enhancing? 

import cv2
import matplotlib.pyplot as plt

img = np.squeeze(opsin_img).astype(np.uint8)
pl.figure();pl.imshow(img)  # kinda interesting

edges = cv2.Canny(img,0,1000)

# this no good
plt.subplot(121),plt.imshow(img,cmap = 'gray')
plt.title('Original Image'), plt.xticks([]), plt.yticks([])
plt.subplot(122),plt.imshow(edges,cmap = 'gray')
plt.title('Edge Image'), plt.xticks([]), plt.yticks([])

plt.show()

#%% doesnt work
#blur = cv2.bilateralFilter(opsin_img,9,75,75)
 
#%% try cv2.threshold
import cv2 as cv
ret, thr = cv.threshold(opsin_img,350,1,cv.THRESH_BINARY)

pl.figure(); pl.imshow(np.squeeze(thr))
