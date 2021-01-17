import numpy as np
import pylab as pl
import sys
sys.path.append(r'C:\Users\intrinsic\Desktop\pyRTAOI20180530') # for caiman_func access

from caiman_func.initialisation import initialise

import caiman as cm
from caiman.source_extraction import cnmf as cnmf
from caiman.utils.visualization import view_patches_bar, plot_contours
from copy import deepcopy
from scipy.special import log_ndtr

#%% Load two example images
ds_factor = 1.5

opsin_img = cm.load(r'C:\Users\intrinsic\Desktop\pyRTAOI20180530\samples\example2\20171229_OG245_s-026_Cycle00001_Ch1_000001.ome.tif',
                    subindices = slice(0,1,None))
calcium_img = cm.load(r'C:\Users\intrinsic\Desktop\pyRTAOI20180530\samples\example2\20171229_OG245_s-025_Cycle00001_Ch2_000001.ome.tif',
                      subindices = slice(0,1,None))

opsin_img = opsin_img[np.newaxis,:,:].resize(1. / ds_factor, 1. / ds_factor)
calcium_img = calcium_img[np.newaxis,:,:,].resize(1. / ds_factor, 1. / ds_factor)

#%% Display the images
import matplotlib.pyplot as pl

pl.figure(); pl.imshow(np.squeeze(opsin_img), cmap='gray')
pl.figure(); pl.imshow(np.squeeze(calcium_img), cmap='gray')

#%% Extract binary masks from the images

# gSig = blocksize and it has to be gSig % 2 == 1 and > 1
# default inputs: (Y, min_area_size=30, min_hole_size=15, gSig=5, expand_method='closing' / dilation , selem=np.ones((3, 3)))
# gsig default of 5, more cells spotted at 7
# min hole size doesn't make much difference but area size does - it removes the smallest comps found

from caiman.base.rois import extract_binary_masks_from_structural_channel, nf_match_neurons_in_binary_masks

d1, d2 = opsin_img.shape[-2:]

A1, mr1 = extract_binary_masks_from_structural_channel(opsin_img, gSig=7, min_area_size=100, min_hole_size=15) 
mask1 = np.reshape(np.array(A1.max(axis=1)), (d1, d2), order='F').astype('int')  # change mean to max cause it's a binary mask
# A.shape = (262144, 581) for c1v1 image: 512x512 and 581 elements/cells

A2, mr2 = extract_binary_masks_from_structural_channel(calcium_img, gSig=9, min_area_size=100) # not that good here4
mask2 = np.reshape(np.array(A2.max(axis=1)), (d1, d2), order='F').astype('int')

#%%
pl.figure();pl.imshow(mask1)
pl.figure();pl.imshow(mask2)

#%% Find intersection between the masks
from cv2 import bitwise_and

inter = bitwise_and(mask1, mask2)
pl.figure();pl.imshow(inter)

#%% Compare detected cells with the c1v1 mask

# Convert A from onacid to numpy array
if issparse(A):
    A = np.array(A.todense())
else:
    A = np.array(A)


from cv2 import bitwise_and

onacid_mask = (deepcopy(A)>0).astype('int')  # binarise the onacid output mask
opsin_mask = mask1

cell_ok = []
overlap_ratio = []

dims = d1,d2

for cell in range(onacid_mask.shape[-1]):
    cell_mask = (np.reshape(onacid_mask[:,cell], dims, order='F'))
    mask_inx = np.where(cell_mask==1)
    cell_pix = sum(sum(cell_mask == 1))
    
    inter = bitwise_and(opsin_mask, cell_mask)
    inter_pix = sum(sum(inter))
    overlap = inter_pix/cell_pix
    overlap_ratio.append(overlap)
    
    thresh = 0.7
    if overlap <= thresh:
        onacid_mask[:,cell][onacid_mask[:,cell] == 1] = -3
        
    cell_ok.append(overlap > thresh)
    
#    cell_mask[cell_mask > 1e-10] = 1  # binarise the masks
#    check = cell_mask in mask1
    print('cell ' + str(cell) + ' : ' + str(overlap>thresh))

no_opsin = cnm2.N - sum(cell_ok)

# visualise all comps
summed_A = np.hstack((A1, onacid_mask))
summed_mask = np.reshape(np.array(summed_A.sum(axis=1)), dims, order='F')
pl.figure();pl.imshow(summed_mask)
pl.colorbar()

