### PLOTTING FUNCTIONS ###



#def plot_contours(A, Cn, thr=None, thr_method='max', maxthr=0.2, nrgthr=0.9, display_numbers=True, max_number=None,
#                  cmap=None, swap_dim=False, colors='w', vmin=None, vmax=None, **kwargs):
#    """Plots contour of spatial components against a background image and returns their coordinates
#
#     Parameters:
#     -----------
#     A:   np.ndarray or sparse matrix
#               Matrix of Spatial components (d x K)
#
#     Cn:  np.ndarray (2D)
#               Background image (e.g. mean, correlation)
#
#     thr_method: [optional] string
#              Method of thresholding:
#                  'max' sets to zero pixels that have value less than a fraction of the max value
#                  'nrg' keeps the pixels that contribute up to a specified fraction of the energy
#
#     maxthr: [optional] scalar
#                Threshold of max value
#
#     nrgthr: [optional] scalar
#                Threshold of energy
#
#     thr: scalar between 0 and 1
#               Energy threshold for computing contours (default 0.9)
#               Kept for backwards compatibility. If not None then thr_method = 'nrg', and nrgthr = thr
#
#     display_number:     Boolean
#               Display number of ROIs if checked (default True)
#
#     max_number:    int
#               Display the number for only the first max_number components (default None, display all numbers)
#
#     cmap:     string
#               User specifies the colormap (default None, default colormap)
#
#     Returns:
#     --------
#     Coor: list of coordinates with center of mass, contour plot coordinates and bounding box for each component
#    """

from caiman.base.rois import com
import numpy as np
from scipy.sparse import issparse
from warnings import warn
from matplotlib import pyplot as pl
from past.utils import old_div
import time
from numpy import linalg as LA

from caiman.utils.visualization import get_contours, plot_shapes

#%%
thr=None
thr_method='max'
maxthr=0.2
nrgthr=0.9
display_numbers=True
max_number=None
cmap=None
swap_dim=False
colors='c'  #change color
vmin=None
vmax=None

# inputs = A & dims
# optional input: draw real shapes vs draw circles of similiar radius --> needed for GUI and SLM
# + the other variables from above

A = cnm_init.A.tocsc()

if issparse(A):
    A = np.array(A.todense())
else:
    A = np.array(A)
#
#if swap_dim:
#    Cn = Cn.T
#    print('Swapping dim')
#%%
d1, d2 = dims #np.shape(Cn)
d, nr = np.shape(A)
if max_number is None:
    max_number = nr

if thr is not None:
    thr_method = 'nrg'
    nrgthr = thr
    warn("The way to call utilities.plot_contours has changed. Look at the definition for more details.")

x, y = np.mgrid[0:d1:1, 0:d2:1]

ax = pl.gca()
#if vmax is None and vmin is None:
#    pl.imshow(Cn, interpolation=None, cmap=cmap,
#              vmin=np.percentile(Cn[~np.isnan(Cn)], 1), vmax=np.percentile(Cn[~np.isnan(Cn)], 99))
#else:
#    pl.imshow(Cn, interpolation=None, cmap=cmap,
#              vmin=vmin, vmax=vmax)

coordinates = []
cm = com(A, d1, d2)
#%%
for i in range(np.minimum(nr, max_number)):
    #pars = dict(kwargs)
    if thr_method == 'nrg':
        indx = np.argsort(A[:, i], axis=None)[::-1]
        cumEn = np.cumsum(A[:, i].flatten()[indx]**2)
        cumEn /= cumEn[-1]
        Bvec = np.zeros(d)
        Bvec[indx] = cumEn
        thr = nrgthr

    else:  # thr_method = 'max'
        if thr_method != 'max':
            warn("Unknown threshold method. Choosing max")
        Bvec = A[:, i].flatten()
        Bvec /= np.max(Bvec)
        thr = maxthr

    if swap_dim:
        Bmat = np.reshape(Bvec, dims, order='C') #np.shape(Cn), ---- is the whole Cn reversed as well?
    else:
        Bmat = np.reshape(Bvec, dims, order='F') #np.shape(Cn),
        
    cs = pl.contour(y, x, Bmat, [thr], colors=colors)  #plotting contours! why reversed? for coms? why not reverse coms?
    
   ########### calculate radius for each centre of mass
    
    
    # this fix is necessary for having disjoint figures and borders plotted correctly
    p = cs.collections[0].get_paths()
    v = np.atleast_2d([np.nan, np.nan])
    for pths in p:
        vtx = pths.vertices
        num_close_coords = np.sum(np.isclose(vtx[0, :], vtx[-1, :]))
        if num_close_coords < 2:
            if num_close_coords == 0:
                # case angle
                newpt = np.round(old_div(vtx[-1, :], [d2, d1])) * [d2, d1]
                #import ipdb; ipdb.set_trace()
                vtx = np.concatenate((vtx, newpt[np.newaxis, :]), axis=0)

            else:
                # case one is border
                vtx = np.concatenate((vtx, vtx[0, np.newaxis]), axis=0)
                #import ipdb; ipdb.set_trace()

        v = np.concatenate(
            (v, vtx, np.atleast_2d([np.nan, np.nan])), axis=0)
# store variables - useful or not? use neuron_id as roiidx?
    pars = {}
    pars['CoM'] = np.squeeze(cm[i, :])
    pars['coordinates'] = v
    pars['bbox'] = [np.floor(np.min(v[:, 1])), np.ceil(np.max(v[:, 1])),
                    np.floor(np.min(v[:, 0])), np.ceil(np.max(v[:, 0]))]
    pars['neuron_id'] = i + 1
    coordinates.append(pars)
#%%
if display_numbers:
    for i in range(np.minimum(nr, max_number)):
        if swap_dim:
            ax.text(cm[i, 0], cm[i, 1], str(i + 1), color=colors)
        else:
            ax.text(cm[i, 1], cm[i, 0], str(i + 1), color=colors)

#return coordinates


#%% get contours: returns contours of comps, coms and neuron id --> maybe just use this 
# time both things! this and up
# center of mass here from scipy/ndimage/measurements.py. but same results as the caiman.base.rois function
t0 = time.clock()
coordinates = get_contours(A, dims, thr=0.9)
print(time.clock()-t0)

#%% plot shapes: plots cells detected
from scipy.ndimage.filters import median_filter
Ab = cnm2.Ab # cnm_init does not have Ab as such. A and b separate there
plot_shapes(Ab, dims, num_comps=15, size=(15, 15), comps_per_row=None,
                cmap='viridis', smoother=lambda s: median_filter(s, 3))

#%%
view_patches(Yr, A, C, b, f, d1, d2, YrA=None, secs=1)

#%% remove nan values from contour coordinates
ex = coordinates[0]['coordinates']
ex2 = np.isnan(ex)
good_values = ~ex2.all(1)
ex_filtered = ex[good_values]
#%% calculate radius: e.g. max distance from com to contour
# need this whenever new component is detected to determine radius
comp_radius = np.zeros((cnm_init.A.shape[-1],1))

t0 = time.clock()

# may be quicker to just do com and look at all points in comps rather than just contours? but it is done only once
# or: shape update in the interface too?

for index, cell in enumerate(coordinates):
    # remove nan values from contour coordinates
    cell['coordinates'] = cell['coordinates'][~np.isnan(cell['coordinates']).all(1)] 
    
    cell_coms = np.array([cell['CoM'][1], cell['CoM'][0].copy()])  # reverse x and y for coms
    # need coms from elsewhere because these updated only once every 100 frames
    diff = cell['coordinates'] - cell_coms # may be quicker to do this for all A pixel values instead doing get_contours -- but contours
    # would be nice for plotting maybe too
    dist = LA.norm(diff, axis = 1)
    comp_radius[index] = max(dist) # find the maximum radius
    
print(time.clock()-t0)
#%% OnACID radius finding

no_cells = cnm2.Ab[:, cnm2.gnb:].shape[1]
#%%
# 11th cell only:
one_cell_coords = get_contours(cnm2.Ab[:,cnm2.gnb+10], dims, thr = 0.9)  #wrong neurons id but thats okay?

ex = one_cell_coords[0]['coordinates'] # filter nans
ex2 = np.isnan(ex)
good_values = ~ex2.all(1)
contour1 = ex[good_values]

#%% all cells for comparison
all_coords = get_contours(cnm2.Ab[:,cnm2.gnb:], dims, thr = 0.9)

ex = all_coords[10]['coordinates']
ex2 = np.isnan(ex)
good_values = ~ex2.all(1)
contour2 = ex[good_values]

#%% check they are the same:
same = (contour1 == contour2).all()  # true



#%% do com only on new comp -- possible?
x = com(cnm2.Ab[:,cnm2.gnb+10], dims[0], dims[1])
