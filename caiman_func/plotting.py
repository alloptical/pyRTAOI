#### PLOTTING TOOLS ######

import matplotlib.pyplot as pl
import numpy as np
from scipy.sparse import issparse

import caiman as cm
from caiman.utils.visualization import view_patches_bar, plot_contours

#%% Contour plot on correlation matrix
pl.cla()
A = cnm2.Ab[:, cnm2.gnb:]
# update the contour plot every 1000 frames
crd = cm.utils.visualization.plot_contours(A, Cn, thr=0.9)

#%% Interactive tool to see all cells
view_patches_bar([], A, # or if A is sparse: scipy.sparse.coo_matrix(A.tocsc()[:, :]),
                 C[:, :t], b, f[:,:t],
                 dims[0], dims[1], YrA=YrA[:,:t], img=Cn)


#%% View patches- a bit weird. secs = 0 for interactive, 1 for automatic scrolling through cells
cm.utils.visualization.view_patches([], A, C[:, :t], b[:,:t], f, dims[0], dims[1], YrA=YrA[:, :], secs = 1)

#%% Plot countours - doesn't plot anything?

plot_contours(A, Cn, thr=None, thr_method='max', maxthr=0.2, nrgthr=0.9, display_numbers=True, max_number=None,
                  cmap=None, swap_dim=False, colors='w', vmin=None, vmax=None)

#%% To get and plot coordinates of each cell
coords = cm.utils.visualization.get_contours(A, dims)

cell = 0
# for each cell in coords
shape_ = coords[cell]['coordinates']
pl.figure();pl.plot(*shape_.T)

#%% Recreating 2D spatial info from matrix A
# from view_patches_bar

# if needed
# A = A.todense()
    
d1, d2 = dims
img = np.reshape(np.array(A.mean(axis=1)), (d1, d2), order='F') # 341, 341
pl.figure();pl.imshow(img)
pl.colorbar()  # range of values: 0 (no cell) - 0.007 --> what does it represent?

# Add contours onto the plot
contours_on = False
if contours_on:
    for cell in coords:
        shape_ = cell['coordinates']
        pl.plot(*shape_.T)
        
# Add COMS (centre of mass of cells) onto the plot
coms_on = 0
if coms_on:
    coms = com(A, *dims)
    pl.plot(coms[:,1], coms[:,0], '.r')
