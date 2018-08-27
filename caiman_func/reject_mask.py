"""
Creating a binary reject mask

"""

import numpy as np
import matplotlib.pyplot as pl
import cv2
import caiman as cm

#%% as done before in GUI

def getPixCorr(img,color):
    xcoors = list()
    ycoors = list()

    for i in range(0,img.width()-1):
        for j in range(0,img.height()-1):
            if img.pixel(i,j) == color:
                xcoors.append(i)
                ycoors.append(j)
                
    return xcoors, ycoors
            
def create_pixmap(center,radius):
#				def color():
#					r = random.randrange(0, 255)
#					g = random.randrange(0, 255)
#					b = random.randrange(0, 255)
#					return QColor(r, g, b)

	pixmap = QPixmap(512, 512)
	painter = QPainter()
	painter.begin(pixmap)
	thiscolor = color()
	painter.setBrush(thiscolor)
	painter.drawEllipse(center,radius,radius)
	painter.end()
	thisRgb = thiscolor.rgb()

	return pixmap, thisRgb #, thiscolor
            
com = [100,100]
ROIsize = 20
dims = [512,512]

mask = np.zeros(dims)

currentRadius = ROIsize/2 #self.thisROI.size()/2
currentCenter = com+currentRadius #self.thisROI.pos()+currentRadius
			
print(currentCenter)
print(currentRadius)

pixelmap,color,qcolor = create_pixmap(currentCenter,currentRadius[0])
thisImage = pixelmap.toImage()
xcoors,ycoors = getPixCorr(thisImage,color)

#%% Manual selection
coms = np.array([[100,100],[200,200],[300,300],[400,400],[500,500]])
x, y = com

ROIsize = 20
dims = [512,512]
ds_factor = 2
dims_res = tuple([int(dim/ds_factor) for dim in dims])

#%%
radius = 10
dims_ = dims

if dims_ == dims_res:
    radius = radius/ds_factor

rejected_mask = np.zeros(dims_)

for x,y in coms:
#    ROI_mask_ = cv2.resize(mask,dims_res)
    ROI_mask = np.zeros(dims_)
    xcoords = []
    ycoords = []
    
    if dims_ == dims_res:
        x = round(x/ds_factor)
        y = round(y/ds_factor)

    for pixel_x in range(dims_[0]):
        for pixel_y in range(dims_[1]):
            
            dx = pixel_x - x # round(x/ds_factor)
            dy = pixel_y - y # round(y/ds_factor)
            dist_squared = dx**2 + dy**2
            
            if dist_squared <= radius**2:
                xcoords.append(pixel_x)
                ycoords.append(pixel_y)
                
    ROI_mask[xcoords,ycoords] = 1
    rejected_mask += ROI_mask

if dims_ == dims:
    rejected_mask = cv2.resize(rejected_mask,dims_res)   # shapes more round this way
    

pl.figure();
pl.imshow(rejected_mask)

# use bitwise and to check if a cell overlaps with reject mask --> if yes (>overlap), count it as rejected!

#%% Automatic selection
file = r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\tests on rig\20180822\20180822_OG299_s-001\20180822_OG299_s-001_Cycle00001_Ch2_000001.ome.tif'
#file = r'\\live.rd.ucl.ac.uk\ritd-ag-project-rd00g6-mhaus91\forPat\tests on rig\20180822\20180822_OG299_s-003\20180822_OG299_s-003_Cycle00001_Ch2_000001.ome.tif'

gcamp_img = cm.load(file)
gcamp_img = gcamp_img[np.newaxis,:,:].resize(1. / ds_factor, 1. / ds_factor)

pl.figure();pl.imshow(np.squeeze(gcamp_img), cmap='gray')

#%% create binary mask --> this not v good for filled cells
from caiman.base.rois import extract_binary_masks_from_structural_channel

d1, d2 = dims_res

A_rej, mr_rej = extract_binary_masks_from_structural_channel(gcamp_img, gSig=3, min_area_size=30, min_hole_size=15)
mask_rej = np.reshape(np.array(A_rej.max(axis=1)), (d1, d2), order='F').astype('int')  # change mean to max cause it's a binary mask

pl.figure();pl.imshow(mask_rej)

#%% 8-bit image better
gcamp_img8 = np.squeeze((gcamp_img/256).astype('uint8'))
pl.figure();pl.imshow(gcamp_img8)

#%%
mask_rej = np.zeros_like(gcamp_img8)
mask_rej[np.where(gcamp_img8>10)] = 1  # 20 works better for second example

mask_rej[np.where(gcamp_img>3e3)] = 1

pl.figure(); pl.imshow(np.squeeze(mask_rej))

#%% filtering of small objects - works!
from skimage.morphology import remove_small_objects, remove_small_holes, dilation

mask_rej_bin = gcamp_img8>10
mask_rej_bin = remove_small_objects(mask_rej_bin, min_size=50)  # 10 for zoom = 1.14; 50 for zoom = 2

pl.figure(); pl.imshow(mask_rej_bin)

#%% check if any coms are within the rejected mask
mask_rej_bin_ = mask_rej_bin.astype('int')
[x,y] = np.where(mask_rej_bin_>0)

pl.figure(); pl.imshow(mask_rej_bin_)
pl.plot(y[0],x[0],'r.')    # ---> why is this inverted?

#%% detect blobs with cv2 --> doesn't detect anything atm?
    
#im = cv2.imread(file, cv2.IMREAD_GRAYSCALE)

# Set up the detector with default parameters
#detector = cv2.SimpleBlobDetector()
    
gcamp_img8 = np.squeeze((gcamp_img/256).astype('uint8'))

# Set up the SimpleBlobdetector with specified params
params = cv2.SimpleBlobDetector_Params()

# Filter by Area
#params.filterByArea = True
#params.minArea = 1
#
## Disable unwanted filter criteria params
#params.filterByInertia = False
#params.filterByConvexity = False

   
# Change thresholds
params.minThreshold = 0;
params.maxThreshold = 256;
 
# Filter by Area.
params.filterByArea = True
params.minArea = 20
 
# Filter by Circularity
params.filterByCircularity = 0 # True
params.minCircularity = 0.1
 
# Filter by Convexity
params.filterByConvexity = 0 # True
params.minConvexity = 0.5
 
# Filter by Inertia
params.filterByInertia = 0 #True
params.minInertiaRatio = 0.5



detector = cv2.SimpleBlobDetector_create(params)

# Detect blobs
keypoints = detector.detect(gcamp_img8)


 
## Draw detected blobs as red circles.
# cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS ensures the size of the circle corresponds to the size of blob
im_with_keypoints = cv2.drawKeypoints(gcamp_img8, keypoints, np.array([]), (0,0,255), cv2.DRAW_MATCHES_FLAGS_DRAW_RICH_KEYPOINTS)
 
# Show keypoints
cv2.imshow("Keypoints", im_with_keypoints)
cv2.waitKey(0)
