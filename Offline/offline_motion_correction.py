"""
# motion correct movies with given shifts 
# untested quick draft
"""

import cv2
from skimage.external import tifffile
import numpy as np


def motion_correct_given_shifts(img, shifts):
    """ For using in online realtime scenarios """
    h_i, w_i = img.shape
    M = np.float32([[1, 0, shifts[0]], [0, 1, shifts[1]]])
    new_img = cv2.warpAffine(
        img, M, (w_i, h_i), flags=cv2.INTER_CUBIC, borderMode=cv2.BORDER_REFLECT)
    return new_img

if __name__ == '__main__':
	
	rt_filename = ''
	movie_filename = ''

	# load pyRTAOI output
	with open(rt_filename, 'rb') as input_obj:
		rt_obj = pickle.load(input_obj)
	rt_shifts = rt_obj['shifts']
	ds_factor = rt_obj['ds_factor']
	shifts = rt_shifts*ds_factor

	# load movie
	movie = tifffile.TiffFile(self.p['moviePath'], multifile=True).asarray()

	# get movie dimensions
	num_frames = movie.shape[0]
	frame_dims = movie.shape[1:]

	for i in range(num_frames):
		this_frame = movie[i]
		reg_frame = motion_correct_given_shifts(this_frame,shifts[i])





