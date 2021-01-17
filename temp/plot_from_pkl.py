# convert caiman results saved by pickle to .mat
# will save in the same dir as the same name
import os
import _pickle
import scipy.io
from scipy.sparse import issparse
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
root = tk.Tk()
root.withdraw()
file_full_name = r'D:\DeflectionData\20191029\20191029_OG418_t_0013_rtaoi_DS_2.0_OnlineProc_162308.pkl'
if file_full_name == '':
	file_full_name = filedialog.askopenfilename()
print(file_full_name)
with open(file_full_name, 'rb') as f:
	# file opened
	file_data = _pickle.load(f)

# parse cnmf object
cnm = file_data['cnm2']
A, b = cnm.Ab[:, cnm.gnb:], cnm.Ab[:, :cnm.gnb].toarray()
C, f = cnm.C_on[cnm.gnb:cnm.M], cnm.C_on[:cnm.gnb]
S = cnm.S
t = cnm.t
onlineC = file_data['online_C']
num_frames = C.shape[1]
num_rois = C.shape[0]

# parse params
params = file_data['p']

# plot trace and stim time
plot_range = np.arange(file_data['t_init'],file_data['t_cnm'])


test_trace = np.transpose(onlineC[1:10,plot_range])
filt_trace = np.transpose(scipy.signal.medfilt(test_trace, [1,5]))
plt.figure

plt.plot(test_trace)
plt.plot(filt_trace)

plt.show('hold')

plt.figure
plt.scatter(test_trace,filt_trace)




