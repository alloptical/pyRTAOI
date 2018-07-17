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
root = tk.Tk()
root.withdraw()

def convert2mat(file_full_name = '',save_full_name = ''):
    if file_full_name == '':
        file_full_name = filedialog.askopenfilename()
    print(file_full_name)
    with open(file_full_name, 'rb') as f:
        # file opened
        file_data = _pickle.load(f)
    p = Path(file_full_name)
    file_dir = str(p.parents[0])
    file_name = p.name
    save_name = file_name.replace('.pkl','.mat')
    if save_full_name == '':
        save_full_name = os.path.join(file_dir,save_name)
        
    save_dict = dict()    
    
    save_dict['accepted'] = file_data['accepted']
    save_dict['t_cnm'] = file_data['t_cnm']
    save_dict['coms'] = file_data['coms']
    
    # parse cnmf object
    cnm = file_data['cnm2']
    A, b = cnm.Ab[:, cnm.gnb:], cnm.Ab[:, :cnm.gnb].toarray()
    C, f = cnm.C_on[cnm.gnb:cnm.M], cnm.C_on[:cnm.gnb]
    
    if issparse(A):
        A = np.array(A.todense())
    save_dict['cnm_A'] = A
    save_dict['cnm_b'] = b
    save_dict['cnm_C'] = C #cnm.C
    save_dict['cnm_f'] = f #cnm.f
    save_dict['cnm_N'] = cnm.N
    save_dict['cnm_gnb'] = cnm.gnb
    save_dict['num_frames_init'] = cnm.initbatch


    # cell traces
    save_dict['noisyC'] = cnm.noisyC
#    save_dict['deconvC'] = cnm.C_on
    save_dict['cnm_dims'] = cnm.dims
 
        
        
        
    print('saved as '+save_full_name)
    scipy.io.savemat(save_full_name, mdict = save_dict)

    
if __name__ == '__main__':
    convert2mat()