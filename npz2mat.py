# convert caiman results saved as .npz to .mat
# will save in the same dir as the same name 
import os
import numpy as np
import scipy.io
from scipy.sparse import issparse
import tkinter as tk
from tkinter import filedialog
from pathlib import Path
root = tk.Tk()
root.withdraw()

def npz2mat(file_full_name = '',save_full_name = ''):
    if file_full_name == '':
        file_full_name = filedialog.askopenfilename()
    print(file_full_name)

    file_data = np.load(file_full_name)
    p = Path(file_full_name)
    file_dir = str(p.parents[0])
    file_name = p.name
    save_name = file_name.replace('.npz','.mat')
    if save_full_name == '':
        save_full_name = os.path.join(file_dir,save_name)
        
    save_dict = dict()    
    
    
    # parse cnmf object
    A = file_data['A']
    if issparse(A):
        A = np.array(A.todense())
        
    save_dict['cnm_A'] = A
    save_dict['cnm_b'] = file_data['b']
    save_dict['cnm_f'] = file_data['f']
    
    # cell traces
    save_dict['cnm_C'] = file_data['C']
    save_dict['noisyC'] = file_data['Y_r']
    
    # extracted DF/F values
    save_dict['F_dff'] = file_data['F_dff']
    save_dict['C_df'] = file_data['C_df'] # cleaner signal

    # other
    save_dict['cnm_N'] = file_data['cnm_N']
    save_dict['cnm_gnb'] = file_data['gnb']
    save_dict['cnm_dims'] = file_data['dims']
    save_dict['coms'] = file_data['coms']

    # optional: Cn, thesh_overlap, YrA, F_dff, C_df, F_dff_no_noise
        
    print('saved as '+save_full_name)
    scipy.io.savemat(save_full_name, mdict = save_dict)

    
if __name__ == '__main__':
    npz2mat()