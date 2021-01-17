import os
import scipy.io
import tkinter as tk
import _pickle
from tkinter import filedialog
from pathlib import Path
root = tk.Tk()
root.withdraw()

def pybehav_pkl2mat(file_full_name = '',save_full_name = ''):
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
	print('saving...')
	scipy.io.savemat(save_full_name, mdict = file_data)
	print('saved as '+save_full_name)

if __name__ == '__main__':
	pybehav_pkl2mat()