import numpy as np
import cv2
from copy import deepcopy
import scipy
import matplotlib.pyplot as pl

import time
import sys

if sys.platform == 'win32':
    time = time.clock
else:
    time = time.time

import caiman as cm
from caiman.motion_correction import motion_correct_iteration_fast
from caiman.utils.visualization import view_patches_bar
from caiman.base.rois import com

"""
# create a function for plotting results in real time if needed
def create_frame(cnm2,img_norm,captions):
    A, b = cnm2.Ab[:, cnm2.gnb:], cnm2.Ab[:, :cnm2.gnb].toarray()
    C, f = cnm2.C_on[cnm2.gnb:cnm2.M, :], cnm2.C_on[:cnm2.gnb, :]
    comps_frame = A.dot(C[:,t-1]).reshape(cnm2.dims, order = 'F')*img_norm/np.max(img_norm)   # inferred activity due to components (no background)
    bgkrnd_frame = b.dot(f[:,t-1]).reshape(cnm2.dims, order = 'F')*img_norm/np.max(img_norm)  # denoised frame (components + background)
    if show_residuals:
        all_comps = np.reshape(cnm2.Yres_buf.mean(0),cnm2.dims, order='F')*img_norm/np.max(img_norm)
        all_comps = np.minimum(np.maximum(all_comps*10,0),255)
    else:
        all_comps = (np.array(A.sum(-1)).reshape(cnm2.dims, order = 'F'))                         # spatial shapes
    frame_comp_1 = cv2.resize(np.concatenate([frame_/np.max(img_norm),all_comps*3.],axis = -1),(2*np.int(cnm2.dims[1]*resize_fact),np.int(cnm2.dims[0]*resize_fact) ))
    frame_comp_2 = cv2.resize(np.concatenate([comps_frame*10.,comps_frame+bgkrnd_frame],axis = -1),(2*np.int(cnm2.dims[1]*resize_fact),np.int(cnm2.dims[0]*resize_fact) ))
    frame_pn = np.concatenate([frame_comp_1,frame_comp_2],axis=0).T
    vid_frame = np.repeat(frame_pn[:,:,None],3,axis=-1)
    vid_frame = np.minimum((vid_frame*255.),255).astype('u1')
    cv2.putText(vid_frame,captions[0],(5,20),fontFace = 5, fontScale = 1.2, color = (0,255,0), thickness = 1)
    cv2.putText(vid_frame,captions[1],(np.int(cnm2.dims[0]*resize_fact) + 5,20),fontFace = 5, fontScale = 1.2, color = (0,255,0), thickness = 1)
    cv2.putText(vid_frame,captions[2],(5,np.int(cnm2.dims[1]*resize_fact)  + 20),fontFace = 5, fontScale = 1.2, color = (0,255,0), thickness = 1)
    cv2.putText(vid_frame,captions[3],(np.int(cnm2.dims[0]*resize_fact) + 5 ,np.int(cnm2.dims[1]*resize_fact)  + 20),fontFace = 5, fontScale = 1.2, color = (0,255,0), thickness = 1)
    cv2.putText(vid_frame,'Frame = '+str(t),(vid_frame.shape[1]//2-vid_frame.shape[1]//10,vid_frame.shape[0]-20),fontFace = 5, fontScale = 1.2, color = (0,255,255), thickness = 1)      
    return vid_frame"""

# Run OnACID and optionally plot results in real time
def run_onacid(online_file, init_values, save_results=False, mot_corr=True):
               #Cn_init, cnm_init, img_norm, img_min, dims, Yr, T1 = 10000,  # make them object parameters
               #ds_factor=1, gnb=1, 
               
    img_norm = init_values['img_norm']
    img_min = init_values['img_min'] 
    dims = init_values['dims']
    Yr = init_values['Yr']
    T1 = init_values['T1']
    ds_factor = init_values['ds_factor']
    gnb = init_values['gnb']
    
    Cn_init = init_values['Cn_init']
    cnm_init = init_values['cnm_init']
    
    
    max_shift = np.ceil(10./ds_factor).astype('int')
    cnm2 = deepcopy(cnm_init)
    cnm2.Ab_epoch = []                       # save the shapes at the end of each epoch
    t = 0#cnm2.initbatch                       # current timestep
    tottime = []
    Cn = Cn_init.copy()
    epochs = 1
    
    #plot_contours_flag = True               # flag for plotting contours of detected components at the end of each file
    play_reconstr = False                    # flag for showing video with results online (turn off flags for improving speed)
    save_movie = False                       # flag for saving movie (file could be quite large..)
    #movie_name = folder_name + '/output.avi' # name of movie to be saved
    #resize_fact = 1                        # image resizing factor
    
    """if online_files == 0:                    # check whether there are any additional files
        process_files = fls[:init_files]     # end processing at this file
        init_batc_iter = [initbatch]         # place where to start
        end_batch = T1              
    else:
        process_files = fls[:init_files + online_files]     # additional files
        init_batc_iter = [initbatch] + [0]*online_files     # where to start reading at each file"""
    
    
    shifts = []
    show_residuals = False
    if show_residuals:
        caption = 'Mean Residual Bufer'
    else:
        caption = 'Identified Components'
        
    captions = ['Raw Data','Inferred Activity',caption,'Denoised Data']
    
    """if save_movie and play_reconstr:
        #fourcc = cv2.VideoWriter_fourcc('8', 'B', 'P', 'S') 
        fourcc = cv2.VideoWriter_fourcc(*'XVID')
        out = cv2.VideoWriter(movie_name,fourcc, 30.0, tuple([int(2*x*resize_fact) for x in cnm2.dims]))"""
    
    for iter in range(epochs):    
        """if iter > 0:
            process_files = fls[:init_files + online_files]     # if not on first epoch process all files from scratch
            init_batc_iter = [0]*(online_files+init_files)      """
            
        #for file_count, ffll in enumerate(process_files):  # np.array(fls)[np.array([1,2,3,4,5,-5,-4,-3,-2,-1])]:
            #print('Now processing file ' + ffll)
            #Y_ = cm.load(ffll, subindices=slice(init_batc_iter[file_count],T1,None))
            
        Y_ = cm.load(online_file, subindices=slice(0,T1,None))
        print("film shape: ",Y_.shape)
        
        """if plot_contours_flag:   # update max-correlation (and perform offline motion correction) just for illustration purposes
            if ds_factor > 1:
                Y_1 = Y_.resize(1. / ds_factor, 1. / ds_factor, 1)
            else:
                Y_1 = Y_.copy()                
            if mot_corr:
                templ = (cnm2.Ab.data[:cnm2.Ab.indptr[1]] * cnm2.C_on[0, t - 1]).reshape(cnm2.dims, order='F') * img_norm        
                newcn = (Y_1 - img_min).motion_correct(max_shift, max_shift, template=templ)[0].local_correlations(swap_dim=False)                
                Cn = np.maximum(Cn, newcn)
            else:
                Cn = np.maximum(Cn, Y_1.local_correlations(swap_dim=False))"""
    
        old_comps = cnm2.N                              # number of existing components
        for frame_count, frame in enumerate(Y_):        # now process each file
            if np.isnan(np.sum(frame)):
                raise Exception('Frame ' + str(frame_count) + ' contains nan')
            if t % 100 == 0:
                print('Epoch: ' + str(iter+1) + '. ' + str(t)+' frames have beeen processed in total. '+str(cnm2.N - old_comps)+' new components were added. Total number of components is '+str(cnm2.Ab.shape[-1]-gnb))
                old_comps = cnm2.N
    
            t1 = time()                                 # count time only for the processing part
            frame_ = frame.copy().astype(np.float32)    # 
            if ds_factor > 1:
                frame_ = cv2.resize(frame_, img_norm.shape[::-1])   # downsampling
    
            frame_ -= img_min                                       # make data non-negative
    
            if mot_corr:                                            # motion correct
                templ = cnm2.Ab.dot(cnm2.C_on[:cnm2.M, t - 1]).reshape(cnm2.dims, order='F') * img_norm
                frame_cor, shift = motion_correct_iteration_fast(frame_, templ, max_shift, max_shift)
                shifts.append(shift)
            else:
                templ = None
                frame_cor = frame_
    
            frame_cor = frame_cor / img_norm                        # normalize data-frame
            cnm2.fit_next(t, frame_cor.reshape(-1, order='F'))      # run OnACID on this frame
            tottime.append(time() - t1)                             # store time
    
            t += 1
            
            """if t % 1000 == 0 and plot_contours_flag:
                pl.cla()
                A = cnm2.Ab[:, cnm2.gnb:]
                crd = cm.utils.visualization.plot_contours(A, Cn, thr=0.9)  # update the contour plot every 1000 frames
                pl.pause(1)"""
                
            """if play_reconstr:                                               # generate movie with the results
                vid_frame = create_frame(cnm2,img_norm,captions)             
                if save_movie:
                    out.write(vid_frame)
                cv2.imshow('frame',vid_frame)                
                if cv2.waitKey(1) & 0xFF == ord('q'):
                    break     """                           
    
        print('Cumulative processing speed is ' + str(t/np.sum(tottime))[:5] + ' frames per second.')
    cnm2.Ab_epoch.append(cnm2.Ab.copy())                        # save the shapes at the end of each epoch
            
    """if save_movie:
        out.release()
    cv2.destroyAllWindows()"""
    
    #%  save results (optional)
    #save_results = False
    if save_results:
        np.savez('results_analysis_online_MOT_CORR.npz',
                 Cn=Cn, Ab=cnm2.Ab, Cf=cnm2.C_on, b=cnm2.b, f=cnm2.f,
                 dims=cnm2.dims, tottime=tottime, noisyC=cnm2.noisyC, shifts=shifts)
    
    # extract results from the objects and do some plotting
    A, b = cnm2.Ab[:, cnm2.gnb:], cnm2.Ab[:, :cnm2.gnb].toarray()
    C, f = cnm2.C_on[cnm2.gnb:cnm2.M, t-t//epochs:t], cnm2.C_on[:cnm2.gnb, t-t//epochs:t]
    noisyC = cnm2.noisyC[:,t-t//epochs:t]
    b_trace = [osi.b for osi in cnm2.OASISinstances] if hasattr(cnm2, 'OASISinstances') else [0]*C.shape[0]
    
    
    pl.figure()
    crd = cm.utils.visualization.plot_contours(A, Cn, thr=0.9)
    view_patches_bar(Yr, scipy.sparse.coo_matrix(A.tocsc()[:, :]), C[:, :], b, f,
                     dims[0], dims[1], YrA=noisyC[cnm2.gnb:cnm2.M] - C, img=Cn)
    
    
    #Extract x and y positions of cells 
    coms = com(A, dims[0], dims[1])
    
    return coms, cnm2