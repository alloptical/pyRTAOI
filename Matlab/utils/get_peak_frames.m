function [cell_struct] = get_peak_frames(cell_struct,cmp_fds,save_fd_name,opt)
num_shuf = 300;

peak_frame_range = opt.sta_peak_search_range;

for i = 1:size(cell_struct,2)

    [all_stim1,peak_frames1] = max(cell_struct(i).(cmp_fds{1})(peak_frame_range,:),[],1);
    [all_stim2.peak_frames2 ]= max(cell_struct(i).(cmp_fds{2})(peak_frame_range,:),[],1);

%% test if most of the peak frames are after stimulus onset
%% 

cell_struct(i).([save_fd_name '_zscore']) = (correct_stimulusAUC-mean(shuf_stim_auc))/std(shuf_stim_auc);

end

end

