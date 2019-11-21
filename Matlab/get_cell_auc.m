function [cell_struct] = get_cell_auc(cell_struct,cmp_fds,save_fd_name,opt)
num_shuf = 300;
try
peak_frame_range = opt.sta_peak_search_range;
end
avg_frame_range = opt.sta_avg_frames;
for i = 1:size(cell_struct,2)
    if ~ opt.flag_use_peak
        all_stim1 = mean(cell_struct(i).(cmp_fds{1})(avg_frame_range,:),1);
        all_stim2 = mean(cell_struct(i).(cmp_fds{2})(avg_frame_range,:),1);
    else
        all_stim1 = max(cell_struct(i).(cmp_fds{1})(peak_frame_range,:),[],1);
        all_stim2 = max(cell_struct(i).(cmp_fds{2})(peak_frame_range,:),[],1);

    end
labels = [ones(1,length(all_stim1)),2.*ones(1,length(all_stim2))]';
scores = [ all_stim1  all_stim2]';
scores(scores<0.01)=0;
[~,~,~, correct_stimulusAUC] = perfcurve(labels,scores,2);

% shuffle to get zscore stim auc
shuf_stim_auc = nan(1,num_shuf);
parfor s = 1:num_shuf
    shuf_labels = labels(randperm(length(labels))');
    [~,~,~, shuf_stim_auc(s)] = perfcurve(shuf_labels,scores,1);
end
cell_struct(i).([save_fd_name '_zscore']) = (correct_stimulusAUC-mean(shuf_stim_auc))/std(shuf_stim_auc);
cell_struct(i).(save_fd_name) = correct_stimulusAUC;
end

end

