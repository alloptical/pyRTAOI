function [cell_struct] = get_cell_photo_auc_raw(cell_struct,cmp_fd,ref_fd,opt,varargin)
% compare photostim vs nonphotostim trials, df
% compare frames with stim vs baseline frames
num_shuf = 500;
avg_frame_range = opt.sta_avg_frames;
baseline_avg_frames = opt.baseline_avg_frames;
% for v = 1:numel(varargin)
%     if strcmpi(varargin{v},'baseline_avg_frames')
%         baseline_avg_frames = varargin{v+1};
%     end
% end
for i = 1:size(cell_struct,2)
    all_ref = mean(cell_struct(i).(ref_fd)(avg_frame_range,:),1) - mean(cell_struct(i).(ref_fd)(baseline_avg_frames,:),1);
    all_ref = max(cell_struct(i).(ref_fd)(avg_frame_range,:),[],1) - max(cell_struct(i).(ref_fd)(baseline_avg_frames,:),[],1);

    save_fd_name = [cmp_fd,'_photo_ranksum'];
    all_stim = mean(cell_struct(i).(cmp_fd)(avg_frame_range,:),1) - mean(cell_struct(i).(cmp_fd)(baseline_avg_frames,:),1);
    all_stim = max(cell_struct(i).(cmp_fd)(avg_frame_range,:),[],1) - max(cell_struct(i).(cmp_fd)(baseline_avg_frames,:),[],1);

    num_sample = min([length(all_ref)],length(all_stim));
    num_sample = min([10,num_sample]);
    indices = randperm(length(all_ref));
    indices = indices(1:num_sample);
    ref_sample =all_ref(indices);
    indices = randperm(length(all_stim));
    indices = indices(1:num_sample);
    stim_sample =all_stim(indices);
    if median(stim_sample)<median(ref_sample)
        cell_struct(i).(save_fd_name) = -ranksum(stim_sample,ref_sample);
    else
        cell_struct(i).(save_fd_name) = ranksum(stim_sample,ref_sample);
    end
    
%     
    save_fd_name = [cmp_fd,'_photoauc'];
    labels = [ones(1,length(all_ref)),2.*ones(1,length(all_stim))]';
    scores = [ all_ref  all_stim]';
    [~,~,~, activeAUC] = perfcurve(labels,scores,2);
    cell_struct(i).(save_fd_name) = activeAUC;
    

   
    
end

