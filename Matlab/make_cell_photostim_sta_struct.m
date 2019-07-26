function [cell_struct] = make_cell_photostim_sta_struct(cnm_struct,cell_struct,accepted_idx,photo_stim_frames,photo_sequence_idx,opt)
%% get sta for each cell
% for photoexcitability test 
% ZZ 2019
% photo_stim_frames: all frame indices with photostim
% photo_sequence_idx: target roi index corresponding to each photostim frame
num_shuf = 100;
num_comp = size(cnm_struct,2);
num_cells = numel(accepted_idx);
frames_to_avg = opt.sta_pre_frames+ opt.frames_with_photo+[1:opt.sta_avg_frames];
for i = 1:num_cells
    this_cnm_idx = accepted_idx(i);
    this_cell_trace = cnm_struct(this_cnm_idx).deconvC_full;
    this_stim_frames = photo_stim_frames(photo_sequence_idx == i);
    this_num_trials = numel(this_stim_frames);
    cell_struct(i).cnm_idx = this_cnm_idx;
    cell_struct(i).num_trials = this_num_trials;
    cell_struct(i).sta_amp = 0;
    cell_struct(i).sta_traces = [];
    cell_struct(i).sta_trace = [];

    if(this_num_trials>0)
        % average across trials
        [~,~,~,~,~,cell_struct(i).sta_traces,cell_struct(i).sta_trace] = make_sta_from_traces(this_cell_trace,this_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        cell_struct(i).sta_amp = mean(cell_struct(i).sta_trace(frames_to_avg));
        
        
        % using ROC comparing sta amp to baseline
        % to do: alternatively do what harvey did - compare when the cell is
        % stimulated to when the other cells are stimulated..
        
        this_bs = mean(cell_struct(i).('sta_traces')(:,opt.bs_frame_range,:),2);
        this_active = mean(cell_struct(i).('sta_traces')(:,frames_to_avg,:),2);
        labels = [ones(1,length(this_bs)),2.*ones(1,length(this_active))]';
        scores = [this_bs;this_active];
        [~,~,~,this_auc] = perfcurve(labels,scores,2);
        
        % shuffle to get zscore auc
        shuf_stim_auc = nan(1,num_shuf);
        for s = 1:num_shuf
            shuf_labels = labels(randperm(length(labels))');
            [~,~,~, shuf_stim_auc(s)] = perfcurve(shuf_labels,scores,1);
        end
        this_auc_zscore = (this_auc - mean(shuf_stim_auc))/std(shuf_stim_auc);
        
        cell_struct(i).photo_auc = this_auc;
        cell_struct(i).photo_auc_zscore = this_auc_zscore;
        
        % using a fixed auc threshold
        if  cell_struct(i).photo_auc_zscore > opt.N &&  cell_struct(i).sta_amp > opt.sta_amp_thresh
            cell_struct(i).is_photo = 1;
        else
             cell_struct(i).is_photo = 0;
        end
        
    else
        % add dummies
        cell_struct(i).(['sta_amp']) = 0;
        cell_struct(i).('is_photo') =0;
        cell_struct(i).photo_auc = 0.5;
        cell_struct(i).photo_auc_zscore = 0;
        
        warning(['roi ' num2str(i) ' has no photostim trial!'])
    end
    
end
    
end

% % sta of backgound component (for alignment check)
% [~,~,~,~,~,bg_sta_traces,bg_sta_trace] =...
%     make_sta_from_traces(backgroundC,sens_stim_frames ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
% 
% end

