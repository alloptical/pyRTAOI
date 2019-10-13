function [test_decod_struct] = detect_error_trials(test_decod_struct,this_proj_struct,thresh_struct,lick_struct,photostim_frames,opt,varargin)
% use projection on decoder for detecting
% error trials
% initiation
frame_rate = getOr(opt,'Fs',30);
fds = fields(this_proj_struct);
correct_fd ={fds{contains(fds, '_correct_')}};
incorrect_fd = {fds{contains(fds, '_incorrect_')}};
IF_FRAMEWISE = false;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'IF_FRAMEWISE')
        IF_FRAMEWISE = varargin{v+1};
    end
end
for stim_type = 1:2
    this_fd = ['stim' num2str(stim_type)];
%     this_disc_frame = cd_struct_train.(this_fd).(disc_frames_fd);

    this_correct_fd = correct_fd{contains(correct_fd, ['_' num2str(stim_type)])};
    this_incorrect_fd = incorrect_fd{contains(incorrect_fd, ['_' num2str(stim_type)])};
    correct_traces = this_proj_struct.(this_correct_fd); %[trials, frames]
    incorrect_traces =  this_proj_struct.(this_incorrect_fd);
    
    num_correct_trials = size(correct_traces,1);
    num_incorrect_trials = size(incorrect_traces,1);
    proc_correct_traces = correct_traces-thresh_struct.(this_fd);
    proc_incorrect_traces = incorrect_traces-thresh_struct.(this_fd);
    

    if ~isempty(lick_struct)
    for tr = 1:num_correct_trials
        proc_correct_traces(tr,lick_struct.frame.(['correct_stim' num2str(stim_type)])(tr)+1:end) = nan;
    end
    for tr = 1:num_incorrect_trials
        proc_incorrect_traces(tr,lick_struct.frame.(['incorrect_stim' num2str(stim_type)])(tr)+1:end) = nan;
    end
    end

    [hit_trial_idx,~] = ind2sub(size(proc_incorrect_traces(:,photostim_frames)),find(proc_incorrect_traces(:,photostim_frames)<0));
    hit_trial_idx = unique(hit_trial_idx);
    
    [fa_trial_idx,~] = ind2sub(size(proc_correct_traces(:,photostim_frames)),find(proc_correct_traces(:,photostim_frames)<0));
    fa_trial_idx = unique(fa_trial_idx);
    
    hit_rate = numel(hit_trial_idx)/num_incorrect_trials;
    fa_rate =  numel(fa_trial_idx)/num_correct_trials;
    precision = hit_rate;
    recall = numel(hit_trial_idx)/...
        (numel(hit_trial_idx)+numel(fa_trial_idx));
    fscore = 2*(precision*recall)/(precision+recall);

    % detection time
    this_projs = zeros(size(proc_incorrect_traces));
    this_projs(:,photostim_frames)= proc_incorrect_traces(:,photostim_frames);
    num_trials = size(this_projs,1);
    this_det_frames = nan(num_trials,1);
    for t = 1:length(this_det_frames)
        try
            this_det_frames(t) = find(this_projs(t,:)<0,1);
        end
    end
    this_det_times = this_det_frames./frame_rate;
    if ~isempty(lick_struct)
        this_lick_time = lick_struct.time.(['incorrect_stim' num2str(stim_type)]);
        this_lick_time(isnan(this_det_frames)>0)=nan;
    else
        this_lick_time = [];
    end
    % false alarm time
    this_projs = nan(size(proc_correct_traces));
    this_projs(:,photostim_frames)= proc_correct_traces(:,photostim_frames);
    
    num_trials = size(this_projs,1);
    this_fa_frames = nan(num_trials,1);
    for t = 1:length(this_fa_frames)
        try
            this_fa_frames(t) = find(this_projs(t,:)<0,1);
        end
    end
    this_fa_times = this_fa_frames./frame_rate;

    
    test_decod_struct.stim_frame.(this_correct_fd) = this_fa_frames;
    test_decod_struct.stim_frame.(this_incorrect_fd) = this_det_frames;

    test_decod_struct.(this_fd).hit_rate = hit_rate;
    test_decod_struct.(this_fd).fa_rate = fa_rate;
    test_decod_struct.(this_fd).fscore = fscore;
    test_decod_struct.(this_fd).this_det_times = this_det_times;
    test_decod_struct.(this_fd).this_fa_times = this_fa_times;
    test_decod_struct.(this_fd).this_lick_time = this_lick_time;
    test_decod_struct.(this_fd).photo_enable_frame = photostim_frames;
    test_decod_struct.(this_fd).delay_time = nanmean(this_lick_time) - nanmean(this_det_times);
end


end

