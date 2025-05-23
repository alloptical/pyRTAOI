function [trials,trial_indices,sens_stim_frames,num_stim_type] = sort_trial_types_baseline(trials,caiman_data,FLAG_PYBEHAV_LOADED,IF_GO_NOGO,opt)

if isfield(caiman_data,'stim_frames_caiman')&&(~isempty(caiman_data.stim_frames_caiman))
    sens_stim_frames = caiman_data.sensory_stim_frames+caiman_data.t_init;
    discard_trial_idx = find(sens_stim_frames>caiman_data.t_cnm-opt.sta_post_frames-opt.sta_pre_frames);
    discard_trial_idx = unique(discard_trial_idx);
    sens_stim_frames(discard_trial_idx) = [];
    num_trials = numel(sens_stim_frames);

else
    sens_stim_frames = [];
    discard_trial_idx = [];
    num_trials = numel(trials.stim_type);
end

% discard sens stim frames beyond tottime (deal with senarios when movie is aborted early)

if ~isempty(opt.discard_trials_after)
    discard_trial_idx = [discard_trial_idx, opt.discard_trials_after:numel(sens_stim_frames)];
end
if ~isempty(opt.discard_trials_before)
    discard_trial_idx = [discard_trial_idx, 1:opt.discard_trials_before];
end
if FLAG_PYBEHAV_LOADED
%     trials.miss(trials.firstlick<opt.rw_start_sec) = 0; % discard early response trials from miss trials
%     trials.correct(trials.firstlick<opt.rw_start_sec) = 0; % discard early response trials 
%     trials.incorrect(trials.firstlick<opt.rw_start_sec) = 0; % discard early response trials
%     trials.fa(trials.firstlick<opt.rw_start_sec) = 0; % discard early response trials

    fds = fields(trials);
    try
        for f =1:numel(fds)
            trials.(fds{f})(discard_trial_idx) = [];
        end
    end
end
% trial type
trialOrder = caiman_data.trialOrder;
try
trialOrder(discard_trial_idx) = [];
catch
    disp('no discarded trials')
end
num_stim_type = length(unique(trialOrder)); % orientations/texture

if ~ FLAG_PYBEHAV_LOADED
    % make a dummy trial struct (setting all trials to correct)
    trials.correct = ones(size(sens_stim_frames));
    trials.incorrect = zeros(size(sens_stim_frames));
    trials.stim_type = trialOrder;
    trials.cheated = zeros(size(sens_stim_frames));
    trials.miss = zeros(size(sens_stim_frames));
end


%% Get trial types
trial_indices = struct(); % % get trialtype-outcome indices
for i = 1:num_stim_type
    trial_indices.(['stim_' num2str(i) '_correct']) = find(trials.correct==1&trials.stim_type==i&trials.cheated==0);
    trial_indices.(['stim_' num2str(i) '_incorrect']) = find(trials.incorrect==1&trials.stim_type==i&trials.cheated==0);
    trial_indices.(['stim_' num2str(i) '_correct'])(trial_indices.(['stim_' num2str(i) '_correct'])>num_trials)=[];
    trial_indices.(['stim_' num2str(i) '_incorrect'])(trial_indices.(['stim_' num2str(i) '_incorrect'])>num_trials)=[];
    if  IF_GO_NOGO
        if ~isempty(intersect(opt.go_stim_types,i))
            trial_indices.(['stim_' num2str(i) '_incorrect']) =  find(trials.miss==1&trials.stim_type==i&trials.cheated==0& ...
        (trials.firstlick>opt.rw_win_end_sec|trials.firstlick<opt.withold_win_start_sec|isnan(trials.firstlick)));
            trial_indices.(['stim_' num2str(i) '_correct']) = find(trials.stim_type==i&trials.cheated==0&...
                 (trials.firstlick<opt.rw_win_end_sec&trials.firstlick>opt.withold_win_start_sec&~isnan(trials.firstlick)));

        elseif  ~isempty(intersect(opt.nogo_stim_types,i))
            % no-go stim type correct trials should exclude trials with licks in
            % withold window (use fa as incorrect)
            trial_indices.(['stim_' num2str(i) '_correct']) = find(trials.correct==1&trials.fa==0&trials.stim_type==i&trials.cheated==0);
            trial_indices.(['stim_' num2str(i) '_incorrect']) = find(trials.fa==1&trials.stim_type==i&trials.cheated==0);
        end
    end
end

try
% catch trials
    trial_indices.(['stim_5_miss']) = find((trials.stim_type==3|trials.stim_type==4)&trials.miss==1); 
    trial_indices.(['stim_5_lick']) = find((trials.stim_type==3|trials.stim_type==4)&trials.miss==0); 
end

end

