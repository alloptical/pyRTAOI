function [trial_indices] = sort_trial_types_condition_incl_earlylick(trials,opt)
% trial types for go-nogo condition
% pyrtaoi stim: [1 2 3 4 3 4 5 5 1 2]; % trialOrder in caiman_data
% pyrtaoi var : [1 1 1 1 1 1 1 1 2 2]; % var = 2 are dummy closed-loop trials
% pybehav stim: [1 2 3 3 4 4 3 4 1 2];
% texture:      [1 2 3 3 3 3 3 3 1 2]; % by pybehav
% reward:       [0 2 0 0 2 2 2 0 0 2]; % by pybehav, 2 is go, 0 is no-go
% reward:       [1 2 1 1 2 2 1 2 1 2]; % by pybehav
% target:       [1 2 1 2 1 2 0 0 1 2]; % by pyrtaoi
num_trial_types = 10;
trial_indices = struct(); % % get trialtype-outcome indices
all_pyrtaoi_stimtype =      [1 2 3 4 3 4 5 5 1 2]; % will be used as variation below
all_pybehav_stimtype =      [1 2 3 3 4 4 3 4 1 2];
all_pyrtaoi_var      =      [1 1 1 1 1 1 1 1 2 2];% this gives dummy photostim trials
if opt.IF_GO_NOGO
    all_required_response = [0 2 0 0 2 2 0 2 0 2];
else
    all_required_response = [1 2 1 1 2 2 1 2 1 2];    
end
fds = fields(trials);
% discard first 10 trials for baselining 
% for f = 1:numel(fds)
% trials.(fds{f})(1:10) = 0;
% end

% excluding trials with constent licks/playing with ports
exclud_idx = find(trials.firstlick<opt.withold_win_start_sec&trials.fa==1);
for f = 1:numel(fds)
trials.(fds{f})(exclud_idx) = 0;
end


if ~isfield(trials,'oppo_photostim')
    trials.oppo_photostim = zeros(size(trials.photostim));
end

for v = 1:num_trial_types
    this_rtaoi_stim = all_pyrtaoi_stimtype(v);
    this_pybehav_stim = all_pybehav_stimtype(v);
    
    trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_correct'  ]) = find(trials.firstresponse==all_required_response(v)&trials.stim_type==this_pybehav_stim&trials.trialOrder == this_rtaoi_stim&trials.cheated==0);
    trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_incorrect' ]) = find(trials.correct==0&trials.firstresponse~=all_required_response(v)&trials.firstresponse>0&trials.stim_type==this_pybehav_stim&trials.trialOrder == this_rtaoi_stim&trials.cheated==0);
    trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_miss' ]) = find(trials.miss==1&trials.stim_type==this_pybehav_stim&trials.trialOrder == this_rtaoi_stim &trials.cheated==0& (trials.firstlick>opt.rw_win_end_sec|trials.firstlick<opt.withold_win_start_sec|isnan(trials.firstlick)));
    
    if opt.IF_GO_NOGO
        if ~isempty(intersect(opt.go_stim_types,this_pybehav_stim))
            trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_incorrect' ]) =  find(trials.miss==1&trials.stim_type==this_pybehav_stim&trials.cheated==0&trials.trialOrder == this_rtaoi_stim& ...
                (trials.firstlick>opt.rw_win_end_sec|trials.firstlick<opt.withold_win_start_sec|isnan(trials.firstlick)));
            trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_correct']) = find(trials.stim_type==this_pybehav_stim&trials.cheated==0&trials.trialOrder == this_rtaoi_stim&...
                (trials.firstlick<opt.rw_win_end_sec&trials.firstlick>opt.withold_win_start_sec&~isnan(trials.firstlick)));
            
        elseif ~isempty(intersect(opt.nogo_stim_types,this_pybehav_stim))
            trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_incorrect' ]) = find((trials.fa==1|trials.incorrect==1)&trials.stim_type==this_pybehav_stim&trials.trialOrder == this_rtaoi_stim&trials.cheated==0);
            trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_correct'  ]) = find(trials.correct==1&trials.fa==0&trials.stim_type==this_pybehav_stim&trials.trialOrder == this_rtaoi_stim&trials.cheated==0);
        end
    elseif this_pybehav_stim>2 % deal with catch trials in 2afc sessions (in case pybehav reward settings were wrong)
        trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_correct'  ]) = find(trials.miss<1&trials.firstresponse==all_required_response(v)&trials.stim_type==this_pybehav_stim&trials.trialOrder == this_rtaoi_stim&trials.cheated==0);
        trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_incorrect' ]) = find(trials.miss<1&trials.firstresponse>0&trials.firstresponse~=all_required_response(v)&trials.stim_type==this_pybehav_stim&trials.trialOrder == this_rtaoi_stim&trials.cheated==0);
        trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_miss' ]) = find(trials.firstresponse ==0&trials.stim_type==this_pybehav_stim&trials.trialOrder == this_rtaoi_stim &trials.cheated==0& (trials.firstlick>opt.rw_win_end_sec|trials.firstlick<opt.withold_win_start_sec|isnan(trials.firstlick)));
    end
    this_trials = [ trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_correct'  ]), trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_incorrect'  ])];
    trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_photostim'  ]) = intersect(this_trials,find(trials.stim_type==this_pybehav_stim&trials.trialOrder == this_rtaoi_stim & trials.photostim ==1&trials.trialVar==1&trials.cheated==0));
    trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_nonphotostim'  ]) = intersect(this_trials,find(trials.stim_type==this_pybehav_stim&trials.trialOrder == this_rtaoi_stim & trials.photostim ==0&trials.oppo_photostim ==0&trials.cheated==0));
    trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_dummyphotostim'  ]) =intersect(this_trials, find(trials.stim_type==this_pybehav_stim&trials.trialOrder == this_rtaoi_stim & trials.photostim ==1&trials.trialVar==2&trials.cheated==0));
    try
        trial_indices.(['stim_' num2str(this_pybehav_stim) '_var_' num2str(this_rtaoi_stim) '_oppophotostim'  ]) = intersect(this_trials,find(trials.stim_type==this_pybehav_stim&trials.trialOrder == this_rtaoi_stim & trials.oppo_photostim ==1&trials.cheated==0));       
    end
end
% sort catch trials
% using stim_5 as catch; var as target ensemble
trial_indices.('stim_5_photostim_1') =  find(trials.trialOrder==3&trials.photostim == 1);
trial_indices.('stim_5_photostim_2') =  find(trials.trialOrder==4&trials.photostim == 1);
trial_indices.('stim_5_nonphotostim') =  find(trials.trialOrder==5&trials.photostim == 0);


% performance in photo and nonphoto tirals
outcomes = {'correct','incorrect'};
photo_types = {'photostim','nonphotostim','dummyphotostim'};
stim_types = [1,2];
for s = 1:numel(stim_types)
    ss = stim_types(s);
    for o = 1:numel(outcomes)
        oo = outcomes{o};
        for p = 1:numel(photo_types)
            pp = photo_types{p};
            trial_indices.(['stim_' num2str(ss) '_' pp '_' oo]) = intersect(trial_indices.(['stim_' num2str(ss) '_var_' num2str(ss) '_' pp]),trial_indices.(['stim_' num2str(ss) '_var_' num2str(ss) '_' oo]));
        end
    end
end

% catch trials lick (go/no-go)
trial_indices.('stim_5_photostim_1_lick') = intersect([trial_indices.('stim_3_var_3_incorrect'),trial_indices.('stim_4_var_3_correct')],trial_indices.('stim_5_photostim_1'));
trial_indices.('stim_5_photostim_1_nolick') = intersect([trial_indices.('stim_3_var_3_correct'),trial_indices.('stim_4_var_3_incorrect')],trial_indices.('stim_5_photostim_1'));
trial_indices.('stim_5_photostim_2_lick') = intersect([trial_indices.('stim_3_var_4_incorrect'),trial_indices.('stim_4_var_4_correct')],trial_indices.('stim_5_photostim_2'));
trial_indices.('stim_5_photostim_2_nolick') = intersect([trial_indices.('stim_3_var_4_correct'),trial_indices.('stim_4_var_4_incorrect')],trial_indices.('stim_5_photostim_2'));
trial_indices.('stim_5_nonphotostim_lick') = intersect([trial_indices.('stim_3_var_5_incorrect'),trial_indices.('stim_4_var_5_correct')],trial_indices.('stim_5_nonphotostim'));
trial_indices.('stim_5_nonphotostim_nolick') = intersect([trial_indices.('stim_3_var_5_correct'),trial_indices.('stim_4_var_5_incorrect')],trial_indices.('stim_5_nonphotostim'));

try
trial_indices.('stim_1_photo') = unique([trial_indices.stim_1_photostim_correct,trial_indices.stim_1_photostim_incorrect]);
trial_indices.('stim_1_dummyphoto') = unique([trial_indices.stim_1_dummyphotostim_correct,trial_indices.stim_1_dummyphotostim_incorrect]);
trial_indices.('stim_1_nonphoto') = unique([trial_indices.stim_1_nonphotostim_correct,trial_indices.stim_1_nonphotostim_incorrect]);

trial_indices.('stim_2_photo') = unique([trial_indices.stim_2_photostim_correct,trial_indices.stim_2_photostim_incorrect]);
trial_indices.('stim_2_dummyphoto') = unique([trial_indices.stim_2_dummyphotostim_correct,trial_indices.stim_2_dummyphotostim_incorrect]);
trial_indices.('stim_2_nonphoto') = unique([trial_indices.stim_2_nonphotostim_correct,trial_indices.stim_2_nonphotostim_incorrect]);



end
end

