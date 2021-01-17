function [opt] = init_opt_posthoc(opt)
opt.go_stim_types = [2,4]; % stim types with reward channel == 0 in pybehav
opt.nogo_stim_types = [1,3]; % stim types with reward channel == 0 in pybehav

opt.N = 1.5; % threshold for significant auc
opt.sta_pre_frames = 150;
opt.sta_post_frames = 120;
opt.sta_baseline_frames = 30; % relative to beginning of sta traces
opt.trial_length = 1+opt.sta_pre_frames+opt.sta_post_frames;

opt.sta_stim_frame = 90; % roughly when first contact occurs
opt.gocue_frame = 120; % relative to trial start
opt.stimon_frame = 90; % relative to trial start
opt.end_rw_frame = 180; % end of response window

% time window for first lick, 
% will not be taken as miss trial if first lick occur during withold window
opt.rw_win_end_sec = 5;
opt.withold_win_start_sec = 3;
opt.rw_start_sec = 4; %trials with first lick before this will not be used for analysis; 4 is beginning of response window in pybehav 

% frame indices relative to sta trace
opt.sta_gocue_frame = opt.sta_pre_frames;
opt.sta_trialon_frame = opt.sta_pre_frames-opt.gocue_frame;
opt.sta_avg_frames = [-10:1:-1]+opt.sta_gocue_frame; % different from OnlineProcTex to get photostim effect
opt.sta_peak_search_range =  [-60:1:0]+opt.sta_gocue_frame;

opt.withold_frames_adj = [60 120];

opt.offcell_avg_frames = 30;
opt.sta_amp_thresh = 1;
opt.frame_rate = 30;

opt.flag_use_peak = false; % if false then will use average to get roc

opt.correct_trial_only = false;

% select which online trace to use
opt.IF_USE_FILTC = true;



end

