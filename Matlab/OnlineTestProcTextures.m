% analyse conditioning session in texture task

%% add path - change this for rig
clear all
close all
clc

% ZZ PC
% addpath(genpath('Y:\zzhang\Python\pyRTAOI\Matlab'));
% cd('Y:\zzhang\Python\pyRTAOI\Matlab');

% BRUKER1
 matlab_set_paths_zz

%%  parameters - CHANGE THIS
crop_num_trials = 325; % specify number of trials recorded if aborted half way
IF_CONTROL_SESSION = false;
FLAG_PAQ_LOADED = false;
disp('CHECK SETTINGS BEFORE CONTINEU!')
keyboard
IF_GO_NOGO = false;
IF_USE_PYRTAOI_STIMTYPE = true; % for condition session this will be different from pybehav
opt.IF_GO_NOGO = IF_GO_NOGO;
pyrtaoi_condition_types = [1,2,3,4]; % for checking files only. pyrtaoi will use [1,2] for closed-loop stim_type; 3,4 for photstim tex1 and tex2 ensembles in catch trials, and 5 for catch trial without photostim
pybehav_condition_types = [1,2,3,4];
photo_ensemble_types =    [1,2,3,4]; % catch trial photo types are offset from condition type by 1
easy_trial_idx = 1:10;
% init opt
[opt] = init_opt_posthoc(opt);
% init color
[trial_color] = online_tex_init_color();

%% load CAIMAN data
[caiman_file,caiman_path] = uigetfile('*.mat','Select caiman data');
caiman_data = load(fullfile(caiman_path,caiman_file));
disp(['Loaded file :',fullfile(caiman_path,caiman_file)])
opt.ds_factor = caiman_data.ds_factor;
try
    opt.photo_enable_frame = double(caiman_data.offsetFrames);
    opt.sta_stimon_frame = opt.photo_enable_frame+ opt.sta_baseline_frames;
catch
    opt.photo_enable_frame = [];
    opt.sta_stimon_frame = [];
    warning('no photostim for the session loaded')
end
% %% check processing time per frame
% figure
% plot(caiman_data.tottime)
%% load Pybehavior data
try
    [pb_file,pb_path] = uigetfile('*.mat','Select pybehavior data',caiman_path);
    behavior_data =  load(fullfile(pb_path,pb_file)) ;
    disp(['Loaded file :',fullfile(pb_path,pb_file)])
    
    FLAG_PYBEHAV_LOADED = true;
catch
    FLAG_PYBEHAV_LOADED = false;
end
%% load paq file - optional, check if photostim trial indices are correct
try
[paq_file,paq_path] = uigetfile('*.paq','Select paq file',caiman_path);
[paq_photo_frames,paq_trial_frames,running_speed,aom_volt] = read_paq_file(fullfile(paq_path,paq_file));
[paq_photostim_trial_idx] = get_trials_with_photostim( paq_trial_frames, paq_photo_frames );
disp(['Loaded file :',fullfile(paq_path,paq_file)])
FLAG_PAQ_LOADED = true;
catch
    warning('loading paq error')
end

%% load config file: baseline output struct (file named as 'OutputParams')
[baseline_file,baseline_path] = uigetfile('*.mat','Select baseline OutputParams',caiman_path);
baseline_output = load(fullfile(baseline_path,baseline_file));
disp(['Loaded file :',fullfile(baseline_path,baseline_file)])
decod_struct = baseline_output.output;
norm_weights = decod_struct.trigger_weights;
norm_thresh = decod_struct.trigger_thresh;
trigger_idx = decod_struct.trigger_idx;
target_idx = unique(decod_struct.target_idx);
target_ensembles = decod_struct.target_ensembles;
thresh_sd = decod_struct.thresh_sd;
%% config save path
save_path = [caiman_path filesep 'analysis_files'];
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

fig_save_path = [caiman_path filesep 'figures'];
if ~exist(fig_save_path, 'dir')
    mkdir(fig_save_path)
end
% setup save path and name
opt.output_path = caiman_path;
opt.exp_name = strrep(caiman_file,'.mat','proc');
disp(['analysis files savepath:' save_path])
%% organise data (generate plots for sanity check)
tot_num_trials = min([crop_num_trials,length(caiman_data.trialOrder),numel(behavior_data.results)]);

% get sensory frames
if(~isempty(caiman_data.stim_frames_caiman))
    sens_stim_frames = caiman_data.sensory_stim_frames+caiman_data.t_init;
else
    sens_stim_frames = [];
end

% discard sens stim frames beyond tottime (deal with senarios when movie is aborted early)
% note that here sen_stim_frames are trial-trigger frames to pybehav
sens_stim_frames(sens_stim_frames>caiman_data.t_cnm-opt.sta_post_frames) = [];
tot_num_trials = min([tot_num_trials,length(sens_stim_frames)]);
sens_stim_frames = sens_stim_frames(1:tot_num_trials);
num_trials = length(sens_stim_frames);


% check if trial order matches in behavior and caiman file
if FLAG_PYBEHAV_LOADED
    [trials,odd_trial_idx] = make_trials_struct(behavior_data);
    % discard trials after num_trials
    trials = structfun(@(x)x(1:num_trials),trials,'UniformOutput',false);
    
    % discard early response trials
    trials.miss(trials.firstlick<opt.rw_start_sec) = 0; 
    trials.correct(trials.firstlick<opt.rw_start_sec) = 0; 
    trials.incorrect(trials.firstlick<opt.rw_start_sec) = 0; 
    trials.fa(trials.firstlick<opt.rw_start_sec) = 0;
end
% pyrtaoi trial type
trialOrder = caiman_data.trialOrder(1:num_trials); %
try
trialVar = caiman_data.trialVar(1:num_trials);
catch
    trialVar = ones(1,num_trials);
end
num_stim_type = length(unique(trialOrder)); % orientations/texture
trials.trialOrder = trialOrder;
trials.trialVar = trialVar;


% only get correct trials
if opt.correct_trial_only && FLAG_PYBEHAV_LOADED
    sens_stim_frames = sens_stim_frames(trials.correct==1);
    trialOrder = trialOrder(trials.correct==1);
elseif ~ FLAG_PYBEHAV_LOADED
    % make a dummy trial struct (setting all trials to correct)
    trials.correct = ones(size(sens_stim_frames));
    trials.incorrect = zeros(size(sens_stim_frames));
    trials.stim_type = trialOrder;
    trials.stim_var = ones(size(trialOrder));
    trials.miss = zeros(size(sens_stim_frames));
end


% get trials with photostims
[photostim_trial_idx,num_photo_per_trial] = get_trials_with_photostim( caiman_data.sensory_stim_frames, caiman_data.online_photo_frames );
trials.photostim = zeros(1,num_trials);
trials.photostim(photostim_trial_idx(photostim_trial_idx<=num_trials))=1;

% get trials with photostim of the opposite ensemble
try
[oppo_photostim_trial_idx,num_photo_per_trial] = get_trials_with_photostim( caiman_data.sensory_stim_frames, caiman_data.online_oppo_photo_frames );
trials.oppo_photostim = zeros(1,num_trials);
trials.oppo_photostim(oppo_photostim_trial_idx(oppo_photostim_trial_idx<=num_trials))=1;
catch
    warning('no oppo photostim found')
end

% check if photostim trigger sent 
if FLAG_PAQ_LOADED
   disp(['photostim missmatch trial idx:' num2str(setdiff(photostim_trial_idx,paq_photostim_trial_idx))]);
%    figure; scatter(paq_photo_frames,ones(size(paq_photo_frames)),'MarkerEdgeColor','r'); hold on
%    scatter(caiman_data.online_photo_frames,ones(size(caiman_data.online_photo_frames)),'MarkerEdgeColor','black'); ylim([0.5,1.5])
end

if ~isempty( caiman_data.online_photo_frames)
    photo_stim_frames =  caiman_data.online_photo_frames + caiman_data.t_init;
    photo_stim_frames(photo_stim_frames>sens_stim_frames(end)+opt.trial_length)=[];
    photo_trial_idx = find(trials.photostim==1);
    dummy_trial_idx = find(trials.trialVar==2);
    [dummy_photostim_trial_idx,idx]=intersect(photo_trial_idx,dummy_trial_idx);
    true_photo_trial_idx = setdiff(photo_trial_idx,dummy_trial_idx);
    dummy_photo_stim_frames = photo_stim_frames(idx);
    photo_stim_frames = setdiff(photo_stim_frames,dummy_photo_stim_frames);
    try
        oppo_photo_frames = caiman_data.online_oppo_photo_frames + caiman_data.t_init;
    catch
        oppo_photo_frames = [];
        oppo_trial_idx = find(trials.oppo_photostim==1);
    end
else
    photo_stim_frames = [];
    disp('no photostim found')
end

%% check if frames with photostim were blocked by caiman
if FLAG_PAQ_LOADED
    movie_photo_stim_frames = photo_stim_frames - caiman_data.t_init; % for checking with raw movie
    movie_frames_skipped = caiman_data.frames_skipped;
    
    check_skip_trace = zeros(1,movie_frames_skipped(end));
    check_skip_trace(movie_frames_skipped+1) = 1;
    check_aom_volt = aom_volt(1:length(check_skip_trace));
    figure; hold on
    plot(check_skip_trace);
    stem(movie_photo_stim_frames,2.*ones(size(movie_photo_stim_frames)),'color','r')
    stem(paq_photo_frames,1.5.*ones(size(paq_photo_frames)))
    plot(check_aom_volt)
    
    try
        movie_dummy_photo_stim_frames = dummy_photo_stim_frames - caiman_data.t_init;
        stem(movie_dummy_photo_stim_frames,2.*ones(size(movie_dummy_photo_stim_frames)),'color',tint([1,0, 0],.5))
       

        movie_oppo_photo_frames = oppo_photo_frames - caiman_data.t_init; % for checking with raw movie
        stem(movie_oppo_photo_frames,2.*ones(size(movie_oppo_photo_frames)),'color','black')      
    end

    % check aom volt in photo and dummuy photo trials
    [~,~,~,aom_dummy,~,~,aom_dummy_avg] = make_sta_from_traces(check_aom_volt,dummy_photo_stim_frames- caiman_data.t_init,30,150,1:opt.sta_baseline_frames);
    [~,~,~,aom_photo,~,~,aom_photo_avg] = make_sta_from_traces(check_aom_volt,photo_stim_frames- caiman_data.t_init,30,150,1:opt.sta_baseline_frames);

    figure('name','aom volt check');subplot(1,2,1); plot(aom_dummy','color',[.5 .5 .5]);ylim([0 0.5]); title('control photo trials')
    subplot(1,2,2);plot(aom_photo','color','r');ylim([0 0.5]); title('photo trials')
       
    % check caiman recorded photo frame vs actual photo frame
    pre_frames = 30;
    check_photo_trace = zeros(1,movie_frames_skipped(end));
    check_photo_trace(paq_photo_frames) = 1;
    [~,~,~,paq_photo] = make_sta_from_traces(check_photo_trace,movie_photo_stim_frames,pre_frames,150,1:opt.sta_baseline_frames);
    check_photo_trace = zeros(1,movie_frames_skipped(end));
    check_photo_trace(paq_photo_frames) = 1;
    [~,~,~,paq_dummyphoto] = make_sta_from_traces(check_photo_trace,movie_dummy_photo_stim_frames,pre_frames,150,1:opt.sta_baseline_frames);
    check_photo_trace = zeros(1,movie_frames_skipped(end));
    check_photo_trace(paq_photo_frames) = 1;
    [~,~,~,paq_oppophoto] = make_sta_from_traces(check_photo_trace,movie_oppo_photo_frames,pre_frames,150,1:opt.sta_baseline_frames);
    
    figure('name','photostim delay check','position',[200 200 1200 500]);
    subplot(1,4,1);hold on;imagesc(paq_photo); colormap('gray')
    plot([30 30],ylim,'color','r');xlabel('Frames'); title('Photostim triggers'); axis square
    subplot(1,4,2);hold on;imagesc(paq_dummyphoto); colormap('gray')
    plot([30 30],ylim,'color','r');xlabel('Frames'); title('Dummy Photostim triggers'); axis square   
    subplot(1,4,3);hold on;imagesc(paq_oppophoto); colormap('gray')
    plot([30 30],ylim,'color','r');xlabel('Frames'); title('Oppo Photostim triggers'); axis square 
    
    % get trials with very long delay 
    max_photo_delay = 10; % frames
    binwidth = 3;
    [~,photo_delays] = find(paq_photo == 1); photo_delays = photo_delays-pre_frames;
    [~,dummyphoto_delays] = find(paq_dummyphoto == 1); dummyphoto_delays = dummyphoto_delays-pre_frames;
    [~,oppophoto_delays] = find(paq_oppophoto == 1); oppophoto_delays = oppophoto_delays-pre_frames;
     subplot(1,4,4); hold on
     histogram(photo_delays,'displaystyle','stairs','edgecolor','r','binwidth',binwidth)
     histogram(dummyphoto_delays,'displaystyle','stairs','edgecolor', trial_color.dummyphoto,'binwidth',binwidth)
     histogram(oppophoto_delays,'displaystyle','stairs','edgecolor',[0 0 0],'binwidth',binwidth); axis square
     xlabel('Delay frames'); ylabel('Num. trials')
     delayed_trials_idx = sort([true_photo_trial_idx(photo_delays>max_photo_delay) dummy_photostim_trial_idx(dummyphoto_delays>max_photo_delay) oppo_photostim_trial_idx(oppophoto_delays>max_photo_delay)]);
     delayed_phototrial_idx = sort([true_photo_trial_idx(photo_delays>max_photo_delay) oppo_photostim_trial_idx(oppophoto_delays>max_photo_delay)]);
     disp(['delayed trial idx:' num2str(delayed_trials_idx)])
     if ~isempty(delayed_phototrial_idx) % put these trials into cheated
         first_bug_trial = delayed_trials_idx(1);
         disp(['First delayed trial:' num2str(delayed_trials_idx(1))])
         suptitle(['First delayed trial:' num2str(delayed_trials_idx(1)) ', stim type:' num2str(trials.stim_type(first_bug_trial)) 'Var' num2str(trials.trialVar(first_bug_trial))...
             ' photo:'  num2str(trials.photostim(first_bug_trial)) ' fa:' num2str(trials.fa(first_bug_trial))])
         
         trials.cheated(delayed_phototrial_idx(delayed_phototrial_idx<tot_num_trials)) = 1;
     end
     export_fig([fig_save_path filesep 'PhotodelayCheck_' strrep(caiman_file,'.mat','')],'-png')

end
disp('got trial struct')
%% get contineous miss or licking trials
opt.type_color = [trial_color.('stim1');trial_color.('stim2')];
opt.trial_color = trial_color;
min_num_active_licks = 15; 
trials.active_licks = (~isnan(trials.firstlick))&(~trials.cheated);
nolicking_period = pt_continuousabove(1-trials.active_licks,0,0.1,min_num_active_licks,100000,0);
num_nolicking_period = size(nolicking_period,1);
disp(['Num. consecutive no response periods:' num2str(num_nolicking_period)])

figure
plot(trials.firstlick)
hold on;plot(trials.active_licks);legend('lick time','responded')
disp([' excluding the last no-licking period, OR QUIT DEBUGGING AND SKIP THIS'])
if num_nolicking_period>0
    keyboard
    last_no_licking_trials = nolicking_period(end,1):1:num_trials;
    trials = structfun(@(x)x(setdiff(1:numel(x),last_no_licking_trials+5)),trials,'UniformOutput',false);   
end
tot_num_trials = length(trials.stim_type);
sens_stim_frames = sens_stim_frames(1:tot_num_trials);

%% get trial indices
% trial types for go-nogo condition
% pyrtaoi stim: [1 2 3 4 3 4 5 5 1 2]; % trialOrder in caiman_data
% pyrtaoi var : [1 1 1 1 1 1 1 1 2 2]; % var = 2 are dummy closed-loop trials
% pybehav stim: [1 2 3 3 4 4 3 4 1 2];
% texture:      [1 2 3 3 3 3 3 3 1 2]; % by pybehav
% reward:       [0 2 0 0 2 2 2 0 0 2]; % by pybehav, 2 is go, 0 is no-go
% target:       [1 2 1 2 1 2 0 0 1 2]; % by pyrtaoi

[trial_indices] = sort_trial_types_condition(trials,opt);
if IF_CONTROL_SESSION
    trial_indices = sort_trial_types_control(trials,opt);
end
%exclude odd trials if any
trial_indices = structfun(@(x)setdiff(x,odd_trial_idx),trial_indices,'UniformOutput',false);
disp('sorted trial indices')
%% quantify behavior
% including early licks
[all_trial_indices] = sort_trial_types_condition(trials,opt);

figure('name','animal performance','position',[400 400 2400 600])
for dummy = 1
    if IF_GO_NOGO
        [dp,hr,fa,num_trials_struct] = get_gonogo_performance(trial_indices);
        
        subplot(1,4,1)
        fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(dp),'UniformOutput',false));
        scatter_cmp_conditions(dp,[],...
            1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',0);
        xtickangle(30)
        ylabel('d-prime')
        axis square
        subplot(1,4,2)
        fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(hr),'UniformOutput',false));
        scatter_cmp_conditions(hr,[],...
            1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',0);
        ylim([0,1])
        xtickangle(30)
        ylabel('Hit rate')
        axis square
        
        subplot(1,4,3)
        fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(fa),'UniformOutput',false));
        scatter_cmp_conditions(fa,[],...
            1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',0);
        ylim([0,1])
        xtickangle(30)
        ylabel('FA rate')
        axis square

    else
        [overall,tex1,tex2,num_trials_struct] = get_2afc_performance(all_trial_indices);
        
        subplot(1,4,1)
        fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(overall),'UniformOutput',false));
        scatter_cmp_conditions(overall,[],...
            1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',0);
        xtickangle(30)
        ylim([0,1])
        ylabel('Overall fraction correct')
        axis square
        subplot(1,4,2)
        fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(tex1),'UniformOutput',false));
        scatter_cmp_conditions(tex1,[],...
            1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',0);
        ylim([0,1])
        xtickangle(30)
        ylabel('Tex1 fraction correct')
        axis square
        
        subplot(1,4,3)
        fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(tex2),'UniformOutput',false));
        scatter_cmp_conditions(tex2,[],...
            1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',0);
        ylim([0,1])
        xtickangle(30)
        ylabel('Tex2 fraction correct')
        axis square
        
    end
end


subplot(1,4,4) 
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(num_trials_struct),'UniformOutput',false));
scatter_cmp_conditions(num_trials_struct,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',0);
xtickangle(30)
ylabel('Num. trials')
axis square

export_fig([fig_save_path filesep 'PerformanceSummary_' strrep(caiman_file,'.mat','')],'-png')


%% catch trial behavior
[pc_lick,lick1,lick2,num_catch_trials] = get_catchtrial_performance(trial_indices);
fd_colors = [trial_color.nonphoto;trial_color.photo;trial_color.photo];
figure('name','catch trial performance','position',[400 400 2400 600])

subplot(1,4,1)
scatter_cmp_conditions(pc_lick,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',0);
ylim([0,1])
xtickangle(30)
ylabel('Fraction lick')
axis square

subplot(1,4,2)
scatter_cmp_conditions(lick1,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',0);
ylim([0,1])
xtickangle(30)
ylabel('Fraction lick port1')
axis square

subplot(1,4,3)
scatter_cmp_conditions(lick2,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',0);
ylim([0,1])
xtickangle(30)
ylabel('Fraction lick port2')
axis square

subplot(1,4,4) 
scatter_cmp_conditions(num_catch_trials,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',0);
xtickangle(30)
ylabel('Num. trials')
axis square

export_fig([fig_save_path filesep 'CatchTrialsPerformanceSummary_' strrep(caiman_file,'.mat','')],'-png')

%% generate config file for control session
% change closed-loop trials to photostim or no photostim trials
% only matters for pyrtaoi
if ~IF_CONTROL_SESSION
control_trials = trials;
control_trials.trialOrder(control_trials.photostim==1&control_trials.trialOrder<=2) = 2+control_trials.trialOrder(control_trials.photostim==1&control_trials.trialOrder<=2);
control_trials.trialOrder(control_trials.photostim==0) = 5;

% same stim_type and stim_var pairs but shuffled
shuf_idx = randperm(num_trials-numel(easy_trial_idx));
control_idx = [easy_trial_idx,shuf_idx+numel(easy_trial_idx)];
control_trials = structfun(@(x)x(control_idx),control_trials,'UniformOutput',false);
pyrtaoi_seq = [control_trials.trialOrder; ones(1,num_trials)];
pybehav_seq = [control_trials.stim_type; ones(1,num_trials)];


% make trial sequence file for pyrtaoi and pybehav
save_time = datestr(now,'yyyymmdd_HHMM');
pyrtaoi_file_save_name = [opt.exp_name '_RTAOiPyBehavior_CONTROL_' save_time];
dlmwrite([opt.output_path filesep pyrtaoi_file_save_name '.txt'],pyrtaoi_seq)
disp(['saved as:' opt.output_path filesep pyrtaoi_file_save_name '.txt'] )
pyrtaoi_file_save_name = [opt.exp_name '_PyBehavior_CONTROL_' save_time];
dlmwrite([opt.output_path filesep pyrtaoi_file_save_name '.txt'],pybehav_seq)
disp(['saved as:' opt.output_path filesep pyrtaoi_file_save_name '.txt'] )
end
%% get online trajectory and weights
online_weights = caiman_data.ROIw;
online_thresh = caiman_data.ROIsumThresh;
online_sd = caiman_data.sd_level;
online_bs = caiman_data.bs_level;
online_w = caiman_data.ROIw;
online_traj = caiman_data.online_traj;
disp('got online trajectory and weights')
%% make cnm data structure
cnm_struct = struct();
cnm_dims = double(caiman_data.cnm_dims);
cnm_image = reshape(caiman_data.cnm_b,cnm_dims);
cnm_A = full(caiman_data.cnm_A);
num_frames = caiman_data.t_cnm;
num_comp = size(cnm_A,2);
comp_shape =[cnm_dims(1),cnm_dims(2),num_comp];
cm = com(sparse(double(cnm_A)),cnm_dims(1),cnm_dims(2)); % center of mass


if ~isempty(caiman_data.frames_skipped)
    skip_frames = caiman_data.frames_skipped + caiman_data.t_init;
    tot_frames = num_frames + numel(skip_frames);
    caiman_frames = setdiff([1:tot_frames],skip_frames);
else
    caiman_frames = 1:num_frames;
    tot_frames = num_frames;
end

for i = 1:num_comp
    cnm_struct(i).shape = reshape(cnm_A(:,i),cnm_dims);
    cnm_struct(i).centroid = cm(i,:);
    %     cnm_struct(i).frame_added = find(cnm_struct(i).noisyC >0,1);
    
    % set skipped frames to nan then interpolate
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.online_C(i+1,1:num_frames); % use cnm_C(i,1:num_frames) for offline analysis!
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).onlineC_full = temp_trace;
    
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.cnm_C(i,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).deconvC_full = temp_trace;
    
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.noisyC(i+1,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).noisy_full = temp_trace;
    % use stimon (start of deflection) frame as 'stim frame' for sta traces
    cnm_struct(i).stim_frames = sens_stim_frames+opt.gocue_frame;
    
end

temp_trace = nan(1,tot_frames);
temp_trace(caiman_frames) =  caiman_data.noisyC(1,1:num_frames);
temp_trace = fillmissing(temp_trace,'linear');
backgroundC = temp_trace;
disp('made cnm_struct')

% only get accepted cells
accepted_idx = caiman_data.accepted_idx+1;
num_cells = numel(accepted_idx);


%% plot spatial components and save to cell struct
com_fov = zeros(cnm_dims);
accepted_com_fov = zeros(cnm_dims);
for i = 1:num_comp
    com_fov = com_fov+cnm_struct(i).shape;
end

for i = accepted_idx
    accepted_com_fov = com_fov+cnm_struct(i).shape;
end


cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = [colormap(lines);colormap(lines);colormap(lines)];
close

figure('name','fov','units','normalized','outerposition',[0 0 1 1])
subplot(1,3,1)
imagesc(com_fov)
colormap(gray)
axis square
title('Detected ROIs')

subplot(1,3,2)
imagesc(accepted_com_fov)
colormap(gray)
axis square
title('Accepted ROIs')


subplot(1,3,3)
[CC,jsf] = plot_contours(sparse(double(cnm_A)),cnm_image,cnm_plot_options,1,[],[],[1 1 1]);
colormap(gray)
axis square
title('GCaMP')
export_fig([fig_save_path filesep 'ROIs_' strrep(caiman_file,'.mat','.png')])

cell_struct = struct();
for i = 1:num_cells
    this_idx = accepted_idx(i);
    temp_coords = jsf(this_idx).coordinates;
    lin_idx = zeros(size(temp_coords,1),1);
    
    for t = 1:size(temp_coords,1)
        lin_idx(t) = sub2ind(cnm_dims,temp_coords(t,1),temp_coords(t,2));
    end
    cell_struct(i).contour = CC{this_idx};
    cell_struct(i).lin_coords = lin_idx;
    cell_struct(i).coordinates = jsf(this_idx).coordinates;
    cell_struct(i).pix_values = jsf(this_idx).values;
    cell_struct(i).centroid = jsf(this_idx).centroid;
    cell_struct(i).opsin_positive = 0;
    cell_struct(i).cnm_idx = this_idx;
    cell_struct(i).jsf = jsf(this_idx);
    
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.filt_C(i,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    
    cell_struct(i).filtC = temp_trace;
    try
    cell_struct(i).cnn_score = caiman_data.CNN_predictions(this_idx);
    catch
        cell_struct(i).cnn_score = 1;
    end

end
%% check trigger and target cell traces
plot_frame_range = double([caiman_data.t_init sens_stim_frames(end)+opt.trial_length]);
figure('name','online trigger cell traces','units','normalized','outerposition',[0 0 1 1]); hold on
plot_offset = 10;
cell_count = 1;

this_idx = trigger_idx;
xlim(plot_frame_range)

for i = 1:length(this_idx)
    ii = this_idx(i);
    cell_count = cell_count+1;
    plot(cell_struct(ii).filtC+cell_count*plot_offset,'black');
    text(plot_frame_range(1),cell_count*plot_offset,['Cell ' num2str(ii) ', W ' num2str(online_w(ii))], 'horizontalalignment','right', 'color','black')
    
end
set(gca,'ytick',[])
ylim([-plot_offset plot_offset*(cell_count+1)])
for i = 1:numel(sens_stim_frames)
    this_color = trial_color.(['stim' num2str(trials.stim_type(i) )]);
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.stimon_frame,ylim,'color',this_color) % trial-on + roughly first touch
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.gocue_frame,ylim,'color',this_color,'linestyle',':') % go-cue
end

for i = 1:numel(photo_stim_frames)
    plot([photo_stim_frames(i) photo_stim_frames(i)],ylim,'color',trial_color.photostim)
end

for i = 1:numel(dummy_photo_stim_frames)
    plot([dummy_photo_stim_frames(i) dummy_photo_stim_frames(i)],ylim,'color',trial_color.photostim,'linestyle',':')
end
hold on;plot(online_traj,'color',[.5 .5 .5])
text(plot_frame_range(1),0,['Pop. trajectory'], 'horizontalalignment','right', 'color','black')

export_fig([fig_save_path filesep 'TriggerTrace_' strrep(caiman_file,'.mat','')],'-pdf','-painters')

%% target full traces
figure('name','online target cell traces','units','normalized','outerposition',[0 0 1 1]); hold on
cell_count = 1;
for e = 1:2
    this_idx = target_ensembles{e};
    for i = 1:length(this_idx)
        ii = this_idx(i);
        cell_count = cell_count+1;
        plot(cell_struct(ii).filtC+cell_count*plot_offset,'black');
        text(double(caiman_data.t_init),cell_count*plot_offset,['Cell ' num2str(ii) ', W ' num2str(online_w(ii))], 'horizontalalignment','right', 'color',trial_color.(['stim' num2str(e)]))

    end
end
set(gca,'ytick',[])

for i = 1:numel(sens_stim_frames)
    this_color = trial_color.(['stim' num2str(trials.stim_type(i) )]);
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.stimon_frame,ylim,'color',this_color) % trial-on + roughly first touch
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.gocue_frame,ylim,'color',this_color,'linestyle',':') % go-cue
end

for i = 1:numel(photo_stim_frames)
    plot([photo_stim_frames(i) photo_stim_frames(i)],ylim,'color',trial_color.photostim)
end

for i = 1:numel(dummy_photo_stim_frames)
    plot([dummy_photo_stim_frames(i) dummy_photo_stim_frames(i)],ylim,'color',trial_color.photostim,'linestyle',':')
end

xlim(plot_frame_range)
ylim([0 (cell_count+1)*plot_offset])
export_fig([fig_save_path filesep 'TargetTrace_' strrep(caiman_file,'.mat','')],'-pdf','-painters')


%% plot full traces for all rois
figure('name','all ROI traces','units','normalized','outerposition',[0 0 1 1]); hold on
plot_offset = 5;
stim_cell_count = 1;
non_stim_cell_count = 1;

for i = 1:num_cells
    this_cell_trace = zscore(cnm_struct(cell_struct(i).cnm_idx).deconvC_full);
    plot(this_cell_trace+i*plot_offset,'color','black','linewidth',1.5)
    stim_cell_count = stim_cell_count+1;
end
xlabel('Frames')
ylabel('ROI index')
yticks([ 1:num_cells].*plot_offset)
yticklabels(1:num_cells)
xlim([caiman_data.t_init tot_frames])
ylim([0 num_cells].*plot_offset+5)

% background
plot(backgroundC,'color',[.5 .5 .5],'linestyle',':')

for i = 1:numel(sens_stim_frames)
    this_color = trial_color.(['stim' num2str(trials.stim_type(i) )]);
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.stimon_frame,ylim,'color',this_color) % trial-on + roughly first touch
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.gocue_frame,ylim,'color',this_color,'linestyle',':') % go-cue
end

for i = 1:numel(photo_stim_frames)
    plot([photo_stim_frames(i) photo_stim_frames(i)],ylim,'color',trial_color.photostim) % photostim with power
end

for i = 1:numel(dummy_photo_stim_frames)
    plot([dummy_photo_stim_frames(i) dummy_photo_stim_frames(i)],ylim,'color',trial_color.photostim,'linestyle',':') % galvo stim with no laser power
end
export_fig([fig_save_path filesep 'AllROITrace_' strrep(caiman_file,'.mat','.png')])

%% get cell STA
for i = 1:num_cells
    this_cell_trace = zscore(cnm_struct(cell_struct(i).cnm_idx).deconvC_full);
    this_online_trace = cell_struct(i).filtC; % online trace med filtered
    
    
    this_num_trials = numel(cnm_struct(cell_struct(i).cnm_idx).stim_frames );
    this_sens_stim_frames =  cnm_struct(cell_struct(i).cnm_idx).stim_frames;

    cell_struct(i).sta_amp = 0;
    cell_struct(i).sta_traces = [];
    cell_struct(i).sta_trace = [];

    
    if(this_num_trials>0)
        [~,~,~,cell_struct(i).sta_traces,~,~,cell_struct(i).sta_trace] = make_sta_from_traces(this_cell_trace,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        [~,~,~,cell_struct(i).raw_sta_traces,~,~,cell_struct(i).raw_sta_trace] = make_sta_from_traces(this_online_trace,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        
        cell_struct(i).sta_amp = mean(cell_struct(i).sta_trace(opt.sta_avg_frames));
        
    end
    cell_struct(i).sta_traces = cell_struct(i).sta_traces-min(cell_struct(i).sta_traces(:,1:opt.sta_baseline_frames),[],2);
    
end
disp('got cell_struct sta')

%% Normlalise traces to baseline (from first trials)
% get baseline and std from the first sets of easy trials
% - looks more consistent with traninig session to normalise this way
% doesnt make sense when sd is close to zero!! - just normalise to baseline
% in case intensity drifts...
cell_bs = cell2mat(arrayfun(@(x)mean(reshape(cell_struct(x).raw_sta_traces(easy_trial_idx,1:opt.sta_baseline_frames),[1,numel(easy_trial_idx)*opt.sta_baseline_frames]),2),1:num_cells,'UniformOutput', false));
cell_sd = cell2mat(arrayfun(@(x)std(reshape(cell_struct(x).raw_sta_traces(easy_trial_idx,:),[1,numel(easy_trial_idx)*opt.trial_length])),1:num_cells,'UniformOutput', false));
cell_mean = cell2mat(arrayfun(@(x)mean(cell_struct(x).raw_sta_trace(:)),1:num_cells,'UniformOutput', false));

for i = 1:num_cells
    cell_struct(i).raw_sta_traces = (cell_struct(i).raw_sta_traces - cell_bs(i));
    cell_struct(i).raw_sta_trace = (cell_struct(i).raw_sta_trace - cell_bs(i));
end

disp('normalised cell_struct raw_sta_traces')

%% STA of online trajectory
traj_struct = struct();
temp_trace = nan(1,tot_frames);
temp_trace(caiman_frames) = online_traj(1,1:num_frames);
online_traj = fillmissing(temp_trace,'linear');

[~,~,~,traj_struct.sta_traces,~,~,traj_struct.sta_trace] =...
    make_sta_from_traces(online_traj,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);

% sta of backgound component (for alignment check)
[~,~,~,~,~,bg_sta_traces,bg_sta_trace] =...
    make_sta_from_traces(backgroundC,this_sens_stim_frames ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);

%% Sort STAs by different trial types
trial_types = fields(trial_indices);
num_sta_traces = size(cell_struct(1).sta_traces,1);
raw_cell_struct = struct();

for i = 1:numel(trial_types)
    this_fd = trial_types{i};
    this_idx = trial_indices.(this_fd);
    this_idx = this_idx(this_idx<=min([tot_num_trials,num_sta_traces]));
%     this_idx = setdiff(this_idx,easy_trial_idx);
    if ~isempty(this_idx)
        for c = 1:num_cells
            cell_struct(c).(this_fd) = cell_struct(c).sta_traces( this_idx,:)';
            raw_cell_struct(c).(this_fd) = cell_struct(c).raw_sta_traces( this_idx,:)';
            traj_struct.(this_fd) = traj_struct.sta_traces(this_idx,:);
        end
    end
end
% normalise to trial baseline
fd_names = fields(trial_indices);
for f = 1:numel(fd_names)
    try
    for i = 1:num_cells
        cell_struct(i).(fd_names{f}) =  cell_struct(i).(fd_names{f}) - mean(cell_struct(i).(fd_names{f})(1:opt.sta_baseline_frames,:),1);
    end
    catch
       disp([fd_names{f} ' skipped'])
    end
end
disp('sorted cell_struct sta_traces')

%% Plot STAs traces for certain cell types
figure('name','condition sta traces','units','normalized','outerposition',[0 0 1 1])
stim_fds = {'stim_1_var_1_nonphotostim','stim_2_var_2_nonphotostim'...
    'stim_1_var_1_photostim','stim_2_var_2_photostim'}; % extreme stim types

peak_frame_range = opt.sta_peak_search_range;
avg_frame_range = opt.sta_avg_frames;


plot_cell_idx = sort(target_idx);
plot_num_cells = numel(plot_cell_idx);

num_plot_cols = 8;
num_plot_rows = ceil(plot_num_cells/num_plot_cols);
for ii = 1:plot_num_cells
    subtightplot(num_plot_rows,num_plot_cols,ii)
    i = plot_cell_idx(ii);
    hold on
    % correct trials
    plot(cell_struct(i).(stim_fds{1}),'color',trial_color.stim_1_correct,'linewidth',1)
    plot(cell_struct(i).(stim_fds{2}),'color',trial_color.stim_2_correct,'linewidth',1)
    
    % photostim trials
    plot(cell_struct(i).(stim_fds{3}),':','color',[0 0 0],'linewidth',1)
    plot(cell_struct(i).(stim_fds{4}),':','color',[.5 .5 .5],'linewidth',1)
    
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    axis square
    
    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')
    text(1,0.9,['W '  num2str(online_w(i),'%0.001f') ],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')
    
    box off
    % show post-hoc decoder weights
    try
        text(0.05,.7,['stim-w = '  num2str(stim_norm_weights(trigger_idx==i),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
        text(0.05,.6,['choice-w = '  num2str(choice_norm_weights(trigger_idx==i),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    end
    
    % mark time
    plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':')
    if ~ opt.flag_use_peak
        plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color','black','linewidth',2)
    else
        plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color','black','linewidth',2)
    end
    
    % mark target ensembles
    if( any(target_ensembles{1}==i))
        box on
        this_color = trial_color.(['correct_stim1' ]);
        set(gca,'XColor',this_color,'YColor',this_color,'linewidth',3)      
    end
    
    if( any(target_ensembles{2}==i))
        box on
        this_color = trial_color.(['correct_stim2' ]);
        set(gca,'XColor',this_color,'YColor',this_color,'linewidth',3)
        
    end
    
    
    % modify x ticks
    xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
suptitle(strrep(strrep(caiman_file,'.mat',''),'_',' '))
export_fig([fig_save_path filesep 'STATrace_' strrep(caiman_file,'.mat','.png')])
%% catch trial photostim sta traces
figure('name','condition sta traces','units','normalized','outerposition',[0 0 1 1])
stim_fds = {'stim_5_photostim_1','stim_5_photostim_2'...
    'stim_5_nonphotostim'}; % catch trials

avg_frame_range = opt.sta_avg_frames;

plot_cell_idx = sort(target_idx);
plot_num_cells = numel(plot_cell_idx);

num_plot_cols = 8;
num_plot_rows = ceil(plot_num_cells/num_plot_cols);
for ii = 1:plot_num_cells
    subtightplot(num_plot_rows,num_plot_cols,ii)
    i = plot_cell_idx(ii);
    hold on
    % photo trials
    plot(cell_struct(i).(stim_fds{1}),'color',trial_color.stim_1_correct,'linewidth',1)
    plot(cell_struct(i).(stim_fds{2}),'color',trial_color.stim_2_correct,'linewidth',1)
    
    % catch trials
    plot(cell_struct(i).(stim_fds{3}),'color',[.5 .5 .5],'linewidth',1)
    
    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    axis square
    
    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')
    text(1,0.9,['W '  num2str(online_w(i),'%0.001f') ],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')
    
    box off
    
    
    % mark time
    plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':')
    if ~ opt.flag_use_peak
        plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color','black','linewidth',2)
    else
        plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color','black','linewidth',2)
    end
    
    % mark target ensembles
    if( any(target_ensembles{1}==i))
        box on
        this_color = trial_color.(['correct_stim1' ]);
        set(gca,'XColor',this_color,'YColor',this_color,'linewidth',3)      
    end
    
    if( any(target_ensembles{2}==i))
        box on
        this_color = trial_color.(['correct_stim2' ]);
        set(gca,'XColor',this_color,'YColor',this_color,'linewidth',3)
        
    end
    
    % show post-hoc decoder weights
    try
        text(0.05,.7,['stim-w = '  num2str(stim_norm_weights(trigger_idx==i),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
        text(0.05,.6,['choice-w = '  num2str(choice_norm_weights(trigger_idx==i),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    end
    
    
    % modify x ticks
    xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
suptitle(strrep(strrep(caiman_file,'.mat',''),'_',' '))
export_fig([fig_save_path filesep 'CatchPhotoTrace_' strrep(caiman_file,'.mat','.png')])

%% get photostim sta amp
photo_ensembles = decod_struct.target_ensembles;
photo_types = pyrtaoi_condition_types;
catch_photo_types = [3,4,6,7,8,9]; % catch texture with photostim
catch_nonphoto_type = 5;   % catch texture without photostim

for i = 1:num_cells
    for e = 1:numel(pyrtaoi_condition_types)
        
        this_photo_type = pyrtaoi_condition_types(e);
        this_stim_type = pybehav_condition_types(e);
        this_photo_trials = find(trials.photostim==1&trials.stim_type == this_stim_type&trials.trialVar==1&trials.cheated==0&trials.trialOrder==this_photo_type);
        this_dummy_photo_trials = find(trials.photostim==1&trials.stim_type == this_stim_type&trials.trialOrder==this_photo_type&trials.trialVar==2&trials.cheated==0);
        this_oppo_photo_trials = find(trials.oppo_photostim==1&trials.stim_type == this_stim_type&trials.cheated==0);

        if ~isempty(intersect(this_stim_type,catch_photo_types))       
            % pyrtaoi defines stim type
            this_photo_trials = find(trials.photostim==1&trials.trialOrder==this_photo_type&trials.stim_type>2 &trials.trialVar==1&trials.cheated==0);
            this_control_trials = find(trials.photostim==0&trials.stim_type>2&trials.trialOrder==catch_nonphoto_type&trials.cheated==0);
        elseif ~isempty(this_dummy_photo_trials)         % for dummy photostim trials
            this_control_trials = this_dummy_photo_trials;
        else
            this_control_trials = find(trials.photostim==0&trials.stim_type == this_stim_type&trials.cheated==0);
        end
        
        cell_struct(i).(['sta_amp_nonphoto_' num2str(photo_types(e))]) = mean(mean(cell_struct(i).sta_traces(this_control_trials,opt.sta_avg_frames)));       
        cell_struct(i).(['sta_amp_photo_' num2str(photo_types(e))]) = mean(mean(cell_struct(i).sta_traces(this_photo_trials,opt.sta_avg_frames))); 
        cell_struct(i).(['sta_amp_diffphoto_' num2str(photo_types(e))]) = cell_struct(i).(['sta_amp_photo_' num2str(this_photo_type)]) -cell_struct(i).(['sta_amp_nonphoto_' num2str(this_photo_type)]);
        cell_struct(i).(['sta_trace_nonphoto_' num2str(photo_types(e))]) = cell_struct(i).sta_traces(this_control_trials,:);
        cell_struct(i).(['sta_trace_photo_' num2str(photo_types(e))]) = cell_struct(i).sta_traces(this_photo_trials,:);
        
        if ~isempty(this_oppo_photo_trials)
            cell_struct(i).(['sta_amp_oppophoto_' num2str(photo_types(e))]) = mean(mean(cell_struct(i).sta_traces(this_oppo_photo_trials,opt.sta_avg_frames)));
            cell_struct(i).(['sta_trace_oppophoto_' num2str(photo_types(e))]) = cell_struct(i).sta_traces(this_oppo_photo_trials,:);            
        end

    end
    
end
disp('got photo amp')
%% Plot photostim STA amp on FOV
max_abs_value = 5; % will saturate colorlut at this value if greater than
for e = 1:numel(pyrtaoi_condition_types)
    ee= photo_ensemble_types(e);
    figure('name',['photstim response on fov type' num2str(e)],'units','normalized','outerposition',[0 0 1 1]);
    plot_count = 1;
    ax = subplot(1,3,plot_count);
    if e==1
        [~,zlimit] = plot_value_in_rois( cell_struct, ['sta_amp_photo_' num2str(photo_types(e))],[256 256],ax,'IF_NORM_PIX',1,...
            'IF_CONTOUR',0,'IF_SHOW_OPSIN',0,'target_cell_idx',photo_ensembles{ee},'show_cell_idx',photo_ensembles{ee},'max_abs_value',max_abs_value);
        zlimit(end) = min(zlimit(end),max_abs_value);
        zlimit(1) = max(zlimit(1),-max_abs_value);

        plot_zlimit = zlimit;
    else
        plot_value_in_rois( cell_struct, ['sta_amp_photo_' num2str(photo_types(e))],[256 256],ax,'IF_NORM_PIX',1,...
            'IF_CONTOUR',0,'IF_SHOW_OPSIN',0,'target_cell_idx',photo_ensembles{ee},'zlimit',plot_zlimit,'show_cell_idx',photo_ensembles{ee});
    end
    title(['Stim type:' num2str(photo_types(e)) ', photo+'])
    plot_count= plot_count+1;
    
    ax = subplot(1,3,plot_count);
    plot_value_in_rois( cell_struct, ['sta_amp_nonphoto_' num2str(photo_types(e))],[256 256],ax,'IF_NORM_PIX',1,...
        'IF_CONTOUR',0,'IF_SHOW_OPSIN',0,'zlimit',plot_zlimit,'target_cell_idx',photo_ensembles{ee},'show_cell_idx',photo_ensembles{ee});
    plot_count= plot_count+1;
    title(['Stim type:' num2str(photo_types(e)) ', photo-'])
    
    ax = subplot(1,3,plot_count);
    plot_value_in_rois( cell_struct, ['sta_amp_diffphoto_' num2str(photo_types(e))],[256 256],ax,'IF_NORM_PIX',1,...
        'IF_CONTOUR',0,'IF_SHOW_OPSIN',0,'target_cell_idx',photo_ensembles{ee},'show_cell_idx',photo_ensembles{ee},'zlimit',[-3 3]);

    title(['Stim type:' num2str(photo_types(e)) ', diff'])
    export_fig([fig_save_path filesep 'PhotoSTA_FOV_Stim' num2str(e) strrep(caiman_file,'.mat','.png')])
    
end

%% Plot photostim STA traces
peak_frame_range = opt.sta_peak_search_range;
avg_frame_range = opt.sta_avg_frames;

plot_cell_idx = target_idx;
plot_num_cells = numel(plot_cell_idx);

num_plot_cols = 8;
num_plot_rows = ceil(plot_num_cells/num_plot_cols);
% raw traces
for photo_idx = 1:numel(pyrtaoi_condition_types)
    this_photo_type = pyrtaoi_condition_types(photo_idx);
    figure('name','condition sta traces','units','normalized','outerposition',[0 0 1 1])
    
    for ii = 1:plot_num_cells
        subtightplot(num_plot_rows,num_plot_cols,ii)
        i = plot_cell_idx(ii);
        hold on
        plot(cell_struct(i).(['sta_trace_nonphoto_' num2str(this_photo_type)'])','color',[.5 .5 .5],'linewidth',1)
        plot(cell_struct(i).(['sta_trace_photo_' num2str(this_photo_type)'])','color',[0 0 0],'linewidth',1)
        xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
        axis square
        
        text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')
      
        box off
    
        % mark time
        plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':')
        if ~ opt.flag_use_peak
            plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color','black','linewidth',2)
        else
            plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color','black','linewidth',2)
        end
        
        % mark trigger cell
        if( any(trigger_idx==i))
            box on
            set(gca,'linewidth',3)
            text(0.05,.8,['weight'  num2str(norm_weights(trigger_idx==i),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
        end
        
        % mark target cell in this ensemble
        if( any(photo_ensembles{photo_idx}==i))
            box on
            set(gca,'XColor','r','YColor','r','linewidth',2)
            
        end
        
        % show post-hoc decoder weights
        try
           text(0.05,.7,['stim-w = '  num2str(stim_norm_weights(trigger_idx==i),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
           text(0.05,.6,['choice-w = '  num2str(choice_norm_weights(trigger_idx==i),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
        end
        
        
        % modify x ticks
        xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
        xticks(xaxisvalues)
        xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
        
    end
    suptitle([strrep(strrep(caiman_file,'.mat',''),'_',' ') ' StimType: ' num2str(photo_idx)])
    export_fig([fig_save_path filesep 'Stim' num2str(photo_idx) '_STATrace_' strrep(caiman_file,'.mat','.png')])
end
% shaded traces
x_ticks =[0:1:opt.trial_length-1];
for photo_idx = 1:numel(pyrtaoi_condition_types)
    this_photo_type = pyrtaoi_condition_types(photo_idx);
    figure('name','condition sta traces','units','normalized','outerposition',[0 0 1 1])
    
    for ii = 1:plot_num_cells
        subtightplot(num_plot_rows,num_plot_cols,ii)
        i = plot_cell_idx(ii);
        hold on
        % photostimulated trials in black
        
        this_traces = cell_struct(i).(['sta_trace_photo_' num2str(this_photo_type)']);
        shadedErrorBar(x_ticks,mean(this_traces,1), std(this_traces,[],1),{'color',[0,0,0],'linewidth',2},0.1);
        
        
        this_traces = cell_struct(i).(['sta_trace_nonphoto_' num2str(this_photo_type)']);
        shadedErrorBar(x_ticks,mean(this_traces,1), std(this_traces,[],1),{'color',[0.5,0.5,0.5],'linewidth',2},0.1);
        
        
        xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
        axis square
        
        text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')
        
        
        box off
        
        
        % mark time
        plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':')
        if ~ opt.flag_use_peak
            plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color','black','linewidth',2)
        else
            plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color','black','linewidth',2)
        end
        
        % mark trigger cell
        if( any(trigger_idx==i))
            box on
            set(gca,'linewidth',3)
            text(0.05,.8,['weight'  num2str(norm_weights(trigger_idx==i),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
        end
        
        % mark target cell in this ensemble
        if( any(photo_ensembles{photo_idx}==i))
            box on
            set(gca,'XColor','r','YColor','r','linewidth',2)
            
        end
        
        try
           text(0.05,.7,['stim-w = '  num2str(stim_norm_weights(trigger_idx==i),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
           text(0.05,.6,['choice-w = '  num2str(choice_norm_weights(trigger_idx==i),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
        end
        
        % modify x ticks
        xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
        xticks(xaxisvalues)
        xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
        
    end
    suptitle([strrep(strrep(caiman_file,'.mat',''),'_',' ') ' StimType: ' num2str(photo_idx)])
    export_fig([fig_save_path filesep 'Stim' num2str(photo_idx) '_STATrace_' strrep(caiman_file,'.mat','.png')])
end

%% ensemble average traces
peak_frame_range = opt.sta_peak_search_range;
avg_frame_range = opt.sta_avg_frames;
figure('name','condition sta traces','units','normalized','outerposition',[0 0 .5 1])
% target response in tex1 and tex2 trials
ensemble_cell_indices = {target_ensembles{1},target_ensembles{2}};
ensemble_names = {'Tex1 targets','Tex2 targets'};

% target response in all trial types
ensemble_cell_indices = target_ensembles(photo_ensemble_types);
ensemble_names = cell(1,numel(pyrtaoi_condition_types));
ensemble_names(:) = {'Tar.'};

ylimit = [-1,4];
x_ticks =[0:1:opt.trial_length-1];
plot_count = 1;
num_plot_cols = numel(pyrtaoi_condition_types);
num_plot_rows = numel(ensemble_cell_indices);
for e = 1:numel(ensemble_cell_indices)
    plot_cell_idx = ensemble_cell_indices{e};
    
    for photo_idx = 1:numel(pyrtaoi_condition_types)
        this_photo_type = pyrtaoi_condition_types(photo_idx);
        subplot(num_plot_rows,num_plot_cols,plot_count); hold on;
        
        this_traces = {cell_struct(plot_cell_idx).(['sta_trace_nonphoto_' num2str(this_photo_type)])}';
        this_traces = cell2mat(arrayfun(@(x)permute(x{:},[3 1 2]),this_traces,'UniformOutput',false));
        this_traces = squeeze(mean(this_traces,1));
        if size(this_traces,2)>1
            shadedErrorBar(x_ticks,mean(this_traces,1), std(this_traces,[],1),{'color',[.5 .5 .5],'linewidth',2},0.1);
        end
        
        this_traces = {cell_struct(plot_cell_idx).(['sta_trace_photo_' num2str(this_photo_type)])}';
        this_traces = cell2mat(arrayfun(@(x)permute(x{:},[3 1 2]),this_traces,'UniformOutput',false));
        this_traces = squeeze(mean(this_traces,1));
        if size(this_traces,2)>1
            
            shadedErrorBar(x_ticks,mean(this_traces,1), std(this_traces,[],1),{'color',[0,0,0],'linewidth',2},0.1);
            
        end
        xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
        ylim(ylimit)
        axis square
        
        box off
        
        % mark time
        plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':')
        if ~ opt.flag_use_peak
            plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color','black','linewidth',2)
        else
            plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color','black','linewidth',2)
        end
        
        
        % modify x ticks
        xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
        xticks(xaxisvalues)
        xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
        
        title([ensemble_names{e}, num2str(photo_ensemble_types(e)) ' TrialVar. ' num2str(this_photo_type)])
        plot_count = plot_count+1;
    end
    
end
export_fig([fig_save_path filesep 'TargetEnsemble_STATrace_' strrep(caiman_file,'.mat','.png')])
%%
figure('name','oppo condition sta traces','units','normalized','outerposition',[0 0 .3 1])
ensemble_cell_indices = {target_ensembles{1},target_ensembles{2}};
ensemble_names = {'Tex1 targets','Tex2 targets'};
plot_condition_types = [1,2];
ylimit = [-1,4];
plot_count = 1;
num_plot_cols = numel(plot_condition_types);
num_plot_rows = numel(ensemble_cell_indices);
for e = 1:numel(ensemble_cell_indices)
    plot_cell_idx = ensemble_cell_indices{e};    
    for photo_idx = 1:numel(plot_condition_types)
        
        this_photo_type = plot_condition_types(photo_idx);
        subplot(num_plot_rows,num_plot_cols,plot_count); hold on;
        try
        
        this_traces = {cell_struct(plot_cell_idx).(['sta_trace_nonphoto_' num2str(this_photo_type)])}';
        this_traces = cell2mat(arrayfun(@(x)permute(x{:},[3 1 2]),this_traces,'UniformOutput',false));
        this_traces = squeeze(mean(this_traces,1));
        if size(this_traces,2)>1
            shadedErrorBar(x_ticks,mean(this_traces,1), std(this_traces,[],1),{'color',[.5 .5 .5],'linewidth',2},0.1);
        end
        
        this_traces = {cell_struct(plot_cell_idx).(['sta_trace_oppophoto_' num2str(this_photo_type)])}';
        this_traces = cell2mat(arrayfun(@(x)permute(x{:},[3 1 2]),this_traces,'UniformOutput',false));
        this_traces = squeeze(mean(this_traces,1));
        if size(this_traces,2)>1
            
            shadedErrorBar(x_ticks,mean(this_traces,1), std(this_traces,[],1),{'color',[0,0,0],'linewidth',2},0.1);
            
        end
        xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
        ylim(ylimit)
        axis square
        
        box off
        
        % mark time
        plot([1,1].*opt.sta_gocue_frame, ylim,'color','black','linestyle',':')
        if ~ opt.flag_use_peak
            plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color','black','linewidth',2)
        else
            plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color','black','linewidth',2)
        end
        
        
        % modify x ticks
        xaxisvalues = [0:30:opt.sta_pre_frames+opt.sta_post_frames];
        xticks(xaxisvalues)
        xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
        end
        axis square
        title([ensemble_names{e},' Var. ' num2str(this_photo_type)])
        plot_count = plot_count+1;
    end
    
end
export_fig([fig_save_path filesep 'OppoTargetEnsemble_STATrace_' strrep(caiman_file,'.mat','.png')])


%% Plot trial average as imagesc - photo response
% plot_photo_types = {'_photostim','_nonphotostim','_dummyphotostim'};
plot_photo_types = {'_photostim','_dummyphotostim'};

plot_condition_types = [1,2];
cnn_thresh = min(min(cell2mat({cell_struct(trigger_idx).cnn_score})),0.1);
plot_cell_idx = find(cell2mat({cell_struct(:).cnn_score})>=cnn_thresh);

plot_trial_avg_imagesc(cell_struct,trial_indices,photo_ensembles,plot_cell_idx,plot_photo_types,plot_condition_types,opt)
export_fig([fig_save_path filesep 'TrialAvgStim_STATrace_' strrep(caiman_file,'.mat','.png')])

%% plot trial avg compare outcome
plot_photo_types = {'_photostim_correct','_dummyphotostim_correct','_photostim_incorrect','_dummyphotostim_incorrect'};
plot_condition_types = [1,2];
cnn_thresh = min(min(cell2mat({cell_struct(trigger_idx).cnn_score})),0.1);
plot_cell_idx = find(cell2mat({cell_struct(:).cnn_score})>=cnn_thresh);
plot_trial_avg_imagesc(cell_struct,trial_indices,photo_ensembles,plot_cell_idx,plot_photo_types,plot_condition_types,opt)


%% online trajectory plots
% set field names and colors
test_opt = opt;
% fds_of_interest = {'stim_1_var_1_correct','stim_1_var_1_incorrect','stim_1_var_1_photostim','stim_1_var_1_nonphotostim','stim_1_var_1_dummyphotostim','stim_5_photostim_1',......
%     'stim_2_var_2_correct', 'stim_2_var_2_incorrect','stim_2_var_2_photostim','stim_2_var_2_nonphotostim','stim_2_var_2_dummyphotostim','stim_5_photostim_2'};
fds_of_interest = fields(trial_indices);
test_opt.fd_names = fds_of_interest;
% add colors 
for i = 1:numel(test_opt.fd_names )
    this_fd = test_opt.fd_names {i};
    
    % texture trials
    if contains(this_fd,'stim_1') && contains(this_fd,'_correct')
        trial_color.(this_fd) = trial_color.correct_stim1;
    end
    if contains(this_fd,'stim_1')&& contains(this_fd,'_incorrect')
        trial_color.(this_fd) = trial_color.incorrect_stim1;
    end
    
    if contains(this_fd,'stim_2') && contains(this_fd,'_correct')
        trial_color.(this_fd) = trial_color.correct_stim2;
    end
    if contains(this_fd,'stim_2')&& contains(this_fd,'_incorrect')
        trial_color.(this_fd) = trial_color.incorrect_stim2;
    end
    
    % add photo types
    if contains(this_fd,'_photostim')
        trial_color.(this_fd) = trial_color.photo;
    end
    if contains(this_fd,'_nonphotostim')
        trial_color.(this_fd) = trial_color.nonphoto;
    end
    if contains(this_fd,'_dummyphotostim')
        trial_color.(this_fd) = trial_color.dummyphoto;
    end
    if contains(this_fd,'_oppophotostim')
        trial_color.(this_fd) = trial_color.oppophoto;
    end
    
    % add shade to incorrect trials
    if contains(this_fd,'_photostim_incorrect')
        trial_color.(this_fd) = shade(trial_color.photo,.5);
    end
    if contains(this_fd,'_nonphotostim_incorrect')
        trial_color.(this_fd) = shade(trial_color.nonphoto,.5);
    end
    if contains(this_fd,'_dummyphotostim_incorrect')
        trial_color.(this_fd) = shade(trial_color.dummyphoto,.5);
    end
    
    % miss trials
    if contains(this_fd,'_miss')
        trial_color.(this_fd) = [94, 34, 92]./255;
    end
    
    % catch trials
    if contains(this_fd,'stim_5')
        trial_color.(this_fd) = [.5 .5 .5];
    end
    
end
test_opt.trial_color = trial_color;
disp('added colors')
%% choose trial types to plot
plot_fds = {'stim_1_var_1_correct','stim_1_var_1_incorrect','stim_1_var_1_photostim','stim_1_var_1_nonphotostim','stim_1_var_1_dummyphotostim','stim_5_photostim_1',......
    'stim_2_var_2_correct', 'stim_2_var_2_incorrect','stim_2_var_2_photostim','stim_2_var_2_nonphotostim','stim_2_var_2_dummyphotostim','stim_5_photostim_2'};
plot_num_cols = 3;
IF_PLOT_AVG_ONLY = 1;
IF_NORMALISE = 1;
% online recorded trajectory
xlimit = [30,270];
ylimit = [-1 1];

% ylimit = [-500 500];

plot_fds1 = {'stim_1_var_1_photostim','stim_1_var_1_nonphotostim','stim_1_var_1_dummyphotostim','stim_1_var_1_oppophotostim'};
plot_pop_vectors(traj_struct,plot_fds1,1,test_opt,...
    'noise_thresh',thresh_sd,'plot_ylabel','Projection',...
    'ylimit',ylimit,'xlimit',xlimit,...
    'plot_num_cols',plot_num_cols,'IF_PLOT_RAW_ONLY',0,'IF_NORMALISE',IF_NORMALISE,'IF_PLOT_AVG_ONLY',IF_PLOT_AVG_ONLY)
suptitle('Closed-loop condition trials, online trajectory, tex1')
export_fig([fig_save_path filesep 'OnlineTrajTex1' strrep(caiman_file,'.mat','.png')])

plot_fds2 = {'stim_2_var_2_photostim','stim_2_var_2_nonphotostim','stim_2_var_2_dummyphotostim','stim_2_var_2_oppophotostim'};
plot_pop_vectors(traj_struct,plot_fds2,1,test_opt,...
    'noise_thresh',thresh_sd,'plot_ylabel','Projection',...
    'ylimit',ylimit,'xlimit',xlimit,...
    'plot_num_cols',plot_num_cols,'IF_PLOT_RAW_ONLY',0,'IF_NORMALISE',IF_NORMALISE,'IF_PLOT_AVG_ONLY',IF_PLOT_AVG_ONLY)
suptitle('Closed-loop condition trials, online trajectory, tex2')
export_fig([fig_save_path filesep 'OnlineTrajTex2' strrep(caiman_file,'.mat','.png')])


% control_fds_of_interest = {'stim_5_photostim_1','stim_5_photostim_2','stim_5_nonphotostim'};
% plot_pop_vectors(traj_struct,control_fds_of_interest,1,test_opt,...
%     'noise_thresh',thresh_sd,'plot_ylabel','Projection',...
%     'ylimit',ylimit,...
%     'plot_num_cols',3,'IF_PLOT_RAW_ONLY',0)
% suptitle('Catch condition trials, online trajectory')


%% trajectory computed using raw_sta_traces
[decod_proj_struct] = get_projections(raw_cell_struct(trigger_idx),norm_weights,plot_fds,'bias',-norm_thresh,'IS_CELL_STRUCT',1);
plot_pop_vectors(decod_proj_struct,plot_fds,1,opt,...
    'ylimit',ylimit,'plot_ylabel','Projection','plot_num_cols',plot_num_cols,'IF_PLOT_RAW_ONLY',1,'IF_NORMALISE',IF_NORMALISE)
suptitle('Decoder projections')

%% ============== POST-HOC DECODER ========================
% need to run PosthocProcTexture before this
%% load posthoc decoder
[decod_file,decod_path] = uigetfile('*.mat','Select PosthocDecoders file',caiman_path);
decoders_struct =  load(fullfile(decod_path,decod_file));
decoders_struct = decoders_struct.save_decod_struct;
disp(['Loaded file :',fullfile(decod_path,decod_file)])
%% project to stim decoders
% get projections on stim decoder
% fds_of_interest = {'stim_2_var_2_photostim','test_trial_idx'}
stim_norm_weights = decoders_struct.stim.norm_weights;
stim_norm_thresh = decoders_struct.stim.norm_thresh;
avg_frames = opt.sta_avg_frames;
test_opt.save_name_ext = 'posthoc';
test_opt.save_path = fig_save_path;
% targets_idx_1 = decod_struct.target_ensembles{1};
% [~,idx1] = intersect(trigger_idx,targets_idx_1);
% target_weights_1 = stim_norm_weights(idx1)
% targets_idx_2 = decod_struct.target_ensembles{2};
% [~,idx2] = intersect(trigger_idx,targets_idx_2);
% target_weights_2 = stim_norm_weights(idx2)

ylimit = [-1000 1000];
stim_proj_struct = struct();
[stim_proj_struct] = get_projections(cell_struct(trigger_idx),stim_norm_weights,plot_fds,'proj_struct',stim_proj_struct,'bias',-stim_norm_thresh,'IS_CELL_STRUCT',1);
this_sup_title = ['StimDecoderProj_' test_opt.save_name_ext strrep(caiman_file,'.mat','')];
plot_pop_vectors(stim_proj_struct,plot_fds,1,test_opt,...
        'plot_ylabel','Hit vs CR decoder Projection','sup_title',this_sup_title,...
        'plot_num_cols',plot_num_cols,'IF_PLOT_RAW_ONLY',0,'ylimit',ylimit,'IF_SAVE_PLOT',true)

% compare projection on stim axis trials
figure('name','stim decoder projection','position',[100 100 800 800])
values = structfun(@(x)mean(x(:,avg_frames),2),stim_proj_struct,'UniformOutput',false);
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(values),'UniformOutput',false));
scatter_cmp_conditions(values,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1,'plot_stats',1,'test_type','ttest');
suptitle('Stim decoder projections')
xtickangle(45)
export_fig([fig_save_path filesep 'StimDecoderProj_' test_opt.save_name_ext strrep(caiman_file,'.mat','.png')])


%% get projections on choice decoder
choice_norm_weights = decoders_struct.choice.norm_weights;
choice_norm_thresh = decoders_struct.choice.norm_thresh;
ylimit = [-60 60];

choice_proj_struct = struct();
[choice_proj_struct] = get_projections(cell_struct(trigger_idx),choice_norm_weights,...
    plot_fds,'proj_struct',choice_proj_struct,'bias',-choice_norm_thresh,'IS_CELL_STRUCT',1);
% compare projection on choice axis trials
this_sup_title = ['ChioceDecoderProj_' test_opt.save_name_ext strrep(caiman_file,'.mat','')];
plot_pop_vectors(choice_proj_struct,plot_fds,1,test_opt,...
    'plot_ylabel','Go vs No-GO Projection','sup_title',this_sup_title,...
    'plot_num_cols',plot_num_cols,'IF_PLOT_RAW_ONLY',0,'ylimit',ylimit,'IF_SAVE_PLOT',true)
% export_fig([fig_save_path filesep 'RawChoiceDecoderProj_' test_opt.save_name_ext strrep(caiman_file,'.mat','.png')])

figure('name','choice decoder projection','position',[100 100 800 800])
values = structfun(@(x)mean(x(:,avg_frames),2),choice_proj_struct,'UniformOutput',false);
fd_colors =  cell2mat(cellfun(@(f)getfield(trial_color,f),fields(values),'UniformOutput',false));
scatter_cmp_conditions(values,[],...
    1,fd_colors,'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1,'plot_stats',1,'test_type','ttest');
suptitle('Choice decoder projections')
xtickangle(45)

export_fig([fig_save_path filesep 'ChoiceDecoderProj_' test_opt.save_name_ext strrep(caiman_file,'.mat','.png')])


%% ============================     END    ================================

%% debug
check_fd = 'stim_2_var_2_dummyphotostim';
figure;imagesc(stim_proj_struct.(check_fd)); title('stim proj')
figure;imagesc(choice_proj_struct.(check_fd)); title('choice proj')
trial_indices.(check_fd)
%%
check_trial_idx = 205;
plot_frame_range = double([sens_stim_frames(check_trial_idx) sens_stim_frames(check_trial_idx+1)]);
figure('name','online trigger cell traces'); hold on
plot_offset = 5;
cell_count = 1;

this_idx = trigger_idx;
xlim(plot_frame_range)

for i = 1:length(this_idx)
    ii = this_idx(i);
%     plot(cell_struct(ii).filtC+cell_count*plot_offset,'black');
    plot(cnm_struct(cell_struct(ii).cnm_idx).deconvC_full+cell_count*plot_offset,'black');

    text(plot_frame_range(1),cell_count*plot_offset,['Cell ' num2str(ii) ', W ' num2str(stim_norm_weights(cell_count))], 'horizontalalignment','right', 'color','black')
     cell_count = cell_count+1;

end
set(gca,'ytick',[])
ylim([0 plot_offset*(cell_count+1)])
for i = 1:numel(sens_stim_frames)
    this_color = trial_color.(['stim' num2str(trials.stim_type(i) )]);
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.stimon_frame,ylim,'color',this_color) % trial-on + roughly first touch
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.gocue_frame,ylim,'color',this_color,'linestyle',':') % go-cue
end

for i = 1:numel(photo_stim_frames)
    plot([photo_stim_frames(i) photo_stim_frames(i)],ylim,'color',trial_color.photostim)
end

for i = 1:numel(dummy_photo_stim_frames)
    plot([dummy_photo_stim_frames(i) dummy_photo_stim_frames(i)],ylim,'color',trial_color.photostim,'linestyle',':')
end
hold on;plot(online_traj,'color',[.5 .5 .5])

% plot online trajectory

figure('name','online traj trace'); hold on
plot(online_traj,'color',[.5 .5 .5])
xlim(plot_frame_range)

for i = 1:numel(sens_stim_frames)
    this_color = trial_color.(['stim' num2str(trials.stim_type(i) )]);
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.stimon_frame,ylim,'color',this_color) % stim-on
    plot([sens_stim_frames(i) sens_stim_frames(i)]+opt.gocue_frame,ylim,'color',this_color,'linestyle',':') % go-cue
end

for i = 1:numel(photo_stim_frames)
    plot([photo_stim_frames(i) photo_stim_frames(i)],ylim,'color',trial_color.photostim,'linestyle',':') % stim-on
end

plot(xlim,[0 0],'color','black')


