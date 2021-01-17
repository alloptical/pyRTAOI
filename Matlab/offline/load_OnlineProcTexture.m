%% 2019 ZZ
% TO try:
% 1. identify anchor points and fit spline in PCs

%% load/save paths
data_folder = 'D:\TextureData\pyrtaoi_results\';
save_folder = 'D:\TextureData\pyrtaoi_proc_data\';
trial_full_path = 'W:\forZoe\cb183\20190429\1_2_3_4_5_6_7\trials';
data_name = '20190429_CB183';
save_path = [save_folder filesep data_name];
if ~exist(save_path, 'dir')
    mkdir(save_path)
end

fig_save_path = [save_path filesep 'figures' ]; % figures are saved here
if ~exist(fig_save_path, 'dir')
    mkdir(fig_save_path)
end

output_file_path = [save_path filesep 'output_files'];
if ~exist(output_file_path, 'dir')
    mkdir(output_file_path)
end

%% init opt and color structure
opt = tex_init_opt(); % get default parameters
opt.idx_fields = {'tex','port','outcome'};
opt.save_path = fig_save_path;
opt.output_path = output_file_path;
opt.exp_name =  strrep(data_name,'_','-');
opt.input_data_path = [data_folder filesep data_name];
opt.pc_trials = [0 60];
opt.frame_range = 1:210;
opt.trial_frames = length(opt.frame_range);
opt.fd_names_fa = {'st_correct_stim_1','st_incorrect_stim_1','st_correct_stim_2','st_incorrect_stim_2'}; % using raw trial traces for PCA
opt.all_fd_names_fa = {'st_correct_stim_1','st_incorrect_stim_1','st_correct_stim_2','st_incorrect_stim_2','st_miss_stim_1','st_miss_stim_2'}; % using raw trial traces for PCA
opt.ref_fields =  {'st_correct_stim_1','st_correct_stim_2'};
opt.pre_trial_bs_frames = 20:40; %
opt.Nstd = 2; % used for compare classification accuracy with random shuffle
opt.trialoutcome_peri_trials = 3; % number of trials pre and post current trial used for local trialoutcome analysis
color = tex_init_color(opt);
figure('name','color','position',[400 400 400 400])
hold on
for f = 1:numel(opt.all_fd_names_fa)
    plot([0,.5],[1,1].*f,'linewidth',2,'color',color.(opt.all_fd_names_fa{f}))
end
box off
axis off
legend('Correct stim1','Incorrect stim1','Correct stim2','Incorrect stim2','Miss stim1','Miss stim2','location','EastOutside')

%% load caiman data
[cell_struct,fov_struct] = concat_OnlineProc(opt.input_data_path);
opt.fov_size = size(fov_struct.com_fov);
opt.ds_factor = 512/opt.fov_size(1);
%% load trial struct (from CB)
load(trial_full_path)

%% check motivation
figure('name','performance')
hold on
plot(cumsum(trials.correct),'color',tint([0,1,0],.5),'linewidth',2)
plot(cumsum(trials.incorrect),'color',tint([1,0,0],.5),'linewidth',2)
plot(cumsum(trials.miss),'color',tint([0,0,0],.5),'linewidth',2)
xlabel('Trial Number')
ylabel('Cumulated sum')

% find motivation drop point -
% - use start of the chunck when miss trial keeps increasing
% trial_window = 50;
% pc_miss = movmean((movsum(trials.miss,[trial_window 0])/trial_window),10);
% pc_miss(1:trial_window) = 0; % ignore first trials
% last_perform_trial = find_continue_increase(pc_miss,30)+1;

% - alternatively use 'trial bias' - ask christina how it is calculated
last_perform_trial = find_continue(trials.bias,30,'const');
last_perform_trial = last_perform_trial(end);
plot([last_perform_trial,last_perform_trial],ylim,':','color','black')
legend('correct', 'incorrect','miss','miss increase','location','northwest')


%% make cell structure with sta
train_trial_range = [1, round(last_perform_trial*opt.pc_trials(end)/100)];
test_trial_range = [round(last_perform_trial*opt.pc_trials(end)/100)+1, last_perform_trial];

opt.full_trace_fd = 'devonC_full';
[train_cell_struct,train_num_trials,train_trial_indices,cell_max] = make_texture_sta_struct( cell_struct,trials,train_trial_range,opt, 'IF_TRAINING',1,'IF_NORMALISE',true);
train_lick_struct = make_texture_lick_struct( trials,train_trial_range,opt,'IF_TRAINING',1);

[test_cell_struct,test_num_trials,test_trial_indices] = make_texture_sta_struct( cell_struct,trials,test_trial_range,opt, 'IF_TRAINING',0,'cell_max',cell_max,'IF_NORMALISE',true);
test_lick_struct = make_texture_lick_struct( trials,test_trial_range,opt, 'IF_TRAINING',0);

%% show number of trials
figure('name','number of trials')
this_plot_color = [];
this_fields = fields(train_num_trials);
for  i = 1:numel(this_fields)
    this_plot_color = [this_plot_color;color.(['st_' this_fields{i}])];
end
subplot(2,1,1)
scatter_cmp_conditions(train_num_trials,'Num. trials',1,this_plot_color,'connect_scatter',0,'VeryBriefXlabel',1);
ylabel('Num trials (train)')

subplot(2,1,2)
scatter_cmp_conditions(test_num_trials,'Num. trials',1,this_plot_color,'connect_scatter',0,'VeryBriefXlabel',1);
ylabel('Num trials (test)')

%% get sta amplitude
opt.sta_amp_frames = opt.withold_frames_adj;
[  train_cell_struct ] = get_texture_sta_amp( train_cell_struct,opt.all_fd_names_fa,opt);
[  test_cell_struct ] = get_texture_sta_amp( test_cell_struct,opt.all_fd_names_fa,opt);

%% get baseline correlation
corr_frame_range  = 1:60;
[  corr_struct ] = get_texture_corr( train_cell_struct,cell_idx_struct,corr_frame_range,opt.all_fd_names_fa);
[  test_corr_struct ] = get_texture_corr( test_cell_struct,cell_idx_struct,corr_frame_range,opt.all_fd_names_fa);

%% trial outcome before incorrect trials
trialoutcome_struct = make_texture_trialoutcome_struct( trials,train_trial_indices,opt);
test_trialoutcome_struct = make_texture_trialoutcome_struct( trials,test_trial_indices,opt);
%% get cell identity - stim and port neurons
[cell_idx_struct,outcomeAUC,stimulusAUC,texAUC,pcorrect,pincorrect] = make_cell_idx_struct(train_cell_struct,opt);
% [test_cell_idx_struct,test_outcomeAUC,test_stimulusAUC,test_texAUC,test_pcorrect,test_pincorrect] = make_cell_idx_struct(test_cell_struct,opt);

%% POPULATION ANALYSIS
%% PCA
traj_struct = struct(); %
opt.num_final_frames = 5;
[ traj_struct ] = get_pca_trajs_trials( traj_struct,train_cell_struct,cell_idx_struct,opt,...
    'IF_USE_TRIAL_AVG',false,'IF_USE_WITHOLD_FRAMES',true,'IF_USE_CORRECT_TRIALS',false);
[ traj_struct ] = add_pca_trajs_trials( traj_struct,train_cell_struct,cell_idx_struct,opt,'trace_fd_names',{'st_miss_stim_1','st_miss_stim_2'} );
%% LDA
opt.m_fa_ld = 4;
[ traj_struct ] = get_ld_withold(traj_struct,opt,'traj_field','all_pca_struct' ); % this uses k-fold cross-validation error, withould subsampling

% opt.ld_cmp_fds = {'st_correct_stim_1','st_correct_stim_2'}; %
% [ traj_struct ] = get_ld_framewise(traj_struct,opt,'traj_field','all_pca_struct' );

%% Selectivity
opt.dist_components = 1:4; % using first 3 PCs for selectivity index
[ traj_struct, traj_dist_exc] = get_selectivity_index( traj_struct,opt,'traj_field',traj_field_struct,...
    'IF_USE_WITHOLD_LD',true,'IF_PLOT_TRAJ',false,'color',color,'cmp_fields',opt.all_fd_names_fa); % project to LD before calculating selectivity

[ traj_struct ] =  get_sel_classfier( traj_struct, opt,'IF_DO_LDA',false,'IF_DO_LOG',true);

for idx = 1:numel(opt.idx_fields)
    this_idx_field = opt.idx_fields{idx};
    traj_struct.sel_classif_accuracy.(this_idx_field)  =  get_classfier_dc_time( traj_struct.selectivity_idx.(this_idx_field), [],opt,'threshold',traj_struct.sel_thresh_log.(this_idx_field) );
end

%% Sparse linear regression (see Scholz et al(Leifer) biorxiv 2018)
% find cells that are useful in predicting trial outcome
lasso_idx_field = 'outcome'; % only looking at subsets of tex-relavent cells (i.e. discriminate texture also have higher activity during withold)
lasso_struct = [];
lasso_input = make_lasso_input(train_cell_struct(cell_idx_struct.(lasso_idx_field)),{'st_correct_stim_1','st_incorrect_stim_1','st_correct_stim_2','st_incorrect_stim_2'},opt);
lasso_struct  =  get_cd_classfier( lasso_struct,lasso_input, opt,'IF_FRAMEWISE',false,'MODEL_TYPE','elastic');
lasso_struct.thresh_fix.stim1 = -lasso_struct.coef0.stim1;
lasso_struct.thresh_fix.stim2 = -lasso_struct.coef0.stim2;
lasso_struct.thresh_framewise.stim1 = repmat(lasso_struct.thresh_fix.stim1,[1,opt.trial_frames]);
lasso_struct.thresh_framewise.stim2 = repmat(lasso_struct.thresh_fix.stim2,[1,opt.trial_frames]);
% add cells with non-zero coeffs to cell_idx_struct
cell_idx_struct.lasso_tex1 = cell_idx_struct.(lasso_idx_field)(lasso_struct.coef.stim1~=0);
cell_idx_struct.lasso_tex2 = cell_idx_struct.(lasso_idx_field)(lasso_struct.coef.stim2~=0);
cell_idx_struct.lasso_tex = cell_idx_struct.(lasso_idx_field)(lasso_struct.coef.stim1~=0|lasso_struct.coef.stim2~=0);

% multipy regression coeff
lasso_proj_struct = struct();
for stim_type = 1:2
    this_tex = num2str(stim_type);
    this_coef = lasso_struct.coef.(['stim' this_tex])(lasso_struct.coef.(['stim' this_tex])~=0);
    this_coef0 = lasso_struct.coef0.(['stim' this_tex]);
    lasso_proj_struct = get_coding_direction(train_cell_struct(cell_idx_struct.(['lasso_tex' this_tex])),color,opt,...
        'IF_GET_PROJ_ONLY', true,'cd',this_coef,'this_proj_struct',lasso_proj_struct,'trace_fds',{['st_correct_stim_' this_tex],['st_incorrect_stim_' this_tex]});
end

lasso_struct =  get_classfier_dc_time( lasso_proj_struct, lasso_struct,opt,'IF_FRAMEWISE',false,'IF_REVERSE',false);
lasso_struct = detect_error_by_cd(lasso_struct,lasso_proj_struct,train_lick_struct,opt,lasso_struct,'IF_FRAMEWISE',false,'IF_REVERSE',false);

%% Coding direction
opt.idx_fields =   {'lasso_tex','port','outcome'};
IF_FRAMEWISE = true;
train_cd_struct = struct();
for idx = 1:numel(opt.idx_fields)
    this_idx_field = opt.idx_fields{idx};
    [this_proj_struct,this_cd,this_db,this_proj_struct_framewise,this_framewise_cd] = get_coding_direction(train_cell_struct(cell_idx_struct.(this_idx_field)),color,opt,...
        'ref_fds',{'sta_correct_stim_1','sta_correct_stim_2'},'correct_fds',{'st_correct_stim_1','st_correct_stim_2'});
    cd_struct = [];
    [ cd_struct ] =  get_cd_classfier( cd_struct,this_proj_struct_framewise, opt,'IF_FRAMEWISE',IF_FRAMEWISE );
    [ cd_struct ] =  get_classfier_dc_time( this_proj_struct_framewise, cd_struct,opt,'IF_FRAMEWISE',IF_FRAMEWISE );
    cd_struct = detect_error_by_cd(cd_struct,this_proj_struct_framewise,train_lick_struct,opt,cd_struct,'IF_FRAMEWISE',IF_FRAMEWISE);
    cd_struct.this_cd = this_cd;
    cd_struct.this_framewise_cd = this_framewise_cd;
    cd_struct.this_proj_struct_framewise = this_proj_struct_framewise;
    cd_struct.this_proj_struct = this_proj_struct;
    train_cd_struct.(this_idx_field) = cd_struct;
end
%% Generalised linear model
% predicting trial outcome using pre-trial trialoutcome, baseline correlation and average activity level in withold window fit model
% (see make_glm_input)
stim_type = 1;
[glminput,glmref,norm_glminput,norm_glmref] = make_glm_input(train_cell_struct,train_num_trials,trialoutcome_struct,corr_struct,cell_idx_struct,stim_type,opt);
selected_inputs = 4:5;
mdl =  fitglm(glminput(:,selected_inputs),glmref,'Distribution','binomial','link','logit');
% test model
[test_glminput,test_glmref,test_norm_glminput,test_norm_glmref] = make_glm_input(test_cell_struct,test_num_trials,test_trialoutcome_struct,test_corr_struct,cell_idx_struct,stim_type,opt);
ytest = predict(mdl,test_glminput(:,selected_inputs));

%% PROCESS test data
% apply decoders test trials
opt_test = opt;
opt_test.pc_trials = [60 100];
test_traj_struct = [];
[ test_traj_struct ] = get_pca_trajs_trials( test_traj_struct,test_cell_struct,cell_idx_struct,opt_test,...
    'train_struct',traj_struct,'pca_struct_name',traj_field_struct);
[ test_traj_struct, test_traj_dist_exc] = get_selectivity_index( test_traj_struct,opt,'traj_field',traj_field_struct,...
    'IF_USE_WITHOLD_LD',true,'IF_PLOT_TRAJ',false,'color',color,'train_struct',traj_struct);

for idx = 1:numel(opt.idx_fields)
    this_idx_field = opt.idx_fields{idx};
    test_traj_struct.sel_classif_accuracy.(this_idx_field)  =  get_classfier_dc_time( test_traj_struct.selectivity_idx.(this_idx_field), [],opt,'threshold',traj_struct.sel_thresh_log.(this_idx_field) );
end
%% test lasso regression
test_lasso_input = make_lasso_input(test_cell_struct(cell_idx_struct.(lasso_idx_field)),{'st_correct_stim_1','st_incorrect_stim_1','st_correct_stim_2','st_incorrect_stim_2'},opt);
test_lasso_proj_struct = struct();
for stim_type = 1:2
    this_tex = num2str(stim_type);
    this_coef = lasso_struct.coef.(['stim' this_tex])(lasso_struct.coef.(['stim' this_tex])~=0);
    this_coef0 = lasso_struct.coef0.(['stim' this_tex]);
    test_lasso_proj_struct = get_coding_direction(test_cell_struct(cell_idx_struct.(['lasso_tex' this_tex])),color,opt,...
        'IF_GET_PROJ_ONLY', true,'cd',this_coef,'this_proj_struct',test_lasso_proj_struct,'trace_fds',{['st_correct_stim_' this_tex],['st_incorrect_stim_' this_tex]});
end
[ lasso_struct_test] =  get_classfier_dc_time( test_lasso_proj_struct, lasso_struct,opt,'IF_FRAMEWISE',false);
lasso_struct_test = detect_error_by_cd(lasso_struct_test,test_lasso_proj_struct,test_lick_struct,opt,lasso_struct,'IF_FRAMEWISE',false);


%% test coding direction
test_cd_struct = struct();
for idx = 1:numel(opt.idx_fields)
    this_idx_field = opt.idx_fields{idx};
    cd_struct =  train_cd_struct.(this_idx_field);
    [this_proj_struct_test,~,~,this_proj_struct_framewise_test,~] = get_coding_direction(test_cell_struct(cell_idx_struct.(this_idx_field)),color,opt,...
        'ref_fds',{'sta_correct_stim_1','sta_correct_stim_2'},'correct_fds',{'st_correct_stim_1','st_correct_stim_2'},...
        'cd',cd_struct.this_cd,'framewise_cd',cd_struct.this_framewise_cd);
    [ cd_struct_test] =  get_classfier_dc_time( this_proj_struct_framewise_test, cd_struct,opt,'IF_FRAMEWISE',IF_FRAMEWISE);
    cd_struct_test = detect_error_by_cd(cd_struct_test,this_proj_struct_framewise_test,test_lick_struct,opt,cd_struct,'IF_FRAMEWISE',IF_FRAMEWISE);
    cd_struct_test.this_proj_struct_framewise = this_proj_struct_framewise_test;
    test_cd_struct.(this_idx_field) = cd_struct_test;
end

%% MAKE OUTPUT FILES
opt.target_idx_fd = 'tex1'; % cells to stimulate
opt.trigger_idx_fd = 'outcome'; % cells for readout
pop_params.weights = []; % thess should be determined by population analysis results
pop_params.thresh = [];
pop_params.frames_enable_trigger = [];

[output] = generate_cell_idx_file(cell_struct,cell_idx_struct,[],opt);


%% ============================= PLOTS =====================================
IF_SAVE_PLOTS = false;
%% POPULATION ANALYSIS PLOTS
%% plot trajectories
traj_field_struct = 'all_pca_struct';
plot_check_pcs(traj_struct,opt,color,'IF_SAVE', IF_SAVE_PLOTS,'save_name_ext','train')
plot_fa_trajs_avg( traj_struct,color,opt,'traj_field',traj_field_struct,'IF_SAVE', IF_SAVE_PLOTS,'last_frame',150 )
plot_fa_trajs_raw( traj_struct,color,opt,'traj_field',traj_field_struct,'IF_SAVE', IF_SAVE_PLOTS,'last_frame',150 )

plot_check_pcs(test_traj_struct,opt,color,'IF_SAVE', IF_SAVE_PLOTS,'save_name_ext','test')
plot_fa_trajs_avg( test_traj_struct,color,opt,'traj_field',traj_field_struct,'IF_SAVE', IF_SAVE_PLOTS,'last_frame',150,'save_name_ext','test' )

%% plot selectivity
results = struct(); % dummy
[results,traj_struct ]= plot_selectivity_index( traj_struct,traj_dist_exc,color,opt,results,'IF_SAVE',IF_SAVE_PLOTS);
plot_selectivity_index_raw( traj_struct,color,opt,'classf_type','log','save_name_ext','Log train','IF_SAVE',IF_SAVE_PLOTS);
plot_dist_to_correct( traj_dist_exc,train_lick_struct,traj_struct,color,opt,'save_name_ext',' train','IF_SAVE',IF_SAVE_PLOTS);
traj_struct = plot_shuf_accuracy(traj_struct,opt,'save_name_ext',' train', 'IF_NEW_STRUCT',true);

[results,test_traj_struct ]= plot_selectivity_index( test_traj_struct,test_traj_dist_exc,color,opt,results,'IF_SAVE',IF_SAVE_PLOTS);
% plot_selectivity_index_raw( test_traj_struct,color,opt,'classf_type','log','save_name_ext','Log test','IF_SAVE',IF_SAVE_PLOTS);
plot_dist_to_correct( test_traj_dist_exc,test_lick_struct,traj_struct,color,opt,'save_name_ext',' test','IF_SAVE',IF_SAVE_PLOTS);
test_traj_struct = plot_shuf_accuracy(test_traj_struct,opt,'save_name_ext',' test', 'IF_NEW_STRUCT',true);
plot_hit_fa_rate(test_traj_struct.sel_classif_accuracy,opt,'save_name_ext',' test', 'IF_NEW_STRUCT',true,'IF_SAVE',IF_SAVE_PLOTS);

%% plot coding direction
opt.disc_frames_fd = 'shuf_disc_frames'; % when classifier performs better than mean+3std chance performance
for idx = 1:numel(opt.idx_fields)
    this_idx_field = opt.idx_fields{idx};
    if strcmp(this_idx_field,'stim')
        continue
    end
    cd_struct =  train_cd_struct.(this_idx_field);
    cd_struct_test = test_cd_struct.(this_idx_field);
    
    figure('name','train coding direction plot','units','normalized','outerposition',[0 0 1 1]);
    plot_coding_direction(cd_struct,cd_struct.this_proj_struct_framewise,train_lick_struct,color,opt)
    this_title = ['Train ' opt.exp_name ' framewise cd error detection ' this_idx_field];
    suptitle(strrep(this_title,'_',' '))
    
    figure('name','test coding direction plot','units','normalized','outerposition',[0 0 1 1]);
    plot_coding_direction(cd_struct_test,cd_struct_test.this_proj_struct_framewise,test_lick_struct,color,opt)
    this_title = ['Test ' opt.exp_name ' framewise cd error detection ' this_idx_field];
    suptitle(strrep(this_title,'_',' '))
end

%% plot lasso regression
figure('name','test lasso regression plot','units','normalized','outerposition',[0 0 1 1]);
plot_coding_direction(lasso_struct_test,test_lasso_proj_struct,test_lick_struct,color,opt)

figure('name','train lasso regression plot','units','normalized','outerposition',[0 0 1 1]);
plot_coding_direction(lasso_struct,lasso_proj_struct,train_lick_struct,color,opt)

%% compare correlation (baseline)
this_plot_color = [];
plot_idx_fd = 'tex1';
fd_names = fields(corr_struct.(plot_idx_fd));
for  i = 1:numel(fd_names)
    this_fd_name = strrep(fd_names{i},'_abs_avg','');
    this_plot_color = [this_plot_color;color.(this_fd_name)];
end
figure
scatter_cmp_conditions(corr_struct.(plot_idx_fd),'correlation',1,this_plot_color,'connect_scatter',0,'VeryBriefXlabel',0,'plot_stats',1);
ylabel('Baseline correlation (abs)')
%% ========================== GLM PLOTS ===================================
% check individual predictor significance
figure('name','train individual predictor stats')
for i = 1:num_predictors
subplot(1,num_predictors,i)
values = struct();
values.correct = glminput(glmref==0,i);
values.incorrect = glminput(glmref==1,i);
scatter_cmp_conditions(values,' ',1,[],'connect_scatter',0,'VeryBriefXlabel',0,'plot_stats',1,'add_jitter',1);
ylabel(['tex' num2str(stim_type),' ' strrep(predictor_names{i},'_',' ')])
end

% check individual predictor significance
figure('name','test individual predictor stats')
for i = 1:num_predictors
subplot(1,num_predictors,i)
values = struct();
values.correct = test_glminput(test_glmref==0,i);
values.incorrect = test_glminput(test_glmref==1,i);
scatter_cmp_conditions(values,' ',1,[],'connect_scatter',0,'VeryBriefXlabel',0,'plot_stats',1,'add_jitter',1);
ylabel(['tex' num2str(stim_type),' ' strrep(predictor_names{i},'_',' ')])
end
%%
figure('name','glm output')
subplot(1,2,1)
ytrain = predict(mdl,glminput(:,selected_inputs));
values = struct();
values.correct = ytrain(glmref==0);
values.incorrect = ytrain(glmref==1);
scatter_cmp_conditions(values,' ',1,[],'connect_scatter',0,'VeryBriefXlabel',0,'plot_stats',1,'add_jitter',1);
ylabel('GLM output, trian')

subplot(1,2,2)
values = struct();
values.correct = ytest(test_glmref==0);
values.incorrect = ytest(test_glmref==1);
scatter_cmp_conditions(values,' ',1,[],'connect_scatter',0,'VeryBriefXlabel',0,'plot_stats',1,'add_jitter',1);
ylabel('GLM output, test')

%% ==================== BEHAVIORAL PERFORMANCE PLOTS ======================
figure('name','cumsum performance')
hold on
plot(cumsum(trials.correct),'color',tint([0,1,0],.5),'linewidth',2)
plot(cumsum(trials.incorrect),'color',tint([1,0,0],.5),'linewidth',2)
plot(cumsum(trials.miss),'color',tint([0,0,0],.5),'linewidth',2)
plot(cumsum(trials.fa),'color',[0,0,0],'linewidth',2)

xlabel('Trial Number')
ylabel('Cumulated sum')
legend('correct', 'incorrect','miss','FA','location','northwest')
%%
figure('name','behavior full traces');
subplot(2,1,1)
hold on
trial_outcomes = {'miss','correct','incorrect','fa'};
for  i = 1:numel(trial_outcomes)
    trial_start_vector = zeros(1,tot_frames);
    trial_start_vector(trials.trial_start_frames(trials.(trial_outcomes{i})==1)) = 1;
    stem(trial_start_vector.*(num_init_comp*plot_offset),'MarkerEdgeColor','none','MarkerFaceColor','none','color',color.([trial_outcomes{i},'_trial']),'linewidth',1)
end

% mark training and text trials
plot(trials.trial_start_frames(train_trial_range),[0 0],'linewidth',2,'color','black')
plot(trials.trial_start_frames(test_trial_range),[0 0],'linewidth',2,'color',[.5 .5 .5])

% mark stim type
trial_start_vector = zeros(1,tot_frames);
trial_start_vector(trials.trial_start_frames(trials.stim_type==1))=-200;
stem(trial_start_vector,'MarkerEdgeColor','none','MarkerFaceColor','none','color',color.stim1,'linewidth',1)
trial_start_vector = zeros(1,tot_frames);
trial_start_vector(trials.trial_start_frames(trials.stim_type==2))=-200;
stem(trial_start_vector,'MarkerEdgeColor','none','MarkerFaceColor','none','color',color.stim2,'linewidth',1)

ylabel('all trials')

subplot(2,1,2)
hold on
trial_outcomes = {'correct','incorrect','miss','fa'};
for  i = 1:numel(trial_outcomes)
    trial_start_vector = zeros(1,tot_frames);
    trial_start_vector(trials.trial_start_frames(trials.(trial_outcomes{i})==1&trials.bias == 0)) = 1;
    stem(trial_start_vector.*(num_init_comp*plot_offset),'MarkerEdgeColor','none','MarkerFaceColor','none','color',color.([trial_outcomes{i},'_trial']),'linewidth',1)
end

% mark training and text trials
plot(trials.trial_start_frames(train_trial_range),[0 0],'linewidth',2,'color','black')
plot(trials.trial_start_frames(test_trial_range),[0 0],'linewidth',2,'color',[.5 .5 .5])

% mark stim type
trial_start_vector = zeros(1,tot_frames);
trial_start_vector(trials.trial_start_frames(trials.stim_type==1))=-200;
stem(trial_start_vector,'MarkerEdgeColor','none','MarkerFaceColor','none','color',color.stim1,'linewidth',1)
trial_start_vector = zeros(1,tot_frames);
trial_start_vector(trials.trial_start_frames(trials.stim_type==2))=-200;
stem(trial_start_vector,'MarkerEdgeColor','none','MarkerFaceColor','none','color',color.stim2,'linewidth',1)

ylabel('exclu. biased')

tot_num_correct = numel(find(trials.correct==1));
tot_num_incorrect = numel(find(trials.incorrect==1));
tot_num_miss = numel(find(trials.miss==1));
tot_num_fa = numel(find(trials.fa==1));

suptitle(['Num. correct:' num2str(tot_num_correct), '  incorrect:' num2str(tot_num_incorrect),'  miss:' num2str(tot_num_miss),'  FA:' num2str(tot_num_fa)])

%% performance rate
figure('name','performance overview')
subplot(2,1,1)
perform = struct();
trial_outcomes = {'miss','correct','incorrect','fa'};
this_trial_range = train_trial_range;
this_plot_color = [];
for  i = 1:numel(trial_outcomes)
    for stim_type = 1:2
        this_stim_num_trials = numel(find(trials.stim_type(this_trial_range(1):this_trial_range(end)) == stim_type));
        this_fd_name = [trial_outcomes{i} '_stim_' num2str(stim_type)];
        perform.(this_fd_name) = ...
            numel(find(trials.(trial_outcomes{i})(this_trial_range(1):this_trial_range(end))==1 & trials.stim_type(this_trial_range(1):this_trial_range(end)) == stim_type))/this_stim_num_trials;
        this_plot_color = [this_plot_color;color.(['st_' this_fd_name])];
    end
end
scatter_cmp_conditions(perform,'performance',1,this_plot_color,'connect_scatter',0,'VeryBriefXlabel',1);
ylabel('Fraction of trials (train)')

subplot(2,1,2)
perform = struct();
trial_outcomes = {'miss','correct','incorrect','fa'};
this_trial_range = test_trial_range;
this_plot_color = [];
for  i = 1:numel(trial_outcomes)
    for stim_type = 1:2
        this_stim_num_trials = numel(find(trials.stim_type(this_trial_range(1):this_trial_range(end)) == stim_type));
        this_fd_name = [trial_outcomes{i} '_stim_' num2str(stim_type)];
        perform.(this_fd_name) = ...
            numel(find(trials.(trial_outcomes{i})(this_trial_range(1):this_trial_range(end))==1 & trials.stim_type(this_trial_range(1):this_trial_range(end)) == stim_type))/this_stim_num_trials;
        this_plot_color = [this_plot_color;color.(['st_' this_fd_name])];
    end
end
scatter_cmp_conditions(perform,'performance',1,this_plot_color,'connect_scatter',0,'VeryBriefXlabel',1);
ylabel('Fraction of trials (test)')

%% trial outcome around incorrect trials
figure
subplot(1,2,1)
hold on
fd_names = fields(trialoutcome_struct.('st_incorrect_stim_1'));
for  i = 1:numel(fd_names)
    this_fd_name = fd_names{i};
    this_traces = trialoutcome_struct.('st_incorrect_stim_1').(this_fd_name);
    this_trace = nanmean(this_traces,2);
    if contains(this_fd_name,'fa')
        plot(this_trace,'color',color.(this_fd_name),'linestyle','--','linewidth',2);
    else
        plot(this_trace,'color',color.(this_fd_name));
    end
end
ylabel('Fraction of trials')

subplot(1,2,2)
hold on
fd_names = fields(trialoutcome_struct.('st_incorrect_stim_2'));
for  i = 1:numel(fd_names)
    this_fd_name = fd_names{i};
    this_traces = trialoutcome_struct.('st_incorrect_stim_2').(this_fd_name);
    this_trace = nanmean(this_traces,2);
    if contains(this_fd_name,'fa')
        plot(this_trace,'color',color.(this_fd_name),'linestyle','--','linewidth',2);
    else
        plot(this_trace,'color',color.(this_fd_name));
    end
end
ylabel('Fraction of trials')
legend(cellfun(@(x)strrep(x,'_',' '),fd_names,'UniformOutput',false))

%% ====================  SIMPLE CALCIUM TRACES PLOTS ===========================
%% plot FOV
figure('name','fov','position',[100 100 1200 800])
subplot(1,2,1)
imagesc(fov_struct.com_fov)
colormap(gray)
axis square

subplot(1,2,2)
[CC,jsf] = plot_contours(sparse(double(fov_struct.cnm_A)),fov_struct.cnm_image,fov_struct.cnm_plot_options,1,[],[],[1 1 1]);

%% plot full calcium traces
figure('name','full traces'); hold on
plot_offset = 5;
tot_frames = length(cell_struct(1).deconvC_full);
num_init_comp = size(cell_struct,2);
% plot accepted components
for i = 1:num_init_comp
    %     plot(cell_struct(i).noisyC_full+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
    plot(cell_struct(i).deconvC_full+ double(i*plot_offset),'color',[.5 .5 .5],'linewidth',1.5)
end
xlabel('Frames')
ylabel('ROI index')
yticks([ 1:num_init_comp].*plot_offset)
yticklabels(1:num_init_comp)
xlim([0 tot_frames])

%% plot as point cloud
figure
xyzpoints = cell2mat({train_cell_struct(cell_idx_struct.('tex')).('sta_correct_stim_1')}');
% pcshow(xyzpoints)

%% plot AUCs
% how does cell response change over time
% some stim/port cells do not stay stim/port - what changed?
figure;
hold on
scatter(stimulusAUC,test_stimulusAUC)
plot([0 1],[0,1])
xlabel('Train stim AUC')
ylabel('Test stim AUC')

figure;
hold on
scatter(pcorrect,test_pcorrect)
plot([0 1],[0,1])
xlabel('Train p correct')
ylabel('Test p correct')
%% comapre stim and choice AUC
color.port = [78 145 252]./255;
color.stim = [255 166 50]./255;
bin_width = 0.05;
figure('position',[100 100 800 800]);

% outcome AUC for stimulus 1 trials
plt_outcomeAUC = test_outcomeAUC{1};
plt_stimulusAUC = test_stimulusAUC;
subplot(2,3,1)
hold on
scatter(plt_stimulusAUC,plt_outcomeAUC,'MarkerEdgeColor','black')
[correlation,p] = corrcoef(plt_stimulusAUC,plt_outcomeAUC);

scatter(plt_stimulusAUC(cell_idx_struct.('stim')),plt_outcomeAUC(cell_idx_struct.('stim')),'MarkerEdgeColor', color.stim)
scatter(plt_stimulusAUC(cell_idx_struct.('port')),plt_outcomeAUC(cell_idx_struct.('port')),'MarkerEdgeColor', color.port)
xlabel('Stimulus AUC (Choice1)')
ylabel('Choice AUC (Stim1)')
xlim([0 1])
ylim([0 1])

text(0.1,1,['r = ' num2str(correlation(1,2),'%10.3f') ' p = '  num2str(p(1,2),'%10.3f')],...
    'units','normalized', 'horizontalalignment','left', 'color','black')

axis square

subplot(2,3,2)
hold on
histogram(plt_stimulusAUC,'FaceColor',[.5 .5 .5],'EdgeColor','none','BinWidth',bin_width)
histogram(plt_stimulusAUC(cell_idx_struct.('stim')),'DisplayStyle','Stairs','EdgeColor',color.stim,'BinWidth',bin_width)
histogram(plt_stimulusAUC(cell_idx_struct.('port')),'DisplayStyle','Stairs','EdgeColor',color.port,'BinWidth',bin_width)
xlabel('Stimulus AUC(Choice1)')
ylabel('Num ROIs')
xlim([0 1])
axis square


subplot(2,3,3)
hold on
histogram(plt_outcomeAUC,'FaceColor',[.5 .5 .5],'EdgeColor','none','BinWidth',bin_width)
histogram(plt_outcomeAUC(cell_idx_struct.('stim')),'DisplayStyle','Stairs','EdgeColor',color.stim,'BinWidth',bin_width)
histogram(plt_outcomeAUC(cell_idx_struct.('port')),'DisplayStyle','Stairs','EdgeColor',color.port,'BinWidth',bin_width)
xlabel('Choice AUC(Stim1)')
ylabel('Num ROIs')
xlim([0 1])
axis square


% outcome AUC for stimulus 2 trials
plt_outcomeAUC = test_outcomeAUC{2};
plt_stimulusAUC = test_stimulusAUC;

subplot(2,3,4)
hold on
scatter(plt_stimulusAUC,plt_outcomeAUC,'MarkerEdgeColor','black')
[correlation,p] = corrcoef(plt_stimulusAUC,plt_outcomeAUC);

scatter(plt_stimulusAUC(cell_idx_struct.('stim')),plt_outcomeAUC(cell_idx_struct.('stim')),'MarkerEdgeColor', color.stim)
scatter(plt_stimulusAUC(cell_idx_struct.('port')),plt_outcomeAUC(cell_idx_struct.('port')),'MarkerEdgeColor', color.port)
xlabel('Stimulus AUC(Choice2)')
ylabel('Choice AUC(Stim2)')
xlim([0 1])
ylim([0 1])

text(0.1,1,['r = ' num2str(correlation(1,2),'%10.3f') ' p = '  num2str(p(1,2),'%10.3f')],...
    'units','normalized', 'horizontalalignment','left', 'color','black')

axis square

subplot(2,3,5)
hold on
histogram(plt_stimulusAUC,'FaceColor',[.5 .5 .5],'EdgeColor','none','BinWidth',bin_width)
histogram(plt_stimulusAUC(cell_idx_struct.('stim')),'DisplayStyle','Stairs','EdgeColor',color.stim,'BinWidth',bin_width)
histogram(plt_stimulusAUC(cell_idx_struct.('port')),'DisplayStyle','Stairs','EdgeColor',color.port,'BinWidth',bin_width)
xlabel('Stimulus AUC(Choice2)')
ylabel('Num ROIs')
xlim([0 1])
axis square


subplot(2,3,6)
hold on
histogram(plt_outcomeAUC,'FaceColor',[.5 .5 .5],'EdgeColor','none','BinWidth',bin_width)
histogram(plt_outcomeAUC(cell_idx_struct.('stim')),'DisplayStyle','Stairs','EdgeColor',color.stim,'BinWidth',bin_width)
histogram(plt_outcomeAUC(cell_idx_struct.('port')),'DisplayStyle','Stairs','EdgeColor',color.port,'BinWidth',bin_width)
xlabel('Choice AUC(Stim2)')
ylabel('Num ROIs')
xlim([0 1])
axis square


this_title = [opt.exp_name ' Stimulus and Choice AUC'];
suptitle(strrep(this_title,'_',' '))

%% compare coding dimension in trainig and test sets
temp_cd_struct = struct();
for idx = 1:numel(opt.idx_fields)
    this_idx_field = opt.idx_fields{idx};
    [this_proj_struct,this_cd,this_db,this_proj_struct_framewise,this_framewise_cd] = get_coding_direction(test_cell_struct(cell_idx_struct.(this_idx_field)),color,opt,...
        'ref_fds',{'sta_correct_stim_1','sta_correct_stim_2'},'correct_fds',{'st_correct_stim_1','st_correct_stim_2'});
    cd_struct = [];
    cd_struct.this_cd = this_cd;
    cd_struct.this_framewise_cd = this_framewise_cd;
    cd_struct.this_proj_struct_framewise = this_proj_struct_framewise;
    temp_cd_struct.(this_idx_field) = cd_struct;
end

figure;
for idx = 1:numel(opt.idx_fields)
    this_idx_field = opt.idx_fields{idx};
    subplot(4,numel(opt.idx_fields),idx); hold on
    plot(cell_idx_struct.(this_idx_field),train_cd_struct.(this_idx_field).this_cd,'-o','color','black')
    plot(cell_idx_struct.(this_idx_field),temp_cd_struct.(this_idx_field).this_cd,'-o','color',[.5 .5 .5]')
    plot(xlim,[0 0],'-','color','r')
    ylabel('CD')
    title(this_idx_field)
    
    
    subplot(4,numel(opt.idx_fields),idx+numel(opt.idx_fields)); hold on
    plot(cell_idx_struct.(this_idx_field),outcomeAUC{1}(cell_idx_struct.(this_idx_field)),'-o','color','black')
    plot(cell_idx_struct.(this_idx_field),test_outcomeAUC{1}(cell_idx_struct.(this_idx_field)),'-o','color',[.5 .5 .5]')
    plot(xlim,[.5 .5],'-','color','r')
    ylabel('outcome AUC 1')
    
    subplot(4,numel(opt.idx_fields),idx+numel(opt.idx_fields)*2); hold on
    plot(cell_idx_struct.(this_idx_field),outcomeAUC{2}(cell_idx_struct.(this_idx_field)),'-o','color','black')
    plot(cell_idx_struct.(this_idx_field),test_outcomeAUC{2}(cell_idx_struct.(this_idx_field)),'-o','color',[.5 .5 .5]')
    plot(xlim,[.5 .5],'-','color','r')
    ylabel('outcome AUC 2')
    
    subplot(4,numel(opt.idx_fields),idx+numel(opt.idx_fields)*3); hold on
    plot(cell_idx_struct.(this_idx_field),stimulusAUC(cell_idx_struct.(this_idx_field)),'-o','color','black')
    plot(cell_idx_struct.(this_idx_field),test_stimulusAUC(cell_idx_struct.(this_idx_field)),'-o','color',[.5 .5 .5]')
    plot(xlim,[.5 .5],'-','color','r')
    ylabel('stimulus AUC')
    
    
    xlabel('global ROI idx')
end
%% plot sta traces - compare train and test
correct_trace_names  = {'st_correct_stim_1','st_correct_stim_2'};
plot_idx_field = {'port'};
plot_offset = 2;
ylimit = [-1 3];
for i = 1:numel(plot_idx_field)
    this_idx_field = plot_idx_field{i};
    x_ticks = opt.frame_range;
    num_cells = numel(cell_idx_struct.(this_idx_field));
    num_cols = num_cells;
    num_rows = 2;
    
    figure('name','sta traces(shadederror)','units','normalized','outerposition',[0 0 1 1])
    for c = 1:num_cells
        
        subtightplot(num_rows,num_cols,c)
        for p = 1:numel(correct_trace_names)
            this_plot_fd = correct_trace_names{p};
            hold on
            this_traces = cell2mat({train_cell_struct(cell_idx_struct.(this_idx_field)(c)).(this_plot_fd)}')'+(p-1)*plot_offset;
            shadedErrorBar(x_ticks,mean(this_traces,1),...
                std(this_traces,[],1),{'color',color.(this_plot_fd),'linewidth',2},0.1)
            %             plot(this_traces'+(p-1)*plot_offset,'color',color.(this_plot_fd));
            
            this_incorrect_fd = strrep(this_plot_fd,'_correct_','_incorrect_');
            this_traces = cell2mat({train_cell_struct(cell_idx_struct.(this_idx_field)(c)).(this_incorrect_fd)}')'+(p-1)*plot_offset;
            shadedErrorBar(x_ticks,mean(this_traces,1),...
                std(this_traces,[],1),{'color',color.(this_incorrect_fd),'linewidth',2},0.1)
            %             plot(this_traces'+(p-1)*plot_offset,'color',color.(this_incorrect_fd));
            
        end
        
        text(0.05,1,['gROI ' num2str(cell_idx_struct.(this_idx_field)(c)) ],'units','normalized', 'horizontalalignment','left', 'color','black')
        text(0.05,.9,['outcome AUC1 ' num2str(outcomeAUC{1}(cell_idx_struct.(this_idx_field)(c)),'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')
        text(0.05,.8,['outcome AUC2 '  num2str(outcomeAUC{2}(cell_idx_struct.(this_idx_field)(c)),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
        text(0.05,.7,['stimulus AUC ' num2str(stimulusAUC(cell_idx_struct.(this_idx_field)(c)),'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')
        
        ylim(ylimit)
        ylabel([this_idx_field ' ROI ' num2str(c)])
        xlim([x_ticks(1),x_ticks(end)])
        set(gca,'XTick',[])
        plot([1 1].*opt.withold_frames_adj(1),ylim,':','color',[.5 .5 .5],'linewidth',1)
        plot([1 1].*opt.gocue_frame_adj,ylim,'color',[.5 .5 .5],'linewidth',1)
        
        subtightplot(num_rows,num_cols,c+num_cols)
        for p = 1:numel(correct_trace_names)
            this_plot_fd = correct_trace_names{p};
            hold on
            this_traces = cell2mat({test_cell_struct(cell_idx_struct.(this_idx_field)(c)).(this_plot_fd)}')'+(p-1)*plot_offset;
            shadedErrorBar(x_ticks,mean(this_traces,1),...
                std(this_traces,[],1),{'color',color.(this_plot_fd),'linewidth',2},0.1)
            %             plot(this_traces'+ (p-1)*plot_offset,'color',color.(this_plot_fd));
            
            
            this_incorrect_fd = strrep(this_plot_fd,'_correct_','_incorrect_');
            this_traces = cell2mat({test_cell_struct(cell_idx_struct.(this_idx_field)(c)).(this_incorrect_fd)}')'+(p-1)*plot_offset;
            shadedErrorBar(x_ticks,mean(this_traces,1),...
                std(this_traces,[],1),{'color',color.(this_incorrect_fd),'linewidth',2},0.1)
            %             plot(this_traces'+(p-1)*plot_offset,'color',color.(this_incorrect_fd));
            
        end
        
        text(0.1,1,['gROI ' num2str(cell_idx_struct.(this_idx_field)(c))],...
            'units','normalized', 'horizontalalignment','left', 'color','black')
        
        text(0.05,.9,['outcome AUC1 ' num2str(test_outcomeAUC{1}(cell_idx_struct.(this_idx_field)(c)),'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')
        text(0.05,.8,['outcome AUC2 '  num2str(test_outcomeAUC{2}(cell_idx_struct.(this_idx_field)(c)),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
        text(0.05,.7,['stimulus AUC ' num2str(test_stimulusAUC(cell_idx_struct.(this_idx_field)(c)),'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')
        
        ylim(ylimit)
        ylabel([this_idx_field ' ROI ' num2str(c)])
        xlim([x_ticks(1),x_ticks(end)])
        set(gca,'XTick',[])
        plot([1 1].*opt.withold_frames_adj(1),ylim,':','color',[.5 .5 .5],'linewidth',1)
        plot([1 1].*opt.gocue_frame_adj,ylim,'color',[.5 .5 .5],'linewidth',1)
        
    end
end
%% plot all traces - grouped by trial
trace_names  = {{'st_correct_stim_1','st_incorrect_stim_1'},{'st_correct_stim_2','st_incorrect_stim_2'}};
plot_idx_field = {'port','outcome'};
ylimit = [0 3];
for i = 1:numel(plot_idx_field)
    this_idx_field = plot_idx_field{i};
    x_ticks = opt.frame_range;
    num_cells = numel(cell_idx_struct.(this_idx_field));
    num_cols = numel(trace_names);
    num_rows = 2;
    
    figure('name','cell avg traces(shadederror)','units','normalized','outerposition',[0 0 1 1])
    this_struct = train_cell_struct(cell_idx_struct.(this_idx_field));
    for p = 1:numel(trace_names)
        subplot(num_rows,num_cols,p)
        hold on
        for tr = 1:numel(trace_names{p})
            this_plot_fd = trace_names{p}{tr};
            num_trials = size(this_struct(1).(this_plot_fd),2);
            num_frames = size(this_struct(1).(this_plot_fd),1);
            this_traces_mat = nan(num_cells,num_trials,num_frames);
            for c = 1:num_cells
                this_traces_mat(c,:,:) = cell2mat({this_struct(c).(this_plot_fd)}')'; %[trials,frame]
            end
            this_traces = squeeze(nanmean(this_traces_mat,1));
%             shadedErrorBar(x_ticks,mean(this_traces,1),...
%                 std(this_traces,[],1),{'color',color.(this_plot_fd),'linewidth',2},0.1)
            plot(this_traces','color',color.(this_plot_fd))
        end
        xlim([x_ticks(1),x_ticks(end)])
        ylim(ylimit)
        set(gca,'XTick',[])
        ylabel('Cell avg.')
        title([ 'train ' strrep(strrep(this_plot_fd, 'st_correct_stim_','tex '),'st_incorrect_stim_', 'tex ')])
        axis square
        plot([1 1].*opt.withold_frames_adj(1),ylim,':','color',[.5 .5 .5],'linewidth',1)
        plot([1 1].*opt.gocue_frame_adj,ylim,'color',[.5 .5 .5],'linewidth',1)
    end
    
    
    this_struct = test_cell_struct(cell_idx_struct.(this_idx_field));
    for p = 1:numel(trace_names)
        subplot(num_rows,num_cols,p+num_cols)
        hold on
        for tr = 1:numel(trace_names{p})
            this_plot_fd = trace_names{p}{tr};
            num_trials = size(this_struct(1).(this_plot_fd),2);
            num_frames = size(this_struct(1).(this_plot_fd),1);
            this_traces_mat = nan(num_cells,num_trials,num_frames);
            for c = 1:num_cells
                this_traces_mat(c,:,:) = cell2mat({this_struct(c).(this_plot_fd)}')'; %[trials,frame]
            end
            this_traces = squeeze(nanmean(this_traces_mat,1));
%             shadedErrorBar(x_ticks,mean(this_traces,1),...
%                 std(this_traces,[],1),{'color',color.(this_plot_fd),'linewidth',2},0.1)
            plot(this_traces','color',color.(this_plot_fd))

        end
        xlim([x_ticks(1),x_ticks(end)])
        ylim(ylimit)
        set(gca,'XTick',[])
        ylabel('Cell avg.')
        title([ 'test ' strrep(strrep(this_plot_fd, 'st_correct_stim_','tex '),'st_incorrect_stim_', 'tex ')])
        axis square
        plot([1 1].*opt.withold_frames_adj(1),ylim,':','color',[.5 .5 .5],'linewidth',1)
        plot([1 1].*opt.gocue_frame_adj,ylim,'color',[.5 .5 .5],'linewidth',1)
    end
    subtitle(['Activity averaged across cells in each trial ' this_idx_field])
end
%% plot all traces - grouped by cell
trace_names  = {'st_correct_stim_1','st_incorrect_stim_1','st_correct_stim_2','st_incorrect_stim_2'};
trace_names  = {'st_correct_stim_1','st_correct_stim_2'};
trace_names  = {'st_correct_stim_1','st_incorrect_stim_1'};
trace_names  = {'st_correct_stim_2','st_incorrect_stim_2'};

plot_idx_field = {'lasso_tex2'};

for i = 1:numel(plot_idx_field)
    this_idx_field = plot_idx_field{i};
    x_ticks = opt.frame_range;
    num_cells = numel(cell_idx_struct.(this_idx_field));
    num_cols = ceil(sqrt(num_cells));
    num_rows = num_cols;
    
    figure('name','train sta traces(shadederror)','units','normalized','outerposition',[0 0 1 1])
    this_struct = train_cell_struct(cell_idx_struct.(this_idx_field));
    for c = 1:num_cells
        subtightplot(num_rows,num_cols,c)
        hold on
        
        for p = 1:numel(trace_names)
            this_plot_fd = trace_names{p};
            this_traces = cell2mat({this_struct(c).(this_plot_fd)}')';
            shadedErrorBar(x_ticks,mean(this_traces,1),...
                std(this_traces,[],1),{'color',color.(this_plot_fd),'linewidth',2},0.1)
        end
        text(0.1,1,['gROI ' num2str(cell_idx_struct.(this_idx_field)(c))],...
            'units','normalized', 'horizontalalignment','left', 'color','black')
        text(0.05,.9,['tex AUC1 ' num2str(texAUC{1}(cell_idx_struct.(this_idx_field)(c)),'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')
        text(0.05,.8,['tex AUC2 '  num2str(texAUC{2}(cell_idx_struct.(this_idx_field)(c)),'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
        text(0.05,.7,['stimulus AUC ' num2str(test_stimulusAUC(cell_idx_struct.(this_idx_field)(c)),'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')

        
        ylabel(['ROI ' num2str(c)])
        xlim([x_ticks(1),x_ticks(end)])
        set(gca,'XTick',[])
        axis square
        plot([1 1].*opt.withold_frames_adj(1),ylim,':','color',[.5 .5 .5],'linewidth',1)
        plot([1 1].*opt.gocue_frame_adj,ylim,'color',[.5 .5 .5],'linewidth',1)
        
    end
    suptitle(['train ', this_idx_field])
    
    
    figure('name','test sta traces(shadederror)','units','normalized','outerposition',[0 0 1 1])
    this_struct = test_cell_struct(cell_idx_struct.(this_idx_field));
    for c = 1:num_cells
        subtightplot(num_rows,num_cols,c)
        hold on
        
        for p = 1:numel(trace_names)
            this_plot_fd = trace_names{p};
            this_traces = cell2mat({this_struct(c).(this_plot_fd)}')';
            shadedErrorBar(x_ticks,mean(this_traces,1),...
                std(this_traces,[],1),{'color',color.(this_plot_fd),'linewidth',2},0.1)
        end
        text(0.1,1,['gROI ' num2str(cell_idx_struct.(this_idx_field)(c))],...
            'units','normalized', 'horizontalalignment','left', 'color','black')
        
        ylabel(['ROI ' num2str(c)])
        xlim([x_ticks(1),x_ticks(end)])
        set(gca,'XTick',[])
        axis square
        plot([1 1].*opt.withold_frames_adj(1),ylim,':','color',[.5 .5 .5],'linewidth',1)
        plot([1 1].*opt.gocue_frame_adj,ylim,'color',[.5 .5 .5],'linewidth',1)
    end
    suptitle(['test ', this_idx_field])
end

%% Show STA on maps (train)
figure('name','sta on fov','position',[200 200 1200 800])
plot_fds = {'st_correct_stim_1','st_incorrect_stim_1','st_miss_stim_1','st_correct_stim_2','st_incorrect_stim_2','st_miss_stim_2'};
plot_row = 2; %
plot_col = 3;
zlimit = [0 10];
for s = 1: numel(plot_fds)
    ax1 = subplot(plot_row,plot_col,s);
    value_field = plot_fds{s};
    this_title = strsplit(value_field,'_');
    this_title = [this_title{2} ' ' this_title{3} ' ' this_title{4}];
    plot_value_in_rois( train_cell_struct, [value_field '_amp'],[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',1,'zlimit',zlimit);
    title(this_title)
end
suptitle('ALL ROIS (train)')


zlimit = [0 10];
figure('name','tex1 sta on fov','position',[200 200 1200 800])
plot_row = 2; %
plot_col = 3;
for s = 1: numel(plot_fds)
    ax1 = subplot(plot_row,plot_col,s);
    value_field = plot_fds{s};
    this_title = strsplit(value_field,'_');
    this_title = [this_title{2} ' ' this_title{3} ' ' this_title{4}];
    plot_value_in_rois( train_cell_struct(cell_idx_struct.tex1), [value_field '_amp'],[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',1,'zlimit',zlimit);
    title(this_title)
end
suptitle('TEX1 ROIS (train)')

figure('name','tex2 sta on fov','position',[200 200 1200 800])
plot_row = 2; %
plot_col = 3;
for s = 1: numel(plot_fds)
    ax1 = subplot(plot_row,plot_col,s);
    value_field = plot_fds{s};
    this_title = strsplit(value_field,'_');
    this_title = [this_title{2} ' ' this_title{3} ' ' this_title{4}];
    plot_value_in_rois( train_cell_struct(cell_idx_struct.tex2), [value_field '_amp'],[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',1,'zlimit',zlimit);
    title(this_title)
end
suptitle('TEX2 ROIS (train)')


figure('name','port sta on fov','position',[200 200 1200 800])
plot_row = 2; %
plot_col = 3;
for s = 1: numel(plot_fds)
    ax1 = subplot(plot_row,plot_col,s);
    value_field = plot_fds{s};
    this_title = strsplit(value_field,'_');
    this_title = [this_title{2} ' ' this_title{3} ' ' this_title{4}];
    plot_value_in_rois( train_cell_struct(cell_idx_struct.port), [value_field '_amp'],[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',1,'zlimit',zlimit);
    title(this_title)
end
suptitle('PORT ROIS (train)')

figure('name','outcome sta on fov','position',[200 200 1200 800])
plot_row = 2; %
plot_col = 3;
for s = 1: numel(plot_fds)
    ax1 = subplot(plot_row,plot_col,s);
    value_field = plot_fds{s};
    this_title = strsplit(value_field,'_');
    this_title = [this_title{2} ' ' this_title{3} ' ' this_title{4}];
    plot_value_in_rois( train_cell_struct(cell_idx_struct.outcome), [value_field '_amp'],[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',1,'zlimit',zlimit);
    title(this_title)
end
suptitle('OUTCOME ROIS (train)')

%% Show STA on maps (test)
figure('name','sta on fov','position',[200 200 1200 800])
plot_fds = {'st_correct_stim_1','st_incorrect_stim_1','st_miss_stim_1','st_correct_stim_2','st_incorrect_stim_2','st_miss_stim_2'};
plot_row = 2; %
plot_col = 3;
zlimit = [0 15];
for s = 1: numel(plot_fds)
    ax1 = subplot(plot_row,plot_col,s);
    value_field = plot_fds{s};
    this_title = strsplit(value_field,'_');
    this_title = [this_title{2} ' ' this_title{3} ' ' this_title{4}];
    plot_value_in_rois( test_cell_struct, [value_field '_amp'],[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',1,'zlimit',zlimit);
    title(this_title)
end
suptitle('ALL ROIS (test)')


zlimit = [0 8];
figure('name','stim sta on fov','position',[200 200 1200 800])
plot_row = 2; %
plot_col = 3;
for s = 1: numel(plot_fds)
    ax1 = subplot(plot_row,plot_col,s);
    value_field = plot_fds{s};
    this_title = strsplit(value_field,'_');
    this_title = [this_title{2} ' ' this_title{3} ' ' this_title{4}];
    plot_value_in_rois( test_cell_struct(cell_idx_struct.stim), [value_field '_amp'],[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',1,'zlimit',zlimit);
    title(this_title)
end
suptitle('STIM ROIS (test)')

figure('name','port sta on fov','position',[200 200 1200 800])
plot_row = 2; %
plot_col = 3;
for s = 1: numel(plot_fds)
    ax1 = subplot(plot_row,plot_col,s);
    value_field = plot_fds{s};
    this_title = strsplit(value_field,'_');
    this_title = [this_title{2} ' ' this_title{3} ' ' this_title{4}];
    plot_value_in_rois( test_cell_struct(cell_idx_struct.port), [value_field '_amp'],[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',1,'zlimit',zlimit);
    title(this_title)
end
suptitle('PORT ROIS (test)')

figure('name','outcome sta on fov','position',[200 200 1200 800])
plot_row = 2; %
plot_col = 3;
for s = 1: numel(plot_fds)
    ax1 = subplot(plot_row,plot_col,s);
    value_field = plot_fds{s};
    this_title = strsplit(value_field,'_');
    this_title = [this_title{2} ' ' this_title{3} ' ' this_title{4}];
    plot_value_in_rois( test_cell_struct(cell_idx_struct.outcome), [value_field '_amp'],[256 256],ax1,'IF_NORM_PIX',0,'IF_CONTOUR',1,'zlimit',zlimit);
    title(this_title)
end
suptitle('OUTCOME ROIS (test)')