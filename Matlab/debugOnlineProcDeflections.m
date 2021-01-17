%% lines deleted after the '=END=' of OnlineProcDeflections
% things tried but did not help



%% Test choice decoder performance
% get projections
decod_struct = choice_struct;
opt.pop_opt = choice_opt;
weights = fa_struct.transmat* decod_struct.B(2:end)';
thresh = decod_struct.thresh_fix;
norm_weights = weights;
norm_thresh = thresh;
[norm_weights,norm_thresh] = get_norm_weights(weights,thresh,fa_struct.mean,fa_struct.std);
test_opt = stim_opt;
test_opt.trial_color = trial_color;
test_opt.fd_names = {'stim_1_var_3_correct','stim_1_var_3_incorrect','stim_1_var_3_miss',...
                     'stim_1_var_4_correct','stim_1_var_4_incorrect','stim_1_var_4_miss',...
                     'stim_1_var_5_correct','stim_1_var_5_incorrect','stim_1_var_5_miss',...
                     'stim_1_var_7_correct','stim_1_var_7_incorrect','stim_1_var_7_miss',...
                     'stim_3_var_2_correct','stim_3_var_2_incorrect','stim_3_var_2_miss'...
                     'stim_2_var_2_correct','stim_2_var_2_incorrect','stim_2_var_2_miss'...
                     'stim_2_var_4_correct','stim_2_var_4_incorrect','stim_2_var_4_miss'...
                     'stim_2_var_5_correct','stim_2_var_5_incorrect','stim_2_var_5_miss'...
                     'stim_2_var_6_correct','stim_2_var_6_incorrect','stim_2_var_6_miss'};
proj_struct = struct();
 [proj_struct] = get_projections(cell_struct(cell_idx_struct.(fa_opt.idx_fields{1})),norm_weights,test_opt.fd_names,'proj_struct',proj_struct,'bias',-norm_thresh,'IS_CELL_STRUCT',1);

plot_pop_vectors(proj_struct,test_opt.fd_names,1,test_opt,...
        'plot_ylabel','Projection','plot_num_cols',3,'IF_PLOT_RAW_ONLY',1)
suptitle(['Choice decoder projections:' strrep(strrep(caiman_file,'.mat',''),'_',' ')])
export_fig([fig_save_path filesep 'DecodProject_' strrep(caiman_file,'.mat','.png')])
%% plot choice decoder accuracy
test_decod_struct = [];
num_compares = round(numel(test_opt.fd_names)/3);
for i = 1:num_compares
    IF_REVERSE = 0;
    this_correct_fd = test_opt.fd_names{3*(i-1)+1};
    this_incorrect_fd = test_opt.fd_names{3*(i-1)+2};
    this_cmp_fds = {this_correct_fd,this_incorrect_fd};
    this_stim_type = strsplit(this_correct_fd,'_');
    this_stim_type = this_stim_type{2};
    % reverse for stim type2
    if strcmp(this_stim_type,'2')||strcmp(this_stim_type,'3')
        IF_REVERSE = 1;
    end
    try
        [  test_decod_struct{i} ] =  get_binary_decoder_disc_time( proj_struct, decod_struct,...
            this_cmp_fds,choice_opt,'IF_FRAMEWISE',choice_opt.IF_FRAMEWISE,'threshold',0,'IF_REVERSE',IF_REVERSE);
        test_decod_struct{i}.correct_fd = this_correct_fd;
        test_decod_struct{i}.incorrect_fd = this_incorrect_fd;
        
        figure('name','decoder performance tests')
        plot_binary_decoder(test_decod_struct{i},choice_opt)
        suptitle(cellfun(@(x)strrep(x,'_',' '),this_cmp_fds, 'UniformOutput', false))
    catch
        test_decod_struct{i} = [];
    end
end


%% test stim decoder performance
opt.pop_opt = stim_opt;
weights = fa_struct.transmat* stim_struct.B(2:end)';
thresh = stim_struct.thresh_fix;
norm_weights = weights;
norm_thresh = thresh;
% [norm_weights,norm_thresh] = get_norm_weights(weights,thresh,fa_struct.mean,fa_struct.std);
test_opt = stim_opt;
test_opt.trial_color = trial_color;
test_opt.fd_names = {'stim_1_var_3_correct','stim_1_var_3_incorrect','stim_1_var_3_miss',...
                     'stim_1_var_4_correct','stim_1_var_4_incorrect','stim_1_var_4_miss',...
                     'stim_1_var_5_correct','stim_1_var_5_incorrect','stim_1_var_5_miss',...
                     'stim_1_var_7_correct','stim_1_var_7_incorrect','stim_1_var_7_miss',...
                     'stim_3_var_2_correct','stim_3_var_2_incorrect','stim_3_var_2_miss'...
                     'stim_2_var_2_correct','stim_2_var_2_incorrect','stim_2_var_2_miss'...
                     'stim_2_var_4_correct','stim_2_var_4_incorrect','stim_2_var_4_miss'...
                     'stim_2_var_5_correct','stim_2_var_5_incorrect','stim_2_var_5_miss'...
                     'stim_2_var_6_correct','stim_2_var_6_incorrect','stim_2_var_6_miss'};

proj_struct = struct();
 [proj_struct] = get_projections(cell_struct(cell_idx_struct.(fa_opt.idx_fields{1})),norm_weights,test_opt.fd_names,'proj_struct',proj_struct,'bias',-norm_thresh,'IS_CELL_STRUCT',1);

plot_pop_vectors(proj_struct,test_opt.fd_names,1,test_opt,...
        'ylimit',[-50 50],'plot_ylabel','Projection','plot_num_cols',3,'IF_PLOT_RAW_ONLY',1)
suptitle('Stim decoder projections')
export_fig([fig_save_path filesep 'DecodProject_' strrep(caiman_file,'.mat','.png')])

%% get choice decoder after projecting to the stim axis
decod_fds = { 'stim_1_var_5_correct','stim_1_var_5_incorrect'};
test_struct = struct();
choice_after_stim_struct =  get_binary_classifier( test_struct,proj_struct, choice_opt,...
    'IF_CROSSVAL',1,'IF_FRAMEWISE',choice_opt.IF_FRAMEWISE,'fd_names',decod_fds);
[ choice_after_stim_struct ] =  get_binary_decoder_disc_time(proj_struct, choice_after_stim_struct, ...
    decod_fds,choice_opt,'IF_FRAMEWISE',choice_opt.IF_FRAMEWISE,'threshold',0);
figure; 
hold on
plot_binary_decoder(choice_after_stim_struct,choice_opt)
suptitle('Choice-stim decoder')
export_fig([fig_save_path filesep 'ChoiceStimDecoderPerform_' strrep(caiman_file,'.mat','.png')])

%% get weights and thresh after the two decoders
weights = fa_struct.transmat* stim_struct.B(2:end)';
thresh = stim_struct.thresh_fix+choice_after_stim_struct.B(2)/choice_after_stim_struct.B(1);
norm_weights = weights;
norm_thresh = thresh;
[choice_proj_struct] = get_projections(cell_struct(cell_idx_struct.(fa_opt.idx_fields{1})),weights,test_opt.fd_names,'bias',-thresh,'IS_CELL_STRUCT',1);

plot_pop_vectors(choice_proj_struct,test_opt.fd_names,1,test_opt,...
       'ylimit',[-20 40],'plot_ylabel','proj to choice-stim','plot_num_cols',3,'IF_PLOT_RAW_ONLY',1)
suptitle('Choice-stim decoder')

%% plot accuracy - select tiral types to condition according to these plots
test_decod_struct = [];
num_compares = round(numel(test_opt.fd_names)/3);
for i = 1:num_compares
    IF_REVERSE = 0;
    this_correct_fd = test_opt.fd_names{3*(i-1)+1};
    this_incorrect_fd = test_opt.fd_names{3*(i-1)+2};
    this_cmp_fds = {this_correct_fd,this_incorrect_fd};
    this_stim_type = strsplit(this_correct_fd,'_');
    this_stim_type = this_stim_type{2};
    % reverse for stim type2
    if strcmp(this_stim_type,'2')||strcmp(this_stim_type,'3')
        IF_REVERSE = 1;
    end
    try
        [  test_decod_struct{i} ] =  get_binary_decoder_disc_time( choice_proj_struct, choice_after_stim_struct,...
            this_cmp_fds,choice_opt,'IF_FRAMEWISE',choice_opt.IF_FRAMEWISE,'threshold',0,'IF_REVERSE',IF_REVERSE);
        test_decod_struct{i}.correct_fd = this_correct_fd;
        test_decod_struct{i}.incorrect_fd = this_incorrect_fd;
        
        figure('name','decoder performance tests')
        plot_binary_decoder(test_decod_struct{i},choice_opt)
        suptitle(cellfun(@(x)strrep(x,'_',' '),this_cmp_fds, 'UniformOutput', false))
    catch
        test_decod_struct{i} = [];
    end
end

%% generate parameters for pyRTAOI population analysis
pop_params = struct();
pop_params.weights = norm_weights;
pop_params.thresh = norm_thresh;
pop_params.frames_enable_trigger = max([opt.sta_stimon_frame,fa_opt.bin_size*(choice_after_stim_struct.shuf_disc_frame-opt.sta_baseline_frames)]);
pop_params.condition_stim_type = [1,1,3,2,2,2]; % CHANGE THIS
pop_params.condition_stim_var  = [5,4,2,2,5,7]; % type and var need to match
pop_params.condition_type = {[105 104], [202 207 205]}; % match target ensembles 100*stim_type + stim_var

%% =================== MAKE OUTPUT FILE FOR PYRTAOI =======================
[output3,save_time3] = generate_cell_idx_file(cell_struct,cell_idx_struct,pop_params,opt);



%%

 %% combine threshold trials 
proj_struct.TS_L1 = [];
proj_struct.TS_L2 = [];
 for f = 1:numel(TS_L1_fds)
      proj_struct.TS_L1 = [ proj_struct.TS_L1;proj_struct.(TS_L1_fds{f})];
 end
 for f = 1:numel(TS_L2_fds)
      proj_struct.TS_L2= [ proj_struct.TS_L1;proj_struct.(TS_L2_fds{f})];
 end
 
choice_after_stim_struct =  get_binary_classifier( test_struct,proj_struct, choice_opt,'IF_CROSSVAL',1,'IF_FRAMEWISE',choice_opt.IF_FRAMEWISE,'fd_names',choice_opt.fd_names);
[choice_proj_struct] = get_projections(proj_struct,choice_after_stim_struct.B(:,2:end)',test_opt.fd_names,'bias',choice_after_stim_struct.B(:,1));

plot_pop_vectors(choice_proj_struct,test_opt.fd_names,1,test_opt,...
        'plot_ylabel','proj to choice-stim','plot_num_cols',3)
   
[ test_struct ] =  get_binary_decoder_disc_time( proj_struct, choice_after_stim_struct,...
    test_opt.fd_names,test_opt,'IF_FRAMEWISE',test_opt.IF_FRAMEWISE,'threshold',0);
figure;
subplot(1,2,1)
plot_proj_traces(proj_struct,test_opt.fd_names,test_opt,...
    'plot_ylabel','Projection')
title('Stim-choice decoder (TS.L1 vs TS.L2)')

hold on
subplot(1,2,2)
hold on
plot_binary_decoder(test_struct,test_opt)



%% =========================== CHECK PLOTS ======================================
% %% Plot STAs for trigger cells (cells to monitor)
% figure('name','trigger sta traces')
% num_plot_cols = 4;
% num_plot_rois = length(cell_idx_struct.(opt.trigger_idx_fd));
% num_plot_rows = ceil(num_plot_rois/num_plot_cols);
% plot_count = 1;
% for ii = 1:num_plot_rois
%     i = cell_idx_struct.(opt.trigger_idx_fd)(ii);
%     subtightplot(num_plot_rows,num_plot_cols,plot_count)
%     % plot traces
%     hold on
%     for t = 1:size(cell_struct(i).sta_traces,1)
%         plot(cell_struct(i).sta_traces(t,:),'color',opt.trial_color(t,:) ,'linewidth',1)
%     end
%     plot(cell_struct(i).sta_trace,'color',[.5 .5 .5],'linewidth',1.5)
%     
%     xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
%     set(gca,'xtick',[],'xcolor',[1 1 1])
%     axis square
%     
%     plot([opt.sta_pre_frames opt.sta_pre_frames],ylim,'color',[.5 .5 .5]) % start of withold window
%     plot([opt.sta_pre_frames opt.sta_pre_frames] + length(opt.withold_frames_adj),ylim,'color',[0 0 0]) % go-cue
% 
%     plot_count = plot_count+1;
%     
%     text(0.05,1,['ROI ' num2str(i)],'units','normalized', 'horizontalalignment','left', 'color','black')
%     text(0.05,.9,['sensory auc ' num2str(cell_struct(i).correct_stimAUC,'%10.2f')],'units','normalized', 'horizontalalignment','left', 'color','black')
%     text(0.05,.8,['zscore auc '  num2str(cell_struct(i).correct_stimAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
%     
% end
% suptitle([opt.trigger_idx_fd ' rois, stim-triggered response'])


% %% Show sensory cells on maps
% figure('name','pref. texture on fov');
% subplot(1,2,1)
% imagesc(com_fov)
% colormap(gray)
% colorbar('location','southoutside');
% axis square
% title('Detected ROIs')
% 
% ax1 = subplot(1,2,2)
% value_field = 'pref_tex';
% plot_value_in_rois( cell_struct, value_field,[256 256],ax1,'colorlut',[[1,1,1];opt.type_color],'IF_NORM_PIX',0,'IF_CONTOUR',1,'IF_SHOW_OPSIN',1,'zlimit',[0 4]);
% set(gca,'Ydir','reverse')
% title('Sensory cells (colored by pref. texture)')
% 
% %% Plot auc distributions
% figure('name','correct stim auc')
% subplot(1,2,1)
% hold on
% histogram(extractfield(cell_struct,'correct_stimAUC'),'facecolor',[.7 .7 .7],'edgecolor','none','BinWidth',.05)
% histogram(extractfield(cell_struct(extractfield(cell_struct,'is_tuned')==1),'correct_stimAUC'),'facecolor','none','edgecolor',[0,0,1],'BinWidth',.05)
% xlabel('Tex response auc')
% axis square
% 
% subplot(1,2,2)
% hold on
% histogram(extractfield(cell_struct,'correct_stimAUC_zscore'),'facecolor',[.7 .7 .7],'edgecolor','none','BinWidth',.5)
% histogram(extractfield(cell_struct(extractfield(cell_struct,'is_tuned')==1),'correct_stimAUC_zscore'),'facecolor','none','edgecolor',[0,0,1],'BinWidth',.5)
% xlabel('Tex response auc (zscore)')
% axis square
% 
% 

