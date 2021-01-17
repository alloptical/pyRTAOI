%% mix of different check plots
shuf_accuracy = cell2mat(stim_struct.shuf_classif_accuracy');
num_shuf = size(shuf_accuracy,2);
num_frames = size(shuf_accuracy,1);
figure
hold on
for f = 1:num_frames
    scatter(ones(1,num_shuf).*f,shuf_accuracy(f,:),'MarkerEdgeColor','white','MarkerFaceColor',[.5 .5 .5])
end
title('stim')
%%
figure;
hold on
plot(stim_struct.framewise_hr,'color','black')
plot(stim_struct.framewise_fa,'color',[.5 .5 .5])
ylim([0 1])
title('stim decoder')

%%
figure;
hold on
plot(choice_struct.framewise_hr,'color','black')
plot(choice_struct.framewise_fa,'color',[.5 .5 .5])
ylim([0 1])
title('choice decoder')

%%
figure
subplot(2,4,5); hold on
histogram(stim1AUC_oppo_zscore,'FaceColor','none','EdgeColor',trial_color.L1,'DisplayStyle','stairs')
histogram(stim2AUC_oppo_zscore,'FaceColor','none','EdgeColor',trial_color.L2,'DisplayStyle','stairs')
legend('stim1 tuned','stim2 tuned')
axis square

subplot(2,4,6); hold on
scatter(stim1AUC_oppo_zscore,stim2AUC_oppo_zscore,'MarkerEdgeColor','black')
[ r,~,~,~,~,~,p ] = plot_fit_linear( stim1AUC_oppo_zscore,stim2AUC_oppo_zscore,[-5 5],'r');
axis square
xlabel('stim1 (oppo) auc')
ylabel('stim2 (oppo) auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])

subplot(2,4,7); hold on
scatter(correct_stimulusAUC_zscore,stim1AUC_oppo_zscore,'MarkerEdgeColor','black')
scatter(correct_stimulusAUC_zscore(cell_idx_struct.stim1_tuned),stim1AUC_oppo_zscore(cell_idx_struct.stim1_tuned),'MarkerEdgeColor',trial_color.L1)
scatter(correct_stimulusAUC_zscore(cell_idx_struct.stim2_tuned),stim1AUC_oppo_zscore(cell_idx_struct.stim2_tuned),'MarkerEdgeColor',trial_color.L2)
[ r,~,~,~,~,~,p ] = plot_fit_linear( correct_stimulusAUC_zscore,stim1AUC_oppo_zscore,[-5 5],'r');
axis square
xlabel('stim auc')
ylabel('stim1 (oppo) auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])

subplot(2,4,8); hold on
scatter(correct_stimulusAUC_zscore,stim2AUC_oppo_zscore,'MarkerEdgeColor','black')
scatter(correct_stimulusAUC_zscore(cell_idx_struct.stim1_tuned),stim2AUC_oppo_zscore(cell_idx_struct.stim1_tuned),'MarkerEdgeColor',trial_color.L1)
scatter(correct_stimulusAUC_zscore(cell_idx_struct.stim2_tuned),stim2AUC_oppo_zscore(cell_idx_struct.stim2_tuned),'MarkerEdgeColor',trial_color.L2)

[ r,~,~,~,~,~,p ] = plot_fit_linear( correct_stimulusAUC_zscore,stim2AUC_oppo_zscore,[-5 5],'r');
axis square
xlabel('stim auc')
ylabel('stim2 (oppo) auc')
title(['r=' num2str(r,'%10.2f') '; P= ' num2str(p,'%10.4f')])

%%
i = 10;
deconv_trace = cnm_struct(cell_struct(i).cnm_idx).deconvC_full;
online_trace = cnm_struct(cell_struct(i).cnm_idx).onlineC_full;
figure;plot(deconv_trace); hold on; plot(online_trace)
%%
%% Plot STAs 
% same amplitude for contra deflection (left side of animal); different for ipsi
cmp_fds = {'stim_2_var_5_incorrect','stim_2_var_5_correct'};
figure('name','baseline session stim sta traces','units','normalized','outerposition',[0 0 1 1])
plot_num_cells = num_cells;
num_plot_cols = 8;
avg_frame_range = opt.sta_avg_frames;
num_plot_rows = ceil(plot_num_cells/num_plot_cols);
for i = 1:plot_num_cells
    subtightplot(num_plot_rows,num_plot_cols,i)
    hold on
    % correct trials

    plot(cell_struct(i).(cmp_fds{1}),'color',trial_color.L1,'linewidth',1)
    plot(cell_struct(i).(cmp_fds{2}),'color',trial_color.L2,'linewidth',1)

    xlim([0 opt.sta_pre_frames+opt.sta_post_frames+1])
    ylim([-2 20])
    axis square
    
    text(1,1,['Cell ' num2str(i) ' (ROI ' num2str(cell_struct(i).cnm_idx) ')'],'units','normalized','color','black','Horizontalalignment','right','VerticalAlignment','top')

    if(~isempty(find(opsin_positive_idx==i)))
        text(1,1,['Cell ' num2str(i)],'units','normalized','color','r','Horizontalalignment','right','VerticalAlignment','top')
    end
    box off
    
    % mark tuned cells
    if( cell_struct(i).is_tuned)
        box on
        this_color = trial_color.(['L' num2str(cell_struct(i).pref)]);
        set(gca,'XColor',this_color,'YColor',this_color,'linewidth',3)
    end
    
    % mark time
    plot([1,1].*opt.sta_stimon_frame, ylim,'color','black','linestyle',':')
    if ~ opt.flag_use_peak
        plot([avg_frame_range(1),avg_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    else
        plot([peak_frame_range(1),peak_frame_range(end)],[0,0],'color',[.5 .5 .5],'linewidth',2)
    end
    
%     % show auc
%     text(0.05,.8,['stim auc '  num2str(cell_struct(i).correct_stimAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
%     text(0.05,.7,['choice auc '  num2str(cell_struct(i).choiceAUC_zscore,'%10.2f') ],'units','normalized', 'horizontalalignment','left', 'color','black')
    
    % modify x ticks
    xaxisvalues = [0:45:opt.sta_pre_frames+opt.sta_post_frames];
    xticks(xaxisvalues)
    xticklabels(arrayfun(@(x){num2str(x)},(xaxisvalues-opt.sta_trialon_frame)./opt.frame_rate))
    
end
suptitle(strrep(strrep(caiman_file,'.mat',''),'_',' '))
%%
cmp_fds = {'stim_1_var_4_correct','stim_2_var_5_correct'};
figure('name','trial amplitude','units','normalized','outerposition',[0 0 1 1])
plot_num_cells = num_cells;
num_plot_cols = 10;
avg_frame_range = opt.sta_avg_frames;
num_plot_rows = ceil(plot_num_cells/num_plot_cols);
for i = 1:plot_num_cells
    subplot(num_plot_rows,num_plot_cols,i)
    values = struct();
    values.stim1 = mean(cell_struct(i).(cmp_fds{1})(avg_frame_range,:),1);
    values.stim2 = mean(cell_struct(i).(cmp_fds{2})(avg_frame_range,:),1);
    scatter_cmp_conditions(values,[],...
        1,[0 0 0;0 0 0],'connect_scatter',0,'VeryBriefXlabel',1,'ShowMeanInXlabel',1,...
        'add_jitter',1,'add_yjitter',0,'plot_stats',1);
    axis square
    ylabel(['Cell ' num2str(i) ' Amp.'])
end
suptitle([strrep(cmp_fds{1},'_','') ' vs ' strrep(cmp_fds{2},'_','')])

%%
cmp_fds = {'stim_1_var_4_incorrect','stim_1_var_4_correct'};
figure('name','trial amplitude','units','normalized','outerposition',[0 0 1 1])
plot_num_cells = num_cells;
num_plot_cols = 10;
avg_frame_range = opt.sta_avg_frames;
num_plot_rows = ceil(plot_num_cells/num_plot_cols);
for i = 1:plot_num_cells
    subplot(num_plot_rows,num_plot_cols,i)
    values = struct();
    values.incorrect = mean(cell_struct(i).(cmp_fds{1})(avg_frame_range,:),1);
    values.correct = mean(cell_struct(i).(cmp_fds{2})(avg_frame_range,:),1);
    scatter_cmp_conditions(values,[],...
        1,[0 0 0;0 0 0],'connect_scatter',0,'VeryBriefXlabel',1,'ShowMeanInXlabel',1,...
        'add_jitter',1,'add_yjitter',0,'plot_stats',1);
    axis square
    ylabel(['Cell ' num2str(i) ' Amp.'])
end
suptitle([strrep(cmp_fds{1},'_','') ' vs ' strrep(cmp_fds{2},'_','')])

%% response difference correct - incorrect
% response difference in stim1 - stim2

