%% 
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
