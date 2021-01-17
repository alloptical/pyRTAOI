function [decod_accuracy,overall] = get_expected_online_decod_accuracy(decod_proj_struct,plot_trial_indices,monitor_frames)
dummy_stim.tex1 = [];

detect_fun1 = @(x)arrayfun(@(y)find(x(y,monitor_frames(1):monitor_frames(end)-1)<0,1,'first'),1:size(x,1),'un',0);
detect_fun2 = @(x)arrayfun(@(y)find(x(y,monitor_frames(1):monitor_frames(end)-1)>0,1,'first'),1:size(x,1),'un',0);

dummy_stim.tex1 = structfun(@(x)detect_fun1(x),decod_proj_struct,'un',0);
dummy_stim.tex2 = structfun(@(x)detect_fun2(x),decod_proj_struct,'un',0);

dummy_indices.stim_1_dummyphoto_correct = plot_trial_indices.('stim_1_correct')(cellfun(@(x)~isempty(x),dummy_stim.tex1.stim_1_correct));
dummy_indices.stim_1_dummyphoto_incorrect = plot_trial_indices.('stim_1_incorrect')(cellfun(@(x)~isempty(x),dummy_stim.tex1.stim_1_incorrect));
dummy_indices.stim_2_dummyphoto_correct = plot_trial_indices.('stim_2_correct')(cellfun(@(x)~isempty(x),dummy_stim.tex2.stim_2_correct));
dummy_indices.stim_2_dummyphoto_incorrect = plot_trial_indices.('stim_2_incorrect')(cellfun(@(x)~isempty(x),dummy_stim.tex2.stim_2_incorrect));

decod_accuracy = struct();
decod_accuracy.CR =1- numel(dummy_indices.stim_1_dummyphoto_correct)./numel(plot_trial_indices.stim_1_correct);
decod_accuracy.FA = numel(dummy_indices.stim_1_dummyphoto_incorrect)./numel(plot_trial_indices.stim_1_incorrect);
decod_accuracy.Hit = 1-numel(dummy_indices.stim_2_dummyphoto_correct)./numel(plot_trial_indices.stim_2_correct);
decod_accuracy.Miss = numel(dummy_indices.stim_2_dummyphoto_incorrect )./numel(plot_trial_indices.stim_2_incorrect);
overall = mean(cell2mat(struct2cell(decod_accuracy)));

end

