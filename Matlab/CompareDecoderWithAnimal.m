%% test stim decoder performance on hard trials
%% get projections on stim decoder
test_opt = choice_opt;
test_weights = stim_norm_weights;
test_thresh = stim_norm_thresh;
test_fd_names = {'stim_1_var_9_correct','stim_1_var_9_incorrect','stim_1_var_9_miss',...
                'stim_1_var_3_correct','stim_1_var_3_incorrect','stim_1_var_3_miss',...
                     'stim_1_var_5_correct','stim_1_var_5_incorrect','stim_1_var_5_miss',...
                     'stim_1_var_7_correct','stim_1_var_7_incorrect','stim_1_var_7_miss',...
                     'stim_1_var_4_correct','stim_1_var_4_incorrect','stim_1_var_4_miss',...%                      
                     'stim_2_var_2_correct','stim_2_var_2_incorrect','stim_2_var_2_miss',...
                     'stim_2_var_5_correct','stim_2_var_5_incorrect','stim_2_var_5_miss',...
                     'stim_2_var_4_correct','stim_2_var_4_incorrect','stim_2_var_4_miss',...
                     'stim_2_var_6_correct','stim_2_var_6_incorrect','stim_2_var_6_miss',...
                     'stim_2_var_9_correct','stim_2_var_9_incorrect','stim_2_var_9_miss'};

all_stim_proj_struct = struct();
[all_stim_proj_struct] = get_projections(cell_struct(cell_idx_struct.(opt.trigger_idx_fd)),test_weights,test_fd_names,...
    'proj_struct',all_stim_proj_struct,'bias',-test_thresh,'IS_CELL_STRUCT',1);

% plot_pop_vectors(all_stim_proj_struct,test_fd_names,1,opt,...
%         'ylimit',[-20 20],'plot_ylabel','Projection','plot_num_cols',6,'IF_PLOT_RAW_ONLY',0)
%% organise animal and decoder performance in the same structure
stim_decod_result = struct();
frames_to_avg = max(opt.sta_stimon_frame-opt.sta_baseline_frames+15,fa_opt.bin_size*(decod_struct.shuf_disc_frame-opt.sta_baseline_frames))+[0 15];
frames_to_avg = frames_to_avg(1):frames_to_avg(end);
for f = 1:numel(test_fd_names)
    fd = test_fd_names{f};
        this_type = strsplit(fd,'_');
        this_stim_type = str2double(this_type{2});
        this_var_type = str2double(this_type{4});
        this_outcome = this_type{end};
        stim_decod_result(f).animal_correct = double(strcmp(this_outcome,'correct'))*size(all_stim_proj_struct.(fd),1);
        stim_decod_result(f).animal_incorrect = double(strcmp(this_outcome,'incorrect'))*size(all_stim_proj_struct.(fd),1);
        stim_decod_result(f).decod_correct = nan;
        stim_decod_result(f).stim_type = this_stim_type;
        stim_decod_result(f).stim_var = this_var_type;
        stim_decod_result(f).num_trials = size(all_stim_proj_struct.(fd),1);

        stim_decod_result(f).proj_val = [];

        
    if ~isempty(all_stim_proj_struct.(fd))
        this_raw = 1+double( mean(all_stim_proj_struct.(fd)(:,frames_to_avg),2)<0);
        stim_decod_result(f).raw = this_raw;
        stim_decod_result(f).decod_correct = numel(find(this_raw ==this_stim_type));
        stim_decod_result(f).proj_val = mean(all_stim_proj_struct.(fd)(:,frames_to_avg),2);
        if size( trial_hist.reward_type.(fd),2)>1
            stim_decod_result(f).prev_reward = trial_hist.reward_type.(fd)(:,pre_num_trials);
        else
            stim_decod_result(f).prev_reward = trial_hist.reward_type.(fd)(pre_num_trials);
        end

    end
end

% collapse
stim_decod_result_collapse = struct();
count = 1;
for f = 1:3:numel(test_fd_names)-1
    stim_decod_result_collapse(count).animal_pc_correct = nansum(cell2mat({stim_decod_result(f:f+2).('animal_correct')}))/nansum(cell2mat({stim_decod_result(f:f+2).('num_trials')}));
    stim_decod_result_collapse(count).decod_pc_correct = nansum(cell2mat({stim_decod_result(f:f+2).('decod_correct')}))/nansum(cell2mat({stim_decod_result(f:f+2).('num_trials')}));
    
    stim_decod_result_collapse(count).decod_correct_in_animal_correct = stim_decod_result(f).('decod_correct')/stim_decod_result(f).('num_trials');
    stim_decod_result_collapse(count).decod_correct_in_animal_incorrect = stim_decod_result(f+1).('decod_correct')/stim_decod_result(f+1).('num_trials');

    stim_decod_result_collapse(count).stim_type = stim_decod_result(f).stim_type;
    stim_decod_result_collapse(count).stim_var = stim_decod_result(f).stim_var;
     [stim_decod_result_collapse(count).left_volt,stim_decod_result_collapse(count).right_volt] = get_volts(stim_decod_result(f).stim_type,stim_decod_result(f).stim_var);

    count = count+1;
end

%% check if the trials used for training look normal
plot_pop_vectors(all_stim_proj_struct,test_opt.fd_names,1,opt,...
        'plot_ylabel','Projection','plot_num_cols',2,'IF_PLOT_RAW_ONLY',0)
figure
values = struct();
values.stim1 = stim_decod_result(contains(test_fd_names,stim_fds{1})==1).proj_val ;
values.stim2 = stim_decod_result(contains(test_fd_names,stim_fds{2})==1).proj_val ;
scatter_cmp_conditions(values,[],...
    1,[0 0 0;0 0 0],'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1);
axis square
%% compare trial history (previous reward)
figure('name','trial history')
count = 1;
num_plots = round(numel(test_fd_names)/3);
for f = 1:3:numel(test_fd_names)-1
    subplot(2,num_plots/2,count)
    values = struct();
    values.correct = stim_decod_result(f).prev_reward;
    values.incorrect = stim_decod_result(f+1).prev_reward;
    scatter_cmp_conditions(values,[],...
        1,[0 0 0;0 0 0],'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,...
        'add_jitter',1,'add_yjitter',1);
    axis square
    ylabel('Prev reward')
    title(['stim ' num2str(stim_decod_result(f).stim_type) ', var ' num2str(stim_decod_result(f).stim_var)])
    count = count+1;
end

%% compare projection on decoder axis in correct and incorrect trials
figure('name','projection on decoder')
count = 1;
num_plots = round(numel(test_fd_names)/3);
for f = 1:3:numel(test_fd_names)-1
    subplot(2,num_plots/2,count)
    values = struct();
    values.correct = stim_decod_result(f).proj_val ;
    values.incorrect = stim_decod_result(f+1).proj_val ;
    scatter_cmp_conditions(values,[],...
        1,[0 0 0;0 0 0],'connect_scatter',0,'BriefXlabel',0,'ShowMeanInXlabel',1,'add_jitter',1);
    axis square
    ylabel('Proj. on decoder')
    title(['stim ' num2str(stim_decod_result(f).stim_type) ', var ' num2str(stim_decod_result(f).stim_var)])
    count = count+1;
end

%% stim type decoder performance in correct trials vs incorrect trials
figure;
num_types = round(numel(test_fd_names)/3);
subplot(3,1,1);hold on
plot(cell2mat({stim_decod_result_collapse.('decod_correct_in_animal_correct')}),'-s','color','black')
plot(cell2mat({stim_decod_result_collapse.('decod_correct_in_animal_incorrect')}),'-s','color',[.5 .5 .5])

plot(cell2mat({stim_decod_result_collapse.('decod_pc_correct')}),'-*','color','black')
plot(cell2mat({stim_decod_result_collapse.('animal_pc_correct')}),'-*','color',[.7 .7 .7])
xlim([0 num_types+1])
ylim([0 1])
plot([0 num_types+1],[.33 .33 ],'r')
legend('correct trials','incorrect trials','decoder','animal','Location','southeastoutside')
ylabel('Decoder accuracy')

subplot(3,1,2); hold on
b1=bar([extractfield(stim_decod_result_collapse,'stim_type')' ...
   extractfield(stim_decod_result_collapse,'stim_var')']);
b1(1).FaceColor=[0 0 0];b1(2).FaceColor=[.5 .5 .5];
legend('stim type', 'stim var','Location','southeastoutside')
xlim([0 num_types+1])
grid on
ylabel('Index')


subplot(3,1,3); hold on
b2=bar([extractfield(stim_decod_result_collapse,'left_volt')' ...
   extractfield(stim_decod_result_collapse,'right_volt')']);
b2(1).FaceColor=[0 0 0];b2(2).FaceColor=[.5 .5 .5];
legend('left volt', 'right volt','Location','southeastoutside')
xlim([0 num_types+1])
ylabel('Deflection amplitude')

xlabel('Stimulus type')
%% compare animal performance with decoder performance
figure;
subplot(1,3,1);hold on
scatter(cell2mat({stim_decod_result_collapse.animal_pc_correct }),cell2mat({stim_decod_result_collapse.decod_pc_correct }),...
    'MarkerEdgeColor','black','MarkerFaceColor','w');
plot([0 1],[0 1],'color',[.5 .5 .5])
axis square
xlabel('Animal performance (%correct)')
ylabel('Decoder performance (%correct)')

subplot(1,3,2);hold on
scatter(cell2mat({stim_decod_result_collapse.left_volt})-cell2mat({stim_decod_result_collapse.right_volt}),cell2mat({stim_decod_result_collapse.animal_pc_correct }),...
    'MarkerEdgeColor','black','MarkerFaceColor','None');
scatter(cell2mat({stim_decod_result_collapse.left_volt})-cell2mat({stim_decod_result_collapse.right_volt}),cell2mat({stim_decod_result_collapse.decod_pc_correct }),...
    'MarkerEdgeColor',[.5 .5 .5],'MarkerFaceColor','None');

axis square
ylabel('Animal performance (%correct)')
xlabel('Deflection difference (left-right)')

subplot(1,3,3);hold on
axis square
ylabel('Decoder performance (%correct)')
xlabel('Deflection difference (left-right)')


