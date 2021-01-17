%% load data
data_folder = 'D:\pyRTAOI_data\stim_at_fixed_frames\GCaMP6f';
file_list = dir(data_folder);
file_list = file_list(cellfun(@(x)contains(x,'.mat'),{file_list.name}));
file_names = {file_list.name};
num_files = numel(file_names);
%% parameters
hsv = colormap(hsv);
hsv = hsv(2:end-3,:);
opt.frame_rate = 30;
opt.N = 2;
opt.all_reset_frames = 30;
opt.pre_exp_frames = 0;
opt.opt.window_size = 60;
opt.sta_pre_frames = 30;
opt.sta_post_frames = 60;
opt.sta_baseline_frames = 30;
opt.window_size = 60;
opt.sta_avg_frames = 30;
opt.sta_thresh = 1;
opt.trial_color = [];
num_trials = 20;
opt.trial_color = zeros(num_trials,3);
opt.tint_factor = 0.7;
opt.hsv = hsv;
close
indices = round(linspace(1,size(hsv,1),num_trials));
for i = 1:num_trials
opt.trial_color(i,:) = tint(hsv(indices(i),:),opt.tint_factor);
end
opt.fig_save_path = 'D:\pyRTAOI_data\stim_at_fixed_frames\GCaMP6f\plots';
%% process data - done skip this
all_data = struct();
for f = 1:num_files
    this_file = [data_folder filesep file_names{f}];
    disp(['Processing ' file_names{f} ' ...'])
    this_save_name = strrep(file_names{f},'.mat','');
    this_caiman_data = load(this_file);
    [ cell_struct ] = proc_OnlineProcFixedFrames( this_caiman_data, this_save_name,opt);
    all_data(f).cell_struct = cell_struct;
    all_data(f).proc_time = this_caiman_data.tottime;
    all_data(f).file_name = this_file;
end
% save('D:\pyRTAOI data\stim_at_fixed_frames\GCaMP6f\procData\20181008_gcamp6f_stimfixedframes_data','all_data')
% save('D:\pyRTAOI data\stim_at_fixed_frames\GCaMP6f\procData\20181008_gcamp6f_stimfixedframes_opt','opt')
%% load from saved struct 
load('D:\pyRTAOI_data\stim_at_fixed_frames\GCaMP6f\procData\20181008_gcamp6f_stimfixedframes_data')
load('D:\pyRTAOI_data\stim_at_fixed_frames\GCaMP6f\procData\20181008_gcamp6f_stimfixedframes_opt')
%% timing
pro_time = cell2mat({all_data.proc_time}).*1000;
figure('name','process time per frame');
hold on
histogram(pro_time,'facecolor',[.5 .5 .5],'edgecolor','none','normalization','probability');
axis square
avg_time = mean(pro_time);
med_time = median(pro_time);
frame_time = 33.3;
plot([frame_time frame_time],[0 0.3],':','color','r','linewidth',2)
plot([avg_time avg_time],[0 0.3],':','color',[.3 .3 .3],'linewidth',2)
plot([med_time med_time],[0 0.3],':','color',[.7 .7 .7],'linewidth',2)

xlabel('Procesing time (ms)')
ylabel('Fraction of frames')
text(1,1,['Avg. proc. time: ' num2str(avg_time,3) 'ms'],'units','normalized', 'horizontalalignment','right','color',[.3 .3 .3])
text(1,.95,['Med. proc. time: ' num2str(med_time,3) 'ms'],'units','normalized', 'horizontalalignment','right','color',[.7 .7 .7])
text(1,.9,'Frame rate: 30 Hz','units','normalized', 'horizontalalignment','right','color','r')
xlim([0 50])
% export_fig  D:\pyRTAOI_data\stim_at_fixed_frames\GCaMP6f\plots\proc_time.pdf -painters 

%% plot opsin map
plot_exp_idx = 1;

figure
ax = subplot(1,2,1)
plot_value_in_rois( all_data(1).cell_struct, [],[256 256],ax,...
    'IF_NORM_PIX',1,'IF_CONTOUR',0,'zlimit',[0 1],'colorlut',flipud(gray));
% set(gca,'Ydir','reverse')

title('all cells')

ax = subplot(1,2,2)
plot_value_in_rois( all_data(1).cell_struct, 'opsin_positive',[256 256],ax,...
    'IF_NORM_PIX',1,'IF_CONTOUR',0,'zlimit',[0 1],'colorlut',flipud(gray));
set(gca,'Ydir','reverse')

title('opsin positive cells')

%% photostim responses
sta_traces = struct();
sta_traces.target = [];
sta_traces.nontarget = [];

all_target_idx = cellfun(@(x)find(extractfield(x,'num_trials')>0),{all_data.cell_struct},'UniformOutput',false)';
all_nontarget_idx = cellfun(@(x)find(extractfield(x,'num_trials')==0),{all_data.cell_struct},'UniformOutput',false)';
all_sta_traces = cellfun(@(x){x.sta_trace}',{all_data.cell_struct},'UniformOutput',false)';

for f = 1:num_files
    sta_traces.target = [sta_traces.target;all_sta_traces{f}(all_target_idx{f})];
    sta_traces.nontarget = [sta_traces.nontarget;all_sta_traces{f}(all_nontarget_idx{f})];
end
sta_traces.target = cell2mat(sta_traces.target')';
sta_traces.nontarget = cell2mat(sta_traces.nontarget')';

figure;
ylmit = [-10 25];
subplot(1,2,1)
hold on
plot(sta_traces.target','color',[.7 .7 .7])
plot(median(sta_traces.target,1),'black','linewidth',2)
axis square
ylim(ylmit)
title(['targeted mcherry+: ', num2str(length(cell2mat(all_target_idx'))),'rois'])
subplot(1,2,2)

hold on
plot(sta_traces.nontarget','color',[.7 .7 .7])
plot(median(sta_traces.nontarget,1),'black','linewidth',2)
ylim(ylmit)
title(['not targeted:'  num2str(length(cell2mat(all_nontarget_idx'))),'rois'])
axis square
% export_fig  D:\pyRTAOI_data\stim_at_fixed_frames\GCaMP6f\plots\all_sta_black_median.pdf -painters 
%% photostim amplitude
sta_amp = struct();
sta_amp.target = [];
sta_amp.nontarget = [];
all_sta_amp = cellfun(@(x)extractfield(x,'sta_amp'),{all_data.cell_struct},'UniformOutput',false);
for f = 1:num_files
    sta_amp.target = [sta_amp.target;all_sta_amp{f}(all_target_idx{f})'];
    sta_amp.nontarget = [sta_amp.nontarget;all_sta_amp{f}(all_nontarget_idx{f})'];
end

figure
hold on
histogram( sta_amp.target,'facecolor',[231 73 24]./255,'edgecolor','none','normalization','count')
histogram( sta_amp.nontarget,'facecolor',[.5 .5 .5],'edgecolor','none','normalization','count')

%% number of cells detected; photo-responsive
num_cells_detected = cell2mat(cellfun(@(x)numel(x),{all_data.cell_struct},'UniformOutput',false));
num_cells_stimulated = cell2mat(cellfun(@(x)numel(find(extractfield(x,'num_trials')>0)),{all_data.cell_struct},'UniformOutput',false));
num_cells_responsive = cell2mat(cellfun(@(x)numel(find(extractfield(x,'is_photo')>0)),{all_data.cell_struct},'UniformOutput',false));
percent_cells_responsive = 100.*num_cells_responsive./num_cells_detected;

figure;
hold on
bar([num_cells_detected;num_cells_stimulated;num_cells_responsive]')
% scatter(1:num_files,percent_cells_responsive)
legend('detected','stimulated','responded')
xlabel('FOV idx')
ylabel('Number of cells')


figure;
subplot(1,2,1)
values = struct();
values.detected = num_cells_detected;
values.stimulated = num_cells_stimulated;
values.responsive = num_cells_responsive;
scatter_cmp_conditions(values,'Num cells',1,[],'connect_scatter',1)

subplot(1,2,2)
values = struct();
values.percent_cells_responsive = percent_cells_responsive;
scatter_cmp_conditions(values,'% responsive',1,[],'connect_scatter',0)
export_fig  D:\pyRTAOI_data\stim_at_fixed_frames\GCaMP6f\plots\num_response_summary.pdf -painters 




