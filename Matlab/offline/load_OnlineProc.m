% load data
caiman_data = load('D:\pyRTAOI data\stim at fixed frames\20180811_OG300_t_0004_DS_2.0_rtaoi_OnlineProc_163822.mat')


%% make data structure
cnm_struct = struct();
cnm_dims = caiman_data.cnm_dims;
cnm_image = reshape(caiman_data.cnm_b,cnm_dims);
cnm_A = full(caiman_data.cnm_A);
num_frames = caiman_data.t_cnm;

num_comp = size(cnm_A,2);
comp_shape =[cnm_dims(1),cnm_dims(2),num_comp];
cm = com(sparse(double(cnm_A)),cnm_dims(1),cnm_dims(2)); % center of mass
for i = 1:num_comp
    cnm_struct(i).shape = reshape(cnm_A(:,i),cnm_dims);
    cnm_struct(i).noisyC = caiman_data.noisyC(i+1,1:num_frames);
    cnm_struct(i).deconvC = caiman_data.cnm_C(i,1:num_frames);
    cnm_struct(i).centroid = cm(i,:);
    cnm_struct(i).frame_added = caiman_data.frame_added()
end
%% plot spatial components
com_fov = zeros(cnm_dims);
binary_fov = zeros(cnm_dims);
for i = 1:num_comp
    com_fov = com_fov+cnm_struct(i).shape;
end

cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = colormap(lines);
subplot(1,2,1)
imagesc(com_fov)
colormap(gray)
axis square
subplot(1,2,2)

plot_contours(sparse(double(cnm_A)),cnm_image,cnm_plot_options,1,[],[],[1 1 1]);
suptitle('Detected ROIs')

%% plot traces
figure; hold on
plot_offset = 5;
for i = 1:num_comp
    plot(cnm_struct(i).noisyC+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
    plot(cnm_struct(i).deconvC+i*plot_offset,'color','black','linewidth',1.5)
end
xlabel('Frames')
ylabel('ROI index')
xlim([500 size(cnm_struct(1).noisyC,2)])

% show frame of detection


%% stim triggered average
cell_struct = struct();
opt.N = 2;
opt.all_reset_frames = 30;
opt.pre_exp_frames = 0;
opt.opt.window_size = 60;
opt.sta_pre_frames = 30;
opt.sta_post_frames = 60;
opt.sta_baseline_frames = 30;
opt.window_size = 60;

stim_frames = caiman_data.photo_stim_frames_caiman+caiman_data.num_frames_init;
figure('name','sta traces')
num_plot_cols = 4;
num_plot_rows = ceil(num_comp/num_plot_cols);
for i = 1:num_comp
    subplot(num_plot_rows,num_plot_cols,i)
    this_cell_trace = cnm_struct(i).deconvC;
    [~,~,~,~,~,cell_struct(i).sta_traces,cell_struct(i).sta_trace] = make_sta_from_traces(this_cell_trace,stim_frames ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
    hold on
    plot(cell_struct(i).sta_traces','color',[.7 .7 .7])
    plot(cell_struct(i).sta_trace,'color',cnm_plot_options.roi_color(i,:),'linewidth',1.5)
    title(['ROI ' num2str(i)])
end


%% check
figure
plot(cnm_struct(3).noisyC)
%% cell-event triggered average

all_cell_indices = 1:num_comp; 

for i = 1:num_comp
    this_cell_trace = cnm_struct(i).deconvC;
    other_cell_traces = cell2mat({cnm_struct(all_cell_indices~=i).deconvC}');
    cell_struct(i).event_frames = get_event_frames_slidingSTD( cnm_struct(i).deconvC,opt.N,opt.all_reset_frames,opt.window_size,0);
    [~,~,~,~,~,~,cell_struct(i).event_trace] = make_sta_from_traces(this_cell_trace,cell_struct(i).event_frames ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
    [~,~,~,~,~,~,cell_struct(i).event_trig_avg] = make_sta_from_traces(other_cell_traces,cell_struct(i).event_frames ,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
end

num_plot_cols = 10;
num_plot_rows = ceil(num_comp/num_plot_cols);
temp_all_values = [cell2mat({cell_struct(:).event_trace})'; cell2mat({cell_struct(:).event_trig_avg}')];
ylimit = [floor(min(temp_all_values(:))), ceil(max(temp_all_values(:)))];
xlimit = [0 opt.sta_pre_frames+opt.sta_post_frames+1];
clear temp_all_values

figure('name','detected events');
for i = 1:num_comp
    subplot(num_plot_rows,num_plot_cols,i)
    hold on;
    plot(cell_struct(i).event_trig_avg','-');
    plot(cell_struct(i).event_trace,'color','black','linewidth',2)
    ylim(ylimit);
    xlim(xlimit);
    title(['ROI ' num2str(i) ',' num2str(numel(cell_struct(i).event_frames)) ' events'])
end


