function [cnm_struct,cell_struct] = get_cell_struct(caiman_data,sens_stim_frames,opt)
IF_USE_FILTC = 1;
try
 IF_USE_FILTC =   opt.IF_USE_FILTC;
end
cnm_struct = struct();
cnm_dims = double(caiman_data.cnm_dims);
cnm_image = reshape(caiman_data.cnm_b,cnm_dims);
cnm_A = full(caiman_data.cnm_A);
num_frames = caiman_data.t_cnm;

num_comp = size(cnm_A,2);
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
    cnm_struct(i).onlineC = caiman_data.online_C(i+1,1:num_frames);
    cnm_struct(i).deconvC = caiman_data.cnm_C(i,1:num_frames);
    cnm_struct(i).centroid = cm(i,:);
    
    % set skipped frames to nan
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.cnm_C(i,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).deconvC_full = temp_trace;
    
end
%% plot spatial components and save to cell struct
com_fov = zeros(cnm_dims);
accepted_com_fov = zeros(cnm_dims);
accepted_idx = caiman_data.accepted_idx+1;
num_cells = numel(accepted_idx);

for i = 1:num_comp
    com_fov = com_fov+cnm_struct(i).shape;
end

for i = accepted_idx
    accepted_com_fov = com_fov+cnm_struct(i).shape;
end


cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = [colormap(lines);colormap(lines);colormap(lines)];
close

figure('name','fov','position',[100 100 1200 500])
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
    cell_struct(i).cnm_idx = this_idx;
    cell_struct(i).jsf = jsf(this_idx);
    try
    cell_struct(i).cnn_prediction = caiman_data.cnn_predictions(this_idx);
    catch
        cell_struct(i).cnn_prediction =1;
    end
    cell_struct(i).filtC = caiman_data.filt_C(i,1:num_frames);
end
%% Get stim triggered average 
for i = 1:num_cells
    this_cell_trace = zscore( cnm_struct(cell_struct(i).cnm_idx).deconvC_full);
    if IF_USE_FILTC
        this_online_trace = cell_struct(i).filtC; % online trace med filtered
    else
        this_online_trace = cnm_struct(cell_struct(i).cnm_idx).onlineC_full; % online trace without median filter
    end

    
    if ~isempty(sens_stim_frames)
    sta_sens_stim_frames =  sens_stim_frames+opt.gocue_frame;
    this_num_trials = numel(sta_sens_stim_frames );
    else
        this_num_trials = 0;
    end
    cell_struct(i).num_trials = this_num_trials;
    cell_struct(i).is_sensory = 0;
    cell_struct(i).is_offcell = 0;
    cell_struct(i).pref_orient = [];
    cell_struct(i).sta_amp = 0;
    cell_struct(i).sta_traces = [];
    cell_struct(i).sta_trace = [];
    cell_struct(i).accepted = 0;
    cell_struct(i).pref_orient = nan;
    if(~isempty(find(accepted_idx==i)))
        cell_struct(i).accepted = 1;
    end
    
    if(this_num_trials>0)
        % average across all stim types
        % using df 
%         [~,~,~,~,~,cell_struct(i).sta_traces,cell_struct(i).sta_trace] = make_sta_from_traces(this_cell_trace,this_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
          % using raw f
        [~,~,~,cell_struct(i).sta_traces,~,~,cell_struct(i).sta_trace] = make_sta_from_traces(this_cell_trace,sta_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);
        [~,~,~,cell_struct(i).raw_sta_traces,~,~,cell_struct(i).raw_sta_trace] = make_sta_from_traces(this_online_trace,sta_sens_stim_frames,opt.sta_pre_frames,opt.sta_post_frames,1:opt.sta_baseline_frames);

            cell_struct(i).sta_amp = mean(cell_struct(i).sta_trace(opt.sta_avg_frames));

        if  cell_struct(i).sta_amp > opt.sta_amp_thresh
            cell_struct(i).is_sensory = 1;
        end        
    end
    
end

end

