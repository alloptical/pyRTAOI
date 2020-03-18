function [output,save_time,full_file_save_path] = generate_cell_idx_file(cell_struct,cell_idx_struct,pop_params,opt)
% generate file for pyRTAOI photostim configuration
% generate centroid image for targets
% save output to new folder
save_time = datestr(now,'yyyymmdd_HHMM');
file_save_name = [opt.exp_name '_OutputParams_' save_time ];
target_img_save_name = [opt.exp_name '_OutputTargets_' save_time ];
trigger_img_save_name = [opt.exp_name '_OutputTriggers_' save_time ];

full_file_save_path = [opt.output_path, filesep file_save_name '.mat' ];

target_idx_fd = opt.target_idx_fd;
trigger_idx_fd = opt.trigger_idx_fd;
output = struct();
output.target_ensembles = [];
% indices (note that ROIlist in pyRTAOI only contains accepted ROIs)
if iscell(target_idx_fd)&&numel(target_idx_fd)>=1 % match stim types with target ensembles if more than one target-idx field is provided
     output.target_ensembles = cellfun(@(x)cell_idx_struct.(x),target_idx_fd,'UniformOutput',false);
     output.target_idx = unique(cell2mat(output.target_ensembles));
else
    output.target_idx = cell_idx_struct.(target_idx_fd);
end
output.trigger_idx = cell_idx_struct.(trigger_idx_fd);

% sensory auc
try
output.target_sensory_auc = extractfield(cell_struct(output.target_idx),'correct_stimAUC');
output.target_sensory_auc_zscore = extractfield(cell_struct(output.target_idx),'correct_stimAUC_zscore');
catch
    warning('non auc found')
end
num_triggers = numel(output.trigger_idx);

% dummy values (for use later)
output.trigger_weights = ones(1,num_triggers)./num_triggers;
output.trigger_thresh = 0.5;
output.trigger_frames = 45:60;
output.thresh_sd =[];
if ~isempty(pop_params)
    output.trigger_weights = pop_params.weights;
    output.trigger_thresh = pop_params.thresh;
    output.trigger_frames = pop_params.frames_enable_trigger;
    output.condition_type = pop_params.condition_type;
%     output.pop_opt = opt.pop_opt;
    output.thresh_sd = pop_params.thresh_sd;
end

output.exp_name = opt.exp_name;
output.output_data_path = full_file_save_path;

try 
    output.fds_of_interest = pop_params.fds_of_interest;   
end

try
target_centroids = cell2mat({cell_struct(output.target_idx).centroid}');
trigger_centroids = cell2mat({cell_struct(output.trigger_idx).centroid}');

target_img = make_centroid_image(target_centroids,opt.fov_size,opt.ds_factor);
trigger_img = make_centroid_image(trigger_centroids,opt.fov_size,opt.ds_factor);

imwrite(target_img,[opt.output_path, filesep target_img_save_name '.tif']);
imwrite(trigger_img,[opt.output_path, filesep trigger_img_save_name '.tif']);
catch e
    disp(['Error saving centroid images:',e.message]);
end
save(full_file_save_path, 'output');

disp(['Output saved to:' full_file_save_path])

end

