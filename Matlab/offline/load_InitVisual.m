%% file opsin positive sensory cells in an initialisation file
% load data
caiman_data = load('Y:\zzhang\Data\20181024\pyrtaoi_results\20181024_L527_t_0011_rtaoi_DS_2.0_OnlineProc_DS_2.0_165202.mat') % - example

%% color lut
hsv = colormap(hsv);
hsv = hsv(2:end-3,:);
close 

%% stim types
num_stim_type = 8; % orientations
num_stim_per_type = 4;
vis_duration = 30; % frames
photo_duration = 90; % frames
vis_stim_frames = {};
init_vis_stim_frames = {};

if(~isempty(caiman_data.stim_frames_caiman))
    sens_stim_frames = caiman_data.sensory_stim_frames+caiman_data.num_frames_init;
else
    sens_stim_frames = [];
end
for i = 1:num_stim_type
    vis_stim_frames{i} = sens_stim_frames(i-1+[1:num_stim_type:num_stim_type*num_stim_per_type]);
    init_vis_stim_frames{i} = vis_stim_frames{i} - caiman_data.num_frames_init;
end
