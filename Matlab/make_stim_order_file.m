function [stimOrder] = make_stim_order_file()
% find trials with photostim
photo_offset = 120; % frames
cue_offset = 150;

% load data
[file,path] = uigetfile('*.paq','Select paq data');
disp(['Loaded file :',fullfile(path,file)])
paq_path = fullfile(path,file);
stim_frames = get_stim_frames(paq_path,'ToSpiral');
go_cue_frames = get_stim_frames(paq_path,'Cue');
trialon_frames = get_stim_frames(paq_path,'TrialTrigger');
num_trials = numel(trialon_frames);

stimOrder = ones(1,num_trials);
for i = 1:num_trials
    
    if find(go_cue_frames>(trialon_frames(i))& go_cue_frames<(trialon_frames(i)+cue_offset))
        stimOrder(i) = 2; % trials with cue
        if find(stim_frames>(trialon_frames(i))& stim_frames<(trialon_frames(i)+photo_offset))
            stimOrder(i) = 3; % trials with photostim and cue
        end
    end
    
end
temp_name = strsplit(file,'.paq');
save_name = [path filesep temp_name{1} '_stimOrder'];
oris = stimOrder;
save(save_name,'oris')

end

