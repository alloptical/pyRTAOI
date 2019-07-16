function [cnm_struct,cnm_image,num_comp,cnm_dims,tot_frames] = make_cnm_struct(caiman_data, varargin)
%% make cnm data structure out of pyrtaoi results 2019
% adapted from OnlineProcVisual
stim_frames = [];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'stim_frames')
        stim_frames = varargin{v+1};
    end
end

cnm_struct = struct();
cnm_dims = caiman_data.cnm_dims;
cnm_image = reshape(caiman_data.cnm_b,cnm_dims);
cnm_A = full(caiman_data.cnm_A);
num_frames = caiman_data.t_cnm;


num_comp = size(cnm_A,2);
cm = com(sparse(double(cnm_A)),cnm_dims(1),cnm_dims(2)); % center of mass

% deal with skipped frames
if ~isempty(caiman_data.frames_skipped)
    skip_frames = caiman_data.frames_skipped + caiman_data.num_frames_init;
    tot_frames = num_frames + numel(skip_frames);
    caiman_frames = setdiff([1:tot_frames],skip_frames);
else
    caiman_frames = 1:num_frames;
    tot_frames = num_frames;
end

for i = 1:num_comp
    cnm_struct(i).shape = reshape(cnm_A(:,i),cnm_dims);
    cnm_struct(i).noisyC = caiman_data.noisyC(i+1,1:num_frames);
    cnm_struct(i).deconvC = caiman_data.cnm_C(i,1:num_frames);
    cnm_struct(i).centroid = cm(i,:);
    cnm_struct(i).frame_added = find(cnm_struct(i).noisyC >0,1);
    
    % set skipped frames to nan
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.cnm_C(i,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).deconvC_full = temp_trace;
    
    temp_trace = nan(1,tot_frames);
    temp_trace(caiman_frames) =  caiman_data.noisyC(i+1,1:num_frames);
    temp_trace = fillmissing(temp_trace,'linear');
    cnm_struct(i).noisyC_full = temp_trace;
    
    % caiman stim frames
    if(~isempty(stim_frames))
        cnm_struct(i).stim_frames = sens_stim_frames(stim_frames>cnm_struct(i).frame_added);
    else
        cnm_struct(i).stim_frames = [];
    end
end


end

