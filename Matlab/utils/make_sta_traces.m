function sta_traces = make_sta_traces(input_traces, stim_frames, pre_frames, post_frames)

sta_traces = zeros(size(input_traces, 1), length(stim_frames), pre_frames+post_frames+1);

for roi = 1:size(input_traces, 1)
    for stim = 1:length(stim_frames)
        try
            sta_traces(roi, stim, :) = input_traces(roi, stim_frames(stim)-pre_frames:stim_frames(stim)+post_frames);
        catch
            sta_traces(roi, stim, :) = NaN;
        end
    end
end