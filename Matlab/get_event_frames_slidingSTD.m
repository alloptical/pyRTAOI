function [ event_frames ] = get_event_frames_slidingSTD( trace,N,all_reset_frames,window_size,pre_exp_frames)
% mimicing simple online event detection

event_frames    = [];
thresh      = [];
idx = 1;
last_trigger = 0;


for i = pre_exp_frames+window_size+1:length(trace)
    thresh(i) = mean(trace(i-window_size:i))+N*std(trace(i-window_size:i));
    
    if (trace(i)>thresh(i)) && (i>last_trigger+all_reset_frames)
        event_frames(idx) = i;
        idx = idx+1;
        last_trigger = i;
    end
    
end

end

