function [frames_with_photo_stim,frames_with_sens_stim,run_speed,aom_volt] = read_paq_file(paq_full_path)
    [paq_data, paq_channels,~,paq_rate ]= paq2lab_silent(paq_full_path);

    % init
    run_speed = [];
    aom_volt = [];
    frames_with_photo_stim = [];
    frames_with_sens_stim = [];
    
    % set paq channels
    paq_toPV_idx = [];
    paq_run_idx = [];
    paq_toSensory_idx = []; % rtaoi output
    paq_Sensory_idx = []; % pybehaviour output
    paq_frame_idx = [];
    paq_aom_idx = [];
    
    
    paq_toPV_idx  = find(not(cellfun('isempty',  strfind(paq_channels,'ToSpiral'))));
    paq_frame_idx  = find(not(cellfun('isempty', strfind(paq_channels,'FrameEnd'))));
    paq_toSensory_idx = find(not(cellfun('isempty',  strfind(paq_channels,'TrialTrigger'))));%
    paq_aom_idx = find(not(cellfun('isempty',  strfind(paq_channels,'AOM'))));%

    paq_run_idx = find(not(cellfun('isempty', strfind(paq_channels,'Running'))));
    sample_rate = paq_rate;
    
    % frame times
    frame_times_sample = pt_continuousabove(paq_data(:,paq_frame_idx),0,1,0,100000,0);
    frame_times_sample = frame_times_sample(:,1);
    
    % get runnning speed around frame times
    if(~isempty(paq_run_idx))
        try
            samples_per_frame = round(mean(diff(frame_times_sample)));
            num_frames = length(frame_times_sample);
            run_speed = zeros(1,num_frames);
            for i = 1:num_frames
                run_speed(i) = mean(paq_data(frame_times_sample(i)-samples_per_frame:frame_times_sample(i),paq_run_idx));
            end
        catch
            warning('no running speed found')
            
        end
    end
    
    % get aom control voltage around frame times
    if(~isempty(paq_aom_idx))
        try
            samples_per_frame = round(mean(diff(frame_times_sample)));
            num_frames = length(frame_times_sample);
            aom_volt = zeros(1,num_frames);
            for i = 1:num_frames
                aom_volt(i) = max(paq_data(frame_times_sample(i)-samples_per_frame:frame_times_sample(i),paq_aom_idx));
            end
        catch
            warning('aom volt not found')
            
        end
    end
    
    % photo stim times
    if(~isempty(paq_toPV_idx))
        try
            stim_times_sample = pt_continuousabove(paq_data(:,paq_toPV_idx),0,0.1,0,100000,1000);
            stim_times_sample = stim_times_sample(:,1);
            
            % get delay between frame-end and toPV trigger (spirals start immediately upon toPV trigger!)
            frame_before_sec = [];
            frame_to_stim_msec = [];
            for i  =1:length(stim_times_sample)
                frame_before_sec(i) =  frame_times_sample(find(frame_times_sample<stim_times_sample(i), 1 , 'last'));
                frame_to_stim_msec(i) =  stim_times_sample(i)-frame_before_sec(i);
            end
            frame_to_stim_msec = frame_to_stim_msec*1000/sample_rate;
            
            % get stim frames
            frame_times_sec = frame_times_sample./ sample_rate;
            stim_times_sec = stim_times_sample./ sample_rate;
            for i  =1:length(stim_times_sample)
                frames_with_photo_stim(i) = find(frame_times_sec>stim_times_sec(i), 1 , 'first');
            end
        catch
            warning('no photostim found')
        end
        
    end
    
    % sensory stim times
    if(~isempty(paq_toSensory_idx))
        try
        sens_stim_times_sample = pt_continuousabove(paq_data(:,paq_toSensory_idx),0,1,0,100000,1000);

        
        sens_stim_times_sample = sens_stim_times_sample(:,1);

        
        % get sensory stim frames
        frame_times_sec = frame_times_sample./ sample_rate;
        sens_stim_times_sample = sens_stim_times_sample./ sample_rate;
        for i  =1:length(sens_stim_times_sample)
            frames_with_sens_stim(i) = find(frame_times_sec>sens_stim_times_sample(i), 1 , 'first');
        end
        catch
            warning('no sensory stim found')
        end
    end
   % end of reading paqs

end

