% PackIO_OutputArrayBuilder
% Builds a voltage output array of repeating triggers for PackIO
% Lloyd Russell May 2014

% variables to customise
filename       = 'LLOYD_VOLTAGE_OUTPUT_AO_NUMBERCELLS_20REPS_FAST';             % filename for saved output file
exp_dur_secs   = 10*20;                                                    % experiment duration in seconds
stim_rate_hz   = 0.05;                                                     % stim rate in hz
trig_dur_msecs = 50;                                                       % trigger duration in milliseconds
sample_rate_hz = 1000;                                                     % output device sample rate in hz
jitter_msecs   = 200;                                                        % jitter, randomise stim intervals. use 0 for no jitter. TOTAL jitter, plus/minus half this amount.

% conversions
exp_length     = exp_dur_secs * sample_rate_hz;
num_trigs      = ceil(stim_rate_hz * exp_dur_secs);                        % calculate number of stim in entire experiment
trig_length    = trig_dur_msecs/1000 * sample_rate_hz;                     % calculate length of trigger pulse at device sample rate
jitter_length  = jitter_msecs/1000 * sample_rate_hz;                       % calculate length of possible jitter at device sample rate

% calculate trigger times/intervals
post_trig_intervals = (floor(1/stim_rate_hz*sample_rate_hz))-(jitter_length/2) + floor(rand(num_trigs-1, 1)*jitter_length); % takes the minimum stim time, adds a random amount of jitter do give ranged with mean centred at dfined stim interval
trig_times          = [0; cumsum(post_trig_intervals)]+1;                  % first stim at time 0, cumulative sum to get absolute times of next stims given by the intervals. plus one to allow correct indexing

% build array
out_array = zeros(exp_length, 1);                                          % initialise the output array
for i = 1:num_trigs
    index = (trig_times(i):((trig_times(i)+trig_length)-1));               % build an 'index array', trigger times to set high
    out_array(index) = 5;
end

% save output
fid = fopen([filename '.dat'],'w','l');
fwrite(fid,out_array,'double');
fclose(fid);

% visualise
figure; plot(out_array);
