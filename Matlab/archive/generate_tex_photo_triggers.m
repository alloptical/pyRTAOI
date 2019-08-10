function [output_triggers] = generate_tex_photo_triggers(save_path,all_trials_num_target,opt)
num_trials = getOr(opt,'tot_num_trials',100);
pc_stim_trials = getOr(opt,'pc_stim_trials',100);        % percentage of trials with photostimulation
trial_interval = getOr(opt,'trial_interval',10);        % interval between trialon triggers, i.e., trial length (sec)
stim_start_time =  getOr(opt,'stim_start_time',3);      % time delay between trialon trigger and photostim trigger (sec)
slm_switch_time =  getOr(opt,'slm_switch_time',2);      % time delay between trialon trigger and slm trigger(sec, at least 5 ms before photostim)
sample_rate_hz = getOr(opt,'sample_rate_hz',1000);   
trig_dur_msecs = getOr(opt,'trig_dur_msecs',5);         %impulse duration in milliseconds
trig_amp_volts = getOr(opt,'trig_amp_volts',5);
save_name_ext =  getOr(opt,'save_name_ext','');
trig_length = round(trig_dur_msecs/1000*sample_rate_hz);
exp_length = (num_trials+2)*trial_interval*sample_rate_hz;
trig_times = (1:num_trials).*trial_interval.*sample_rate_hz;

trialon_array = zeros(exp_length, 1);   % impulse array (starts of wavelets)
stimon_array = zeros(exp_length, 1);
slm_array = zeros(exp_length, 1);
aom_array = zeros(exp_length, 1);
aom_single_pulse_waveform = zeros(1,round(1/opt.photostim_frequency*sample_rate_hz));
aom_single_pulse_waveform(1:round(opt.duration_per_stim/1000*sample_rate_hz)) = 1;
aom_pulse_waveform = repmat(aom_single_pulse_waveform,1,opt.num_stim_per_trial);

% for safety
aom_pulse_waveform =[zeros(1,10) aom_pulse_waveform];
aom_pulse_waveform(end-1:end) = 0;

% omit photostim triggers randomly
if pc_stim_trials<100
    rd = randperm(num_trials);
    num_stim_trials = round(num_trials*pc_stim_trials/100);
    stim_trial_idx = rd(1:num_stim_trials);

else
    num_stim_trials = num_trials;
    stim_trial_idx = 1:num_trials;
end

slm_first_trial = true;
for i = 1:num_trials
    index = (trig_times(i):((trig_times(i)+trig_length)-1));
    trialon_array(index) = trig_amp_volts;
    this_num_targets = all_trials_num_target(i);
    this_aom_volt = mw2pv(opt.mw_per_cell *this_num_targets)/1000; 

    if ismember(i,stim_trial_idx)
        stimon_array(index + sample_rate_hz*stim_start_time) = trig_amp_volts;
        aom_array(trig_times(i)+[1:length(aom_pulse_waveform)] + sample_rate_hz*stim_start_time) = aom_pulse_waveform.*this_aom_volt;
        if slm_first_trial
            slm_first_trial = false; % omit the first trigger to SLM
        else
            slm_array(index + sample_rate_hz*slm_switch_time) = trig_amp_volts;
        end
    end
end

figure('name','triggers generated'); hold on;
plot(slm_array,'b'); 
plot(stimon_array,'r'); 
plot(trialon_array,'color','black')
plot(aom_array,'c')
legend('toSLM','toSpiral','TrialOn')

trialon_save_name = ['TrialOn_volt_output' num2str(num_trials) '_x_pulses_isi_' num2str(trial_interval) '_' save_name_ext '.dat'];
stimon_save_name = ['toSpiral_volt_output' num2str(num_trials) '_x_pulses_isi_' num2str(trial_interval)  '_' save_name_ext '.dat'];
slm_save_name = ['toSLM_volt_output' num2str(num_trials) '_x_pulses_isi_' num2str(trial_interval)   '_' save_name_ext '.dat'];
aom_save_name = ['toAOM_volt_output' num2str(num_trials) '_x_pulses_isi_' num2str(trial_interval)   '_' save_name_ext '.dat'];

% save to structure
output_triggers = struct();
output_triggers.TrialOn = trialon_array;
output_triggers.toSpiral = stimon_array;
output_triggers.toSLM = slm_array;
output_triggers.toAOM = aom_array;


% save output files

if ~exist(save_path, 'dir')
mkdir(save_path)
disp(['trigger generateor has made a new save path:' save_path])
end

fid = fopen([save_path filesep trialon_save_name],'w','l');
fwrite(fid,trialon_array,'double');
fclose(fid);

fid = fopen([save_path filesep stimon_save_name] ,'w','l');
fwrite(fid,stimon_array,'double');
fclose(fid);

fid = fopen([save_path filesep slm_save_name] ,'w','l');
fwrite(fid,slm_array,'double');
fclose(fid);

fid = fopen([save_path filesep aom_save_name] ,'w','l');
fwrite(fid,aom_array,'double');
fclose(fid);

disp(['Trigger files saved to ' save_path])
end

