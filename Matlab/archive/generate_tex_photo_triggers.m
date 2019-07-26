function [output_triggers] = generate_tex_photo_triggers(save_path,opt)
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

% omit photostim triggers randomly
if pc_stim_trials<100
    rd = randperm(num_trials);
    num_stim_trials = round(num_trials*pc_stim_trials/100);
    stim_trial_idx = rd(1:num_stim_trials);

else
    num_stim_trials = num_trials;
    stim_trial_idx = 1:num_trials;
end

for i = 1:num_trials
    index = (trig_times(i):((trig_times(i)+trig_length)-1));
    trialon_array(index) = trig_amp_volts;
    if ismember(i,stim_trial_idx)
        stimon_array(index + sample_rate_hz*stim_start_time) = trig_amp_volts;
        slm_array(index + sample_rate_hz*slm_switch_time) = trig_amp_volts;
    end
end

figure('name','triggers generated'); hold on;
plot(slm_array,'b'); 
plot(stimon_array,'r'); 
plot(trialon_array,'color','black')
legend('toSLM','toSpiral','TrialOn')

trialon_save_name = ['TrialOn_volt_output' num2str(num_trials) '_x_pulses_isi_' num2str(trial_interval) '_' save_name_ext '.dat'];
stimon_save_name = ['toSpiral_volt_output' num2str(num_trials) '_x_pulses_isi_' num2str(trial_interval)  '_' save_name_ext '.dat'];
slm_save_name = ['toSLM_volt_output' num2str(num_trials) '_x_pulses_isi_' num2str(trial_interval)   '_' save_name_ext '.dat'];

% save to structure
output_triggers = struct();
output_triggers.TrialOn = trialon_array;
output_triggers.toSpiral = stimon_array;
output_triggers.toSLM = slm_array;


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

disp(['Trigger files saved to ' save_path])
end

