function [color] = online_tex_init_color()
% {'st_correct_stim1','st_incorrect_stim1','st_correct_stim2','st_incorrect_stim2'}
linecolor = [255,153,187;237,23,94;153,186,255;20,81,204; 255 51 255; 153 0 153;51 255 153;0 153 76]./255;
trace_fds = {'correct_stim1','incorrect_stim1','correct_stim2','incorrect_stim2',...
    'correct_stim3','incorrect_stim3','correct_stim4','incorrect_stim4'};
mod_trace_fds = {'stim_1_correct','stim_1_incorrect','stim_2_correct','stim_2_incorrect'};
smooth_trace_fds = {'st_correct_smooth_deconv_1','st_incorrect_smooth_deconv_1','st_correct_smooth_deconv_2','st_incorrect_smooth_deconv_2'};

init_fds = {'st_correct_stim_1','st_incorrect_stim_1','st_correct_stim_2','st_incorrect_stim_2',...
    'st_correct_stim_3','st_incorrect_stim_3','st_correct_stim_4','st_incorrect_stim_4'};

for f = 1:numel(init_fds)
    color.(init_fds{f}) = linecolor(f,:);
end

for f = 1:numel(trace_fds)
    color.(trace_fds{f}) = linecolor(f,:);
end
for f = 1:numel(mod_trace_fds)
    color.(mod_trace_fds{f}) = linecolor(f,:);
end

for f = 1:numel(smooth_trace_fds)
    color.(smooth_trace_fds{f}) = linecolor(f,:);
end
color.correct = [0 0 0];
color.incorrect = [.5 .5 .5];
color.stim1 = mean([color.correct_stim1;color.incorrect_stim1],1);
color.stim2 = mean([color.correct_stim2;color.incorrect_stim2],1);
color.stim3 = [.5 .5 .5]; % catch stimulus
color.stim4 = [.5 .5 .5]; % catch
color.stim5 = [.5 .5 .5]; % catch
color.stim_5_miss = [.5 .5 .5]; % catch
color.stim_5_lick = [0 0 0]; % catch


color.port1 = mean([color.correct_stim1;color.incorrect_stim2],1);
color.port2 = mean([color.correct_stim2;color.incorrect_stim1],1);
color.stim5 = [.5 .5 .5];
color.stim6 = [0 0 0];
color.stim3 = [0.5 0.5 0.5];

color.correct_trial = tint([0,1,0],.5);
color.incorrect_trial = tint([1,0,0],.5);
color.miss_trial = tint([0,0,0],.5);
color.fa_trial = [0 0 0];
color.photo = [255 60 10]./255;
color.photostim = [255 72 36]./255;

color.L1 = [0 0 0];
color.L2 = [.5 .5 .5];

color.trigger = [103, 232, 235]./255;
color.target = color.photo;

color.go_trials = [0 0 0];
color.nogo_trials = [.5 .5 .5];

% performance plot
color.all = [0 0 0];
color.photo = color.photostim;
color.nonphoto = [.5 .5 .5];
color.dummyphoto = tint(color.photo,0.5);

end

