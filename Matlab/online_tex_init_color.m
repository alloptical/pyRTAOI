function [color] = online_tex_init_color()
% {'st_correct_stim1','st_incorrect_stim1','st_correct_stim2','st_incorrect_stim2'}
linecolor = [255,153,187;237,23,94;153,186,255;20,81,204]./255;
trace_fds = {'correct_stim1','incorrect_stim1','correct_stim2','incorrect_stim2'};
init_fds = {'st_correct_stim_1','st_incorrect_stim_1','st_correct_stim_2','st_incorrect_stim_2'};



for f = 1:numel(init_fds)
    color.(init_fds{f}) = linecolor(f,:);
end

for f = 1:numel(trace_fds)
    color.(trace_fds{f}) = linecolor(f,:);
end

color.correct = [0 0 0];
color.incorrect = [.5 .5 .5];
color.stim1 = mean([color.correct_stim1;color.incorrect_stim1],1);
color.stim2 = mean([color.correct_stim2;color.incorrect_stim2],1);
color.port1 = mean([color.correct_stim1;color.incorrect_stim2],1);
color.port2 = mean([color.correct_stim2;color.incorrect_stim1],1);

color.correct_trial = tint([0,1,0],.5);
color.incorrect_trial = tint([1,0,0],.5);
color.miss_trial = tint([0,0,0],.5);
color.fa_trial = [0 0 0];
color.photo = [255 60 10]./255;

end

