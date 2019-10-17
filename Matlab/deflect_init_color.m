function [color] = deflect_init_color()
% {'st_correct_stim1','st_incorrect_stim1','st_correct_stim2','st_incorrect_stim2'}
linecolor = [255,153,187;237,23,94;153,186,255;20,81,204]./255;
trace_fds = {'TS_L1','ipsi_L1','TS_L2','contra_L2'};
color.gocue = [255 153 51]./255;
color.correct = [0 0 0];
color.incorrect = [.5 .5 .5];

color.correct_trial = tint([0,1,0],.5);
color.incorrect_trial = tint([1,0,0],.5);
color.miss_trial = tint([0,0,0],.5);
color.fa_trial = [0 0 0];
color.photo = [255 60 10]./255;


for f = 1:numel(trace_fds)
    color.(trace_fds{f}) = linecolor(f,:);
end
color.ipsi_L2 = color.contra_L2;
color.contra_L1 = color.ipsi_L1;
color.L1 = (color.contra_L1+color.TS_L1).*.5;
color.L2 = (color.contra_L2+color.TS_L2).*.5;

% stim types defined by deflection_experiment_trial_seq
% stim 1 - left discrimination trials (10 conditions)
% stim 2 - right discrimination trials (10 conditions)
% stim 3 - threshold trials (+ nogo catch)
% stim 4 - sensory + photostim (10 conditions)
% stim 5 - motor only (no cue)
% stim 6 - cue only (no motor)

% get colors lut for trial types and variation
color.stim1 = color.L1;
color.stim2 = color.L2;
color.stim3 = (color.L1+color.L2).*0.5;
color.stim4 = color.photo;
color.stim5 = [.3 .3 .3];
color.stim6 = color.gocue;








end

% %%
% figure; hold on
% for c = 1:size(linecolor,1)
%     plot([0,1],[c c],'color',linecolor(c,:),'linewidth',2)
% end
% ylim([0 c+1])