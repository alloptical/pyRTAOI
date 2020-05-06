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
color.stim3 = [.5 .5 .5]; % catch 
color.stim4 = [.5 .5 .5]; % catch
color.stim5 = [.5 .5 .5]; % catch
color.stim_5_miss = [.5 .5 .5]; % catch 
color.stim_5_lick = [0 0 0]; % catch
color.stim_5_nolick = color.stim_5_miss ; % catch


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
color.lick = [0 0 0];
color.nolick = [.5 .5 .5];
color.hold = [.5 .5 .5];
color.stim1_trials = color.stim1;
color.stim2_trials = color.stim2;
color.hit = color.st_correct_stim_2;
color.cr = color.st_correct_stim_1;
color.miss = color.st_incorrect_stim_2;
color.fa = color.st_incorrect_stim_1;
color.all = [.5 .5 .5];
color.touch = [0 0 0];
color.tex1 = color.cr;
color.tex2 = color.hit;

% performance plot
color.all = [0 0 0];
color.photo = color.photostim;
color.nonphoto = [.5 .5 .5];
color.dummyphoto = tint(color.photo,0.5);
color.photostim = color.photostim;
color.nonphotostim = [.5 .5 .5];
color.dummyphotostim = tint(color.photo,0.5);
for s = 1:4
    color.(['stim_' num2str(s) '_photo']) = color.photo;
    color.(['stim_' num2str(s) '_nonphoto']) = color.nonphoto;
    color.(['stim_' num2str(s) '_dummyphoto']) = color.dummyphoto;  
end

% added for posthoc
photo_types = {'_nonphoto','_photo','_dummyphoto'};
catch_photo_types = {'_nonphoto','_photo_1','_photo_2'};

for s = 1:2
    for p = 1:numel(photo_types)
        this_color = color.(strrep(photo_types{p},'_',''));
        color.(['stim_' num2str(s) photo_types{p}]) = this_color;
        color.(['stim_' num2str(s) photo_types{p} '_correct']) = this_color;
        color.(['stim_' num2str(s) photo_types{p} '_incorrect']) = shade(this_color,.5);
    end
        color.(['stim_' num2str(s) '_correct_ctr']) = color.(['correct_stim'  num2str(s)]);
        color.(['stim_' num2str(s) '_incorrect_ctr']) = color.(['incorrect_stim'  num2str(s)]);

end

s = 5;% catch type
for p = 1:numel(catch_photo_types)
    this_color = color.(strrep(strrep(catch_photo_types{p},'_photo','stim'),'_',''));
    color.(['stim_' num2str(s) catch_photo_types{p}]) = this_color;
    color.(['stim_' num2str(s) catch_photo_types{p} '_lick']) = this_color;
    color.(['stim_' num2str(s) catch_photo_types{p} '_nolick']) = tint(this_color,.5);
end

color.stim_1_photo = color.photo;
color.stim_1_nonphoto = color.nonphoto;
color.stim_1_dummyphoto = color.dummyphoto;

color.stim_2_photo = color.photo;
color.stim_2_nonphoto = color.nonphoto;
color.stim_2_dummyphoto = color.dummyphoto;
end

