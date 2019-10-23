function [trial_seq] = generate_deflection_trial_seq(loop_stim_types,loop_stim_vars,num_loops,save_dir,fname)
% genenrate trial order file for pybehavior to load
% adapted from experiment_trial_seq.m from Oli

% stim 1 - left discrimination trials (10 conditions)
% stim 2 - right discrimination trials (10 conditions)
% stim 3 - threshold trials (+ nogo catch)
% stim 4 - sensory + photostim (10 conditions)
% stim 5 - motor only (no cue) - THESE ARE NO CUE AND NO MOTOR TRIALS
% stim 6 - cue only (no motor)

stims = cell(1,num_loops);
vars = cell(1,num_loops);

for loop = 1:num_loops

    if loop == 1
        stims{loop} = [1,2,1,2,1,2,1,2,1,2];%2,1,2,1,2,1,2,1,2]; % first set of trials are 100% easy trials in interleaved order
        vars{loop} =  [9 9 9 9 9 9 9 9 9 9];%,6,3,6,3,6,3,6,3,6];      
    else
        % baseline imaging
       stims{loop} =  loop_stim_types; % 20191018 - changed this slightly so its just the diagonal psychometric curve
       vars{loop} =   loop_stim_vars;   
       
        
        while true
            rand_ind = randperm(length(stims{loop}));
            stimtemp = stims{loop}(rand_ind);
            if all(abs((diff(diff(stimtemp))))>0) % diff diff checks for > 3 in a row
                vars{loop} = vars{loop}(rand_ind);
                stims{loop} = stimtemp;
                break
            end                
        end    
    end
end

trial_seq(1,:) = horzcat(stims{:});
trial_seq(2,:) = horzcat(vars{:});


% adding easy trials
easyLeftTrial = 10:6:10000;
easyRightTrial = 13:6:10000;

stim_count = 1;

try
    for i = 1:900
        if any(easyLeftTrial== i) 
            final_stim(i) = 1; 
            final_var(i) = 9;%easyLeftvar(c1);c1 = c1+1;
        elseif any(easyRightTrial== i) 
            final_stim(i) = 2; 
            final_var(i) = 9;%easyRightvar(c2);c2=c2+1;
        else
            final_stim(i) = trial_seq(1,stim_count); 
            final_var(i) = trial_seq(2,stim_count); stim_count = stim_count +1;
        end
    end 
catch
end

clear trial_seq

trial_seq(1,:) = final_stim;
trial_seq(2,:) = final_var;

%

dlmwrite([save_dir filesep fname '.txt'],trial_seq)
disp(['tot. number of trials:' num2str(size(trial_seq,2))])
disp(['saved as:' save_dir filesep fname '.txt'] )
end

