function [trial_seq] = generate_texture_trial_seq(loop_stim_types,loop_stim_vars,num_loops,save_dir,fname)
% genenrate trial order file for pybehavior to load
% adapted from experiment_trial_seq.m from Oli

% stim 1 - texture 1 trials 
% stim 2 - texture 2 trials
% stim 3 - catch trials, reward 1
% stim 4 - catch trials, reward 2
% stim 5 - catch trials, no reward



stims = cell(1,num_loops);
vars = cell(1,num_loops);

for loop = 1:num_loops

    if loop == 1
        stims{loop} = [1,2,1,2,1,2,1,2,1,2];%2,1,2,1,2,1,2,1,2]; % first set of trials are 100% easy trials in interleaved order
        vars{loop} =  [1 1 1 1 1 1 1 1 1 1];%,6,3,6,3,6,3,6,3,6];      
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



dlmwrite([save_dir filesep fname '.txt'],trial_seq)
disp(['tot. number of trials:' num2str(size(trial_seq,2))])
disp(['saved as:' save_dir filesep fname '.txt'] )
end

