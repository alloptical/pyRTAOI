function [pybehav_trial_seq,rtaoi_trial_seq] = generate_texture_trial_seq(loop_pybehavstim_types,loop_pybehav_stim_vars,loop_rtaoi_stim_types,loop_rtaoi_stim_vars,num_loops,save_dir,fname)
% genenrate trial order file for pybehavior to load
% adapted from experiment_trial_seq.m from Oli

% stim 1 - texture 1 trials 
% stim 2 - texture 2 trials
% stim 3 - catch trials, reward 1
% stim 4 - catch trials, reward 2
% stim 5 - catch trials, no reward



stims = cell(1,num_loops);
vars = cell(1,num_loops);
rtaoi_stims = cell(1,num_loops);
rtaoi_vars = cell(1,num_loops);

for loop = 1:num_loops

    if loop == 1
        stims{loop} = [1,2,1,2,1,2,1,2,1,2];%2,1,2,1,2,1,2,1,2]; % first set of trials are 100% easy trials in interleaved order
        vars{loop} =  [1 1 1 1 1 1 1 1 1 1];%,6,3,6,3,6,3,6,3,6];   
        rtaoi_stims{loop} = stims{loop};
        rtaoi_vars{loop} = vars{loop};

    else
        % baseline imaging
       stims{loop} =  loop_pybehavstim_types;
       vars{loop} =   loop_pybehav_stim_vars;  
       
       rtaoi_stims{loop} =  loop_rtaoi_stim_types;
       rtaoi_vars{loop} =   loop_rtaoi_stim_vars;  
     
               
        while true
            rand_ind = randperm(length(stims{loop}));
            stimtemp = stims{loop}(rand_ind);
            rtaoi_stimtemp =  rtaoi_stims{loop}(rand_ind);
            if all(abs((diff(diff(stimtemp))))>0) % diff diff checks for > 3 in a row
                vars{loop} = vars{loop}(rand_ind);
                stims{loop} = stimtemp;
                
                rtaoi_vars{loop} =  rtaoi_vars{loop}(rand_ind);
                rtaoi_stims{loop} =   rtaoi_stimtemp;
              
                
                break
            end                
        end    
    end
end

pybehav_trial_seq(1,:) = horzcat(stims{:});
pybehav_trial_seq(2,:) = horzcat(vars{:});
rtaoi_trial_seq(1,:) = horzcat(rtaoi_stims{:});
rtaoi_trial_seq(2,:) = horzcat(rtaoi_vars{:});



dlmwrite([save_dir filesep fname '.txt'],pybehav_trial_seq)
disp(['tot. number of trials:' num2str(size(pybehav_trial_seq,2))])
disp(['saved as:' save_dir filesep fname '.txt'] )
end

