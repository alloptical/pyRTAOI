function [pybehav_trial_seq] = generate_BASELINE_texture_trial_seq(varargin)
% genenrate trial order file for pybehavior to load
% adapted from experiment_trial_seq.m from Oli

% stim 1 - texture 1 trials
% stim 2 - texture 2 trials
% stim 3 - catch texture, reward 1 (or reward 0 for gonogo)
% stim 4 - catch texture, reward 2
% stim 5 - catch texture, no reward


save_dir = [];
fname = [];
loop_pybehavstim_types = [repmat([1,2],[1,9]),3,4]; %10% catch trials
loop_pybehav_stim_vars = ones(size(loop_pybehavstim_types)); % dummy. not using this feature in pybehav
num_loops = 5;

for v = 1:numel(varargin)
    if strcmpi(varargin{v},'save_dir')
        save_dir = varargin{v+1};
    end
    if strcmpi(varargin{v},'fname')
        fname = varargin{v+1};
    end
    if strcmpi(varargin{v},'loop_pybehavstim_types')
        loop_pybehavstim_types = varargin{v+1};
    end
    if strcmpi(varargin{v},'loop_pybehav_stim_vars')
        loop_pybehav_stim_vars = varargin{v+1};
    end
    if strcmpi(varargin{v},'num_loops')
        num_loops = varargin{v+1};
    end
end


stims = cell(1,num_loops);
vars = cell(1,num_loops);

for loop = 1:num_loops
    
    if loop == 1
        stims{loop} = [1,2,1,2,1,2,1,2,1,2];% first set of trials are 100% easy trials in interleaved order
        vars{loop} =  [1 1 1 1 1 1 1 1 1 1];
        
    else
        % baseline imaging
        stims{loop} =  loop_pybehavstim_types;
        vars{loop} =   loop_pybehav_stim_vars;
        
        
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

pybehav_trial_seq(1,:) = horzcat(stims{:});
pybehav_trial_seq(2,:) = horzcat(vars{:});

if isempty(fname)
    save_time = datestr(now,'yyyymmdd_HHMM');
    x = inputdlg({'ID'}, 'Enter Mouse ID',[1 50],{'cb'});
    fname = [x{1} '_baseline_PyBehavior_' save_time ];
end

if isempty(save_dir)
    disp('Select directory to save trial sequence')
    pause(1)
    save_dir = uigetdir();
end
dlmwrite([save_dir filesep fname '.txt'],pybehav_trial_seq)
disp(['tot. number of trials:' num2str(size(pybehav_trial_seq,2))])
disp(['saved as:' save_dir filesep fname '.txt'] )
end

