function [output_trials,odd_trial_idx] = make_trials_struct(behavior_data,varargin)
withold_win_start = 3; 
reward_win_end = 5;
tot_num_trials = [];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'tot_num_trials')
        tot_num_trials = varargin{v+1};
    end
    if strcmpi(varargin{v},'withold_win_start')
        withold_win_start = varargin{v+1};
    end
end
% make trials.m struct from pybehavior output
fd_names = {'stim_type','correct','incorrect','miss','firstresponse','cheated','stim_var','fa'};
num_fd = numel(fd_names);
odd_trial_idx = [];
output_trials = struct();
if isempty(tot_num_trials)
    tot_num_trials = numel(behavior_data.results);
end
this_norm_trials = 1:tot_num_trials;
func = @(x)func_fristresponse(x,withold_win_start,reward_win_end);
for i = 1:num_fd
    this_fd = fd_names{i};
    if isfield(behavior_data.results{2},this_fd)
        this_data = cell(1,tot_num_trials);
        
        if strcmp(this_fd,'firstresponse') 
            % get first response with actual response recorded AFTER WITHOLD WINDOW ONSET (fa trials first response were all zero!)
            this_responses = cellfun(@(x)x.('responses'),behavior_data.results,'UniformOutput',false);
            this_first_response = nan(1,tot_num_trials);
            this_first_lick = nan(1,tot_num_trials);
            [this_first_lick(~cellfun(@(x)isempty(x),this_responses)),this_first_response(~cellfun(@(x)isempty(x),this_responses))] = cellfun(@(x)func(x), this_responses(cellfun(@(x)~isempty(x),this_responses)));
            
            output_trials.(this_fd) = this_first_response(1:tot_num_trials);
            output_trials.('firstlick') = this_first_lick(1:tot_num_trials);
      
            continue;
        else
            this_data(this_norm_trials) = cellfun(@(x)extractfield(x,this_fd),behavior_data.results(this_norm_trials),'UniformOutput',false);
            this_data(odd_trial_idx) = {nan};           
        end
%         try
%         this_data = cellfun(@(x)extractfield(x,this_fd),behavior_data.results,'UniformOutput',false);
%         catch
%                     this_data = cellfun(@(x)extractfield(x,this_fd),behavior_data.results(1:tot_num_trials-1),'UniformOutput',false);
%         end
        if iscell(this_data)
            checktype = cellfun(@(x)class(x),this_data,'UniformOutput',false);
            [types,~,whichcell]=unique(checktype);
            if contains(types,'cell')
                this_data(this_norm_trials) = cellfun(@(x)cell2mat(x),this_data(this_norm_trials),'UniformOutput',false);
                
            end
            [this_type_idx]= mode(whichcell);

            this_odd_trial = find(whichcell~= this_type_idx);
            odd_trial_idx =unique( [odd_trial_idx,this_odd_trial]);
            this_norm_trials = setdiff(1:tot_num_trials,odd_trial_idx);
            
            this_data_mat = zeros(1,tot_num_trials);
            this_data_mat(this_norm_trials) = cell2mat(this_data(this_norm_trials));
            this_data_mat(this_odd_trial) = nan;

            if numel(types)>1
                warning(['some trial is odd! trial ' num2str(this_odd_trial') ' was excluded in ' this_fd])
            end
            
        end
        output_trials.(this_fd) = this_data_mat;
    else
        warning(['field not found in behavior: ' this_fd ', skipped'])
    end
end



% exclude odd trials
odd_trial_idx = unique(odd_trial_idx);
normal_trial_idx = setdiff(1:tot_num_trials,odd_trial_idx);
if ~isempty(odd_trial_idx)
    disp(['odd trials indices: ' num2str(odd_trial_idx') ]) 
%     output_trials = structfun(@(x)x(normal_trial_idx),output_trials,'UniformOutput',false);
end
%%
% figure
% test_mat = [output_trials.correct; output_trials.incorrect; output_trials.fa; output_trials.miss];
% imagesc(test_mat)
% unique(sum(test_mat,1))
%%
end
function [firstlicktime,firstresponseport] = func_fristresponse(x,withold_win_start,reward_win_end)
firstlicktime = x(1,find(x(1,:)>withold_win_start&x(1,:)<reward_win_end,1)); % first response after withold window
firstresponseport = x(2,find(x(1,:)>withold_win_start&x(1,:)<reward_win_end,1)); % first response after withold window

if isempty(firstlicktime)
    firstlicktime = nan;
    firstresponseport = nan;
end

end
