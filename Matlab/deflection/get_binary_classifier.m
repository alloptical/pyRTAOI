function [decod_struct] = get_binary_classifier( decod_struct,this_proj_struct, opt,varargin )
% logistic regression for binary classification
% input traces: [trial, time, state]
IF_FRAMEWISE = 0;
IF_CROSSVAL = 0;
IF_MULTI_COMPO = 0;
MODEL_TYPE = 'logistic';
num_fit_shuf = 100;
num_comp = 1;
try
    opt.gocue_frame_adj;
catch
    opt.gocue_frame_adj = opt.gocue_bin;
end
frames_to_avg = getOr(opt,'frames_to_avg',[-10:1:0]+ opt.gocue_frame_adj); % use this frames to train classifier (if IF_FRAMEWISE is false)
frames_to_train =  getOr(opt,'frames_to_train',[-60:1:60] + opt.gocue_frame_adj); % train classifiers only for these frames
fd_names = {'correct','incorrect'};
IF_SHUFFLE = true;
this_struct = this_proj_struct;
fds = fields(this_struct);
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'IF_SHUFFLE')
        IF_SHUFFLE = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'IF_FRAMEWISE')
        IF_FRAMEWISE = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'MODEL_TYPE')
        MODEL_TYPE = varargin{v+1};
    end
    if strcmpi(varargin{v},'IF_CROSSVAL')
        IF_CROSSVAL = varargin{v+1};
    end
    if strcmpi(varargin{v},'field_idx')
        field_idx = varargin{v+1};
    end
    if strcmpi(varargin{v},'fd_names')
        fds = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'num_fit_shuf')
        num_fit_shuf = varargin{v+1};
    end
end



% initiation
first_fd = fds{1};
second_fd = fds{2};
coef_shuf = {};
coef0_shuf = {};
all_accur_shuf = [];
if length(size(this_struct.(first_fd)))>2 % check if data has multi-components
    IF_MULTI_COMPO = 1;
    num_comp = size(this_struct.(first_fd),3);
end

this_correct_fd = first_fd;
this_incorrect_fd = second_fd;


data.correct= [];
data.incorrect = [];

trial_length = opt.trial_length;

if IF_FRAMEWISE
    accur_shuf = nan(trial_length,num_fit_shuf);
    b_shuf = nan(trial_length,num_comp+1,num_fit_shuf);
else
    accur_shuf = nan(1,num_fit_shuf);
    b_shuf = nan(1,num_comp+1,num_fit_shuf);
end

if ~IF_FRAMEWISE
    data.correct = reshape(mean(this_struct.(this_correct_fd)(:,frames_to_avg,:),2),[],num_comp);
    data.incorrect = reshape(mean(this_struct.(this_incorrect_fd)(:,frames_to_avg,:),2),[],num_comp);
    
else
    data.correct = this_struct.(this_correct_fd);
    data.incorrect = this_struct.(this_incorrect_fd);

end


% if (~strcmp(MODEL_TYPE,'elastic'))&&(~IF_FRAMEWISE)&&(~IF_MULTI_COMPO)
%     data.correct = data.correct(:,frames_to_avg); %[trial, frame]
%     data.incorrect =   data.incorrect(:,frames_to_avg);
% end
num_trial.correct = size(this_struct.((this_correct_fd)),1);
num_trial.incorrect = size( this_struct.((this_incorrect_fd)),1);


%% FRAMEWISE CLASSIFIER
thresh_framewise = nan(1,trial_length);
thresh_shuf = nan(1,num_fit_shuf);

if IF_FRAMEWISE % different classifiers across frames
    
    for fr = frames_to_train
        if IF_MULTI_COMPO
            this_data.correct = squeeze(data.correct(:,fr,:));
            this_data.incorrect = squeeze(data.incorrect(:,fr,:));
            
        else
            this_data.correct = data.correct(:,fr);
            this_data.incorrect =  data.incorrect(:,fr);
        end
        
        %% TRAIN CLASSIFIER
        if IF_SHUFFLE % shuffle trials for a better fit; helps to deal with unbalanced sample size
            [min_num_sample,min_fd_idx] = min([num_trial.correct num_trial.incorrect]);
            
            min_fd = fd_names{min_fd_idx};
            shuf_fd_idx = setdiff([1,2],min_fd_idx);
            shuf_fd = fd_names{shuf_fd_idx};
            max_num_sample = num_trial.(shuf_fd);
            
            this_small_data = this_data.(min_fd);
            this_big_data = this_data.(shuf_fd); %
            this_cmp_type = [ones(1,num_trial.(min_fd)),ones(1,num_trial.(min_fd)).*2]';
            
            if num_fit_shuf>10 % only use parfor when may shuffles
                parfor s = 1:num_fit_shuf
                    num_folds = min(5,floor(min_num_sample/2)); % number of folds for cross-validation
                    this_shuf_idx =  randsample(1:max_num_sample,min_num_sample);
                    this_shuf_sample = this_big_data(this_shuf_idx,:);
                    if shuf_fd_idx ==1 % making sure the first field comes first
                        this_cmp_data = [this_shuf_sample;this_small_data];
                    else
                        this_cmp_data = [this_small_data;this_shuf_sample];
                    end
                    if (IF_CROSSVAL)
                        indices = crossvalind('Kfold',this_cmp_type,num_folds);
                        this_b = nan(size(this_cmp_data,2)+1,num_folds);
                        this_pc_correct = nan(1,num_folds);
                        try
                            for n = 1:num_folds % loop through cross-validation folds
                                test = (indices == n); train = ~test;
                                [b,dev,stats] = mnrfit(this_cmp_data(train,:),this_cmp_type(train)); % Logistic regression
                                this_test_data = this_cmp_data(test,:);
                                this_test_type = this_cmp_type(test);
                                pred = nan(size(this_test_type));
                                for t = 1:size(this_test_data,1)
                                    pred(t) = this_test_data(t,:)* b(2:end)+b(1); % Use unconditional sample sizes, even for 'cond'
                                end
                                %                 [pred,dlo,dhi]= mnrval(b,this_cmp_data(test,:),stats);
                                this_pc_correct(n) = length(intersect(find(pred>0),find(this_cmp_type(test)==1)))/length(find(this_cmp_type(test)==1));
                                this_b(:,n) = b;
                            end
                            %                     this_b_idx = find(this_pc_correct == max(this_pc_correct));
                            %                     if length(this_b_idx)==1
                            %                         B = this_b(:,this_b_idx);
                            %                     else
                            %                         B = mean(this_b(:,this_b_idx),2);
                            %                     end
                            
                            B = mean(this_b,2);
                            b_shuf(fr,:,s) = B;
                            accur_shuf(fr,s)= mean(this_pc_correct);
                        catch
                            continue
                        end
                    else
                        [B,dev,stats] = mnrfit(this_cmp_data,this_cmp_type);
                        b_shuf(fr,:,s) = B;
                    end
                end
            else
                for s = 1:num_fit_shuf
                    num_folds = min(5,floor(min_num_sample/2)); % number of folds for cross-validation
                    this_shuf_idx =  randsample(1:max_num_sample,min_num_sample);
                    this_shuf_sample = this_big_data(this_shuf_idx,:);
                    if shuf_fd_idx ==1 % making sure the first field comes first
                        this_cmp_data = [this_shuf_sample;this_small_data];
                    else
                        this_cmp_data = [this_small_data;this_shuf_sample];
                    end
                    if (IF_CROSSVAL)
                        indices = crossvalind('Kfold',this_cmp_type,num_folds);
                        this_b = nan(size(this_cmp_data,2)+1,num_folds);
                        this_pc_correct = nan(1,num_folds);
                        try
                            for n = 1:num_folds % loop through cross-validation folds
                                test = (indices == n); train = ~test;
                                [b,dev,stats] = mnrfit(this_cmp_data(train,:),this_cmp_type(train)); % Logistic regression
                                this_test_data = this_cmp_data(test,:);
                                this_test_type = this_cmp_type(test);
                                pred = nan(size(this_test_type));
                                for t = 1:size(this_test_data,1)
                                    pred(t) = this_test_data(t,:)* b(2:end)+b(1); % Use unconditional sample sizes, even for 'cond'
                                end
                                %                 [pred,dlo,dhi]= mnrval(b,this_cmp_data(test,:),stats);
                                this_pc_correct(n) = length(intersect(find(pred>0),find(this_cmp_type(test)==1)))/length(find(this_cmp_type(test)==1));
                                this_b(:,n) = b;
                            end
                            %                     this_b_idx = find(this_pc_correct == max(this_pc_correct));
                            %                     if length(this_b_idx)==1
                            %                         B = this_b(:,this_b_idx);
                            %                     else
                            %                         B = mean(this_b(:,this_b_idx),2);
                            %                     end
                            
                            B = mean(this_b,2);
                            b_shuf(fr,:,s) = B;
                            accur_shuf(fr,s)= mean(this_pc_correct);
                        catch
                            continue
                        end
                    else
                        [B,dev,stats] = mnrfit(this_cmp_data,this_cmp_type);
                        b_shuf(fr,:,s) = B;
                    end
                end
            end
        end
    end
    % for the frames not used for traning threshold, use the closest
    % trained thresh
    thresh_framewise(1:frames_to_train(1)-1) = thresh_framewise(frames_to_train(1));
    thresh_framewise(frames_to_train(end)+1:end) = thresh_framewise(frames_to_train(end));
    thresh_fix = mean(thresh_framewise(frames_to_train));
else % train one classifier using frames of interests
    all_types = [ones(1,size(data.correct,1)), ones(1,size(data.incorrect,1)).*2]';
    all_data = [data.correct; data.incorrect];
    %% TRAIN CLASSIFIER
    if IF_SHUFFLE % shuffle trials for a better fit; helps to deal with unbalanced sample size
        [min_num_sample,min_fd_idx] = min([num_trial.correct num_trial.incorrect]);
%         num_folds = min(5,floor(min_num_sample/2)); % number of folds for cross-validation
        
        min_fd = fd_names{min_fd_idx};
        shuf_fd_idx = setdiff([1,2],min_fd_idx);
        shuf_fd = fd_names{shuf_fd_idx};
        max_num_sample = num_trial.(shuf_fd);
        
        this_small_data = data.(min_fd);
        this_big_data = data.(shuf_fd); %

        coef_shuf = {};
        coef0_shuf = {};
        if num_fit_shuf>10 % only use parfor when may shuffles           
            parfor s = 1:num_fit_shuf
                this_shuf_idx =  randsample(1:max_num_sample,min_num_sample);
                num_folds = min(5,floor(min_num_sample/2)); % number of folds for cross-validation
                
                % linear model
                if strcmp(MODEL_TYPE,'logistic')
                    this_shuf_sample = this_big_data(this_shuf_idx,:);
                    if ~IF_MULTI_COMPO
                        if shuf_fd_idx ==1 % making sure the correct field comes first
                            this_cmp_data = [this_shuf_sample(:);this_small_data(:)];
                            this_cmp_type = [shuf_fd_idx.*ones(size(this_shuf_sample(:)));min_fd_idx.*ones(size(this_small_data(:)))];
                        else
                            this_cmp_data = [this_small_data(:);this_shuf_sample(:)];
                            this_cmp_type = [min_fd_idx.*ones(size(this_small_data(:)));shuf_fd_idx.*ones(size(this_shuf_sample(:)))];
                            
                        end
                    else
                        if shuf_fd_idx ==1 % making sure the correct field comes first
                            this_cmp_data = [this_shuf_sample;this_small_data];
                            this_cmp_type = [shuf_fd_idx.*ones(size(this_shuf_sample,1),1);min_fd_idx.*ones(size(this_small_data,1),1)];
                            
                        else
                            this_cmp_data = [this_small_data;this_shuf_sample];
                            this_cmp_type = [min_fd_idx.*ones(size(this_small_data,1),1);shuf_fd_idx.*ones(size(this_shuf_sample,1),1)];
                        end
                    end
                    if (IF_CROSSVAL)
                        if num_folds>2
                            indices = crossvalind('Kfold',this_cmp_type,num_folds);
                            IF_LEAVEONEOUT = 0;
                        else
                            IF_LEAVEONEOUT = 1;
                            num_folds = length(this_cmp_type);
                        end
                        this_b = nan(num_comp+1,num_folds);
                        this_pc_correct = nan(1,num_folds);
                        for n = 1:num_folds % loop through cross-validation folds
                            if IF_LEAVEONEOUT
                                indices = crossvalind('LeaveMOut',this_cmp_type,1);
                            end
                            test = (indices == n); train = ~test;
                            %                         mdl = fitlm(this_cmp_data(train,:),this_cmp_type(train),'RobustOpts','logistic');
                            %                         b = mdl.Coefficients.Estimate;
                            [b] = mnrfit(this_cmp_data(train,:),this_cmp_type(train)); % Logistic regression
                            this_test_data = this_cmp_data(test,:);
                            this_test_type = this_cmp_type(test);
                            pred = nan(size(this_test_type));
                            for t = 1:size(this_test_data,1)
                                pred(t) = this_test_data(t,:)* b(2:end)+b(1); % Use unconditional sample sizes, even for 'cond'
                            end
                            %                 [pred,dlo,dhi]= mnrval(b,this_cmp_data(test,:),stats);
                            this_pc_correct(n) = length(intersect(find(pred>0),find(this_cmp_type(test)==1)))/length(find(this_cmp_type(test)==1));
                            this_b(:,n) = b;
                        end
                        %                     this_b_idx = find(this_pc_correct == max(this_pc_correct));
                        %                     if length(this_b_idx)==1
                        %                         B = this_b(:,this_b_idx);
                        %                     else
                        %                         B = mean(this_b(:,this_b_idx),2);
                        %                     end
                        
                        B = mean(this_b,2);
                        accur_shuf(s)= mean(this_pc_correct);
                        all_accur_shuf = [all_accur_shuf;this_pc_correct];
                    else
                        [B] = mnrfit(this_cmp_data,this_cmp_type);% last category will be zero
                    end
                    thresh_shuf(s) = -B(1);
                    b_shuf(1,:,s) = B;
                    
                elseif strcmp(MODEL_TYPE,'elastic')
                    this_shuf_sample = this_big_data(this_shuf_idx,:);
                    this_cmp_type = [ones(1,min_num_sample),zeros(1,min_num_sample)]';
                    
                    if shuf_fd_idx ==1 % making sure the correct field comes first
                        this_cmp_data = [this_shuf_sample;data.(min_fd)];
                    else
                        this_cmp_data = [data.(min_fd);this_shuf_sample];
                    end
                    
                    [B,stats] = lassoglm(this_cmp_data,this_cmp_type,'binomial','Alpha',0.75,'CV',3,'NumLambda',50);
                    idxLambda1 = stats.IndexMinDeviance;
                    coef = B(:,idxLambda1);
                    coef0 = stats.Intercept(idxLambda1);
                    coef_shuf{s} = coef;
                    coef0_shuf{s} = coef0;
                    
                end
                
            end
        else
            for s = 1:num_fit_shuf
                this_shuf_idx =  randsample(1:max_num_sample,min_num_sample);
                num_folds = min(5,floor(min_num_sample/2)); % number of folds for cross-validation
                
                % linear model
                if strcmp(MODEL_TYPE,'logistic')
                    this_shuf_sample = this_big_data(this_shuf_idx,:);
                    if ~IF_MULTI_COMPO
                        if shuf_fd_idx ==1 % making sure the correct field comes first
                            this_cmp_data = [this_shuf_sample(:);this_small_data(:)];
                            this_cmp_type = [shuf_fd_idx.*ones(size(this_shuf_sample(:)));min_fd_idx.*ones(size(this_small_data(:)))];
                        else
                            this_cmp_data = [this_small_data(:);this_shuf_sample(:)];
                            this_cmp_type = [min_fd_idx.*ones(size(this_small_data(:)));shuf_fd_idx.*ones(size(this_shuf_sample(:)))];
                            
                        end
                    else
                        if shuf_fd_idx ==1 % making sure the correct field comes first
                            this_cmp_data = [this_shuf_sample;this_small_data];
                            this_cmp_type = [shuf_fd_idx.*ones(size(this_shuf_sample,1),1);min_fd_idx.*ones(size(this_small_data,1),1)];
                            
                        else
                            this_cmp_data = [this_small_data;this_shuf_sample];
                            this_cmp_type = [min_fd_idx.*ones(size(this_small_data,1),1);shuf_fd_idx.*ones(size(this_shuf_sample,1),1)];
                        end
                    end
                    if (IF_CROSSVAL)
                        if num_folds>2
                            indices = crossvalind('Kfold',this_cmp_type,num_folds);
                            IF_LEAVEONEOUT = 0;
                        else
                            IF_LEAVEONEOUT = 1;
                            num_folds = length(this_cmp_type);
                        end
                        this_b = nan(num_comp+1,num_folds);
                        this_pc_correct = nan(1,num_folds);
                        for n = 1:num_folds % loop through cross-validation folds
                            if IF_LEAVEONEOUT
                                indices = crossvalind('LeaveMOut',this_cmp_type,1);
                            end
                            test = (indices == n); train = ~test;
                            %                         mdl = fitlm(this_cmp_data(train,:),this_cmp_type(train),'RobustOpts','logistic');
                            %                         b = mdl.Coefficients.Estimate;
                            [b] = mnrfit(this_cmp_data(train,:),this_cmp_type(train)); % Logistic regression
                            this_test_data = this_cmp_data(test,:);
                            this_test_type = this_cmp_type(test);
                            pred = nan(size(this_test_type));
                            for t = 1:size(this_test_data,1)
                                pred(t) = this_test_data(t,:)* b(2:end)+b(1); % Use unconditional sample sizes, even for 'cond'
                            end
                            %                 [pred,dlo,dhi]= mnrval(b,this_cmp_data(test,:),stats);
                            this_pc_correct(n) = length(intersect(find(pred>0),find(this_cmp_type(test)==1)))/length(find(this_cmp_type(test)==1));
                            this_b(:,n) = b;
                        end
                        %                     this_b_idx = find(this_pc_correct == max(this_pc_correct));
                        %                     if length(this_b_idx)==1
                        %                         B = this_b(:,this_b_idx);
                        %                     else
                        %                         B = mean(this_b(:,this_b_idx),2);
                        %                     end
                        
                        B = mean(this_b,2);
                        accur_shuf(s)= mean(this_pc_correct);
                        all_accur_shuf = [all_accur_shuf;this_pc_correct];
                    else
                        [B] = mnrfit(this_cmp_data,this_cmp_type);% last category will be zero
                    end
                    thresh_shuf(s) = -B(1);
                    b_shuf(1,:,s) = B;
                    
                elseif strcmp(MODEL_TYPE,'elastic')
                    this_shuf_sample = this_big_data(this_shuf_idx,:);
                    this_cmp_type = [ones(1,min_num_sample),zeros(1,min_num_sample)]';
                    
                    if shuf_fd_idx ==1 % making sure the correct field comes first
                        this_cmp_data = [this_shuf_sample;data.(min_fd)];
                    else
                        this_cmp_data = [data.(min_fd);this_shuf_sample];
                    end
                    
                    [B,stats] = lassoglm(this_cmp_data,this_cmp_type,'binomial','Alpha',0.75,'CV',3,'NumLambda',50);
                    idxLambda1 = stats.IndexMinDeviance;
                    coef = B(:,idxLambda1);
                    coef0 = stats.Intercept(idxLambda1);
                    coef_shuf{s} = coef;
                    coef0_shuf{s} = coef0;
                    
                end
                
            end            
        end
        thresh_fix = mean(thresh_shuf);
    else
        %         this_lda_sf = fitcdiscr(all_data,all_types,...
        %             'DiscrimType', 'linear',...
        %             'Gamma', 0, ...
        %             'FillCoeffs', 'on',...
        %             'HyperparameterOptimizationOptions',struct('UseParallel',false),...
        %             'Cost',cost); % classifier on shuffled data
        %         L = this_lda_sf.Coeffs(1,2).Linear;
        %         K = this_lda_sf.Coeffs(1,2).Const;
        %         log_thresh_shuf(s) = -K/L;
        [B] = mnrfit(all_data,all_types);
        thresh_fix = -B(1);
    end
    thresh_framewise = repmat(thresh_fix,[1,trial_length]);
end
%% save to structure
decod_struct.thresh_framewise  = thresh_framewise; % this is actually framewise cd threshold, not necassarily by logistic regression
decod_struct.thresh_fix  = thresh_fix;
% decod_struct.stats  = stats;
decod_struct.B  = squeeze(nanmean(b_shuf,3));
decod_struct.b_shuf = b_shuf;
decod_struct.accur_shuf = accur_shuf;
decod_struct.all_accur_shuf = all_accur_shuf;
if ~isempty(coef_shuf)
    decod_struct.coef_shuf = coef_shuf;
    decod_struct.coef0_shuf = coef0_shuf;
    decod_struct.coef = mean(cell2mat(coef_shuf),2);
    decod_struct.coef0 = mean(cell2mat(coef0_shuf));
end


end

