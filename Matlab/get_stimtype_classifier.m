function [classifier] = get_stimtype_classifier(this_struct,windows,this_opt)
% train classifier using frames where the correct trial has entered steady state
% steady states are defined by 'windows'
% random shuffle trials to balance number of trials
% k-fold cross validation for each shuffle
% logistic regression

fds = this_opt.fd_names;
num_comp = this_opt.m;
IF_SHUFFLE = true;
num_fit_shuf = 50;

% initiation
classifier = struct();
data = struct();
num_trial = struct();
correct_fd ={fds{contains(fds, '_correct_')}};
for fd_idx = 1:numel(correct_fd) 
    this_correct_fd = correct_fd{fd_idx};
    data.(this_correct_fd) = this_struct.(this_correct_fd);  % [trial, frame, components]
    [tr_idx.(this_correct_fd)] = find(windows.(this_correct_fd)>0);
    data_fds = fields(data);
    num_trial.(this_correct_fd) = size(this_struct.(this_correct_fd),1);    
end


if IF_SHUFFLE % shuffle trials for a better fit; helps to deal with unbalanced sample size
    [min_num_sample,min_fd_idx] = min( cellfun(@(x)getfield(num_trial,x),correct_fd));
    num_folds = min(5,floor(min_num_sample/2)); % number of folds for cross-validation
    disp(['training classifier, num folds ' num2str(num_folds)])
    tic
    min_fd = data_fds{min_fd_idx};
    shuf_fd_idx = setdiff([1,2],min_fd_idx);
    shuf_fd = data_fds{shuf_fd_idx};
    max_num_sample = num_trial.(shuf_fd);
    shuf_tr_idx = tr_idx.(shuf_fd);
    small_tr_idx = tr_idx.(min_fd);
    
    this_small_data = data.(min_fd);
    this_small_sample= [];
    for i = 1:num_comp
        temp = this_small_data(:,:,i);
        this_small_sample = [this_small_sample,temp(small_tr_idx)];
    end
    this_big_data = data.(shuf_fd); %
    num_comp = size(this_big_data,3);
    for i = 1:num_comp
        temp =  this_big_data(:,:,i);
        temp(shuf_tr_idx) = nan;
        this_big_data(:,:,i) = temp;
    end
    
    weights_shuf = nan(num_comp+1,num_fit_shuf);
    accur_shuf = nan(1,num_fit_shuf);

    
    for s = 1:num_fit_shuf % shuffle to downsample big data
        this_shuf_idx =  randsample(1:max_num_sample,min_num_sample);
        this_shuf_sample = this_big_data(this_shuf_idx,:,:);
        this_shuf_sample = this_shuf_sample(~isnan(this_shuf_sample));
        this_shuf_sample = reshape(this_shuf_sample,[length(this_shuf_sample)/num_comp,num_comp]);
        if shuf_fd_idx ==1 % making sure the correct field comes first
            this_cmp_data = [this_shuf_sample;this_small_sample];
            this_cmp_type = [ones(length(this_shuf_sample),1);ones(length(this_small_sample),1).*2];
            
        else
            this_cmp_data = [this_small_sample;this_shuf_sample];
            this_cmp_type = [ones(length(this_small_sample),1);ones(length(this_shuf_sample),1).*2];
            
        end
        
        % train classifier with cross validation
        indices = crossvalind('Kfold',this_cmp_type,num_folds);
        this_b = nan(num_comp+1,num_folds);
        this_pc_correct = nan(1,num_folds);
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
        B = mean(b,2);
        accur_shuf(s)= mean(this_pc_correct);
        weights_shuf(:,s) = B;
    end
    toc
    disp('... Done')
end
classifier.B = mean(weights_shuf,2);
classifier.accur= mean(accur_shuf);
classifier.B_raw  = weights_shuf;
classifier.accur_raw = accur_shuf;

end

