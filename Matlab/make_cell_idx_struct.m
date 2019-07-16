function [cell_idx_struct,outcomeAUC_zscore,correct_stimulusAUC,correct_texAUC,pcorrect,pincorrect] = make_cell_idx_struct(cell_struct,opt)
%% this is a messy function testing different ways to define stimulus and port cells 
% add fields to cell_idx_struct
num_cells = size(cell_struct,2);
avg_frame_range = opt.withold_frames_adj;
bs_frame_range  = getOr(opt,'bs_frame_range', 20:40); % 
pcorrect = nan(1,num_cells);
pincorrect = nan(1,num_cells);
pbaseline = nan(1,num_cells);
IfFlip =  nan(1,num_cells);
this_bs = {};

% AUC initiation
num_shuf = 100;
N = getOr(opt,'N',1.5);

outcomeAUC{1} = nan(1,num_cells);
outcomeAUC{2} = nan(1,num_cells);
outcomeAUC_zscore = outcomeAUC; % if cell preferrably respond in correct or incorrect trials for each texture

correct_texAUC = outcomeAUC; % if texture evoked activity is lager than baseline in correct trials
incorrect_texAUC = outcomeAUC; % if texture evoked activity is lager than baseline in incorrect trials

correct_stimulusAUC = nan(1,num_cells);
correct_stimulusAUC_zscore = nan(1,num_cells); % if cell preferrably respond to one texture in correct trials
incorrect_stimulusAUC_zscore = nan(1,num_cells); % if cell preferrably respond to one texture in incorrect trials

%% using inferred spike rate
for i = 1:num_cells

    % correct stimulus AUC
    all_stim1 = [mean(cell_struct(i).('st_correct_stim_1')(avg_frame_range,:),1)];
    all_stim2 = [mean(cell_struct(i).('st_correct_stim_2')(avg_frame_range,:),1)];
    labels = [ones(1,length(all_stim1)),2.*ones(1,length(all_stim2))]';
    scores = [ all_stim1  all_stim2]';
    [~,~,~, correct_stimulusAUC(i)] = perfcurve(labels,scores,1);
    
    %% shuffle to get zscore stim auc
    shuf_stim_auc = nan(1,num_shuf);
    for s = 1:num_shuf
        shuf_labels = labels(randperm(length(labels))');
        [~,~,~, shuf_stim_auc(s)] = perfcurve(shuf_labels,scores,1);
    end
    correct_stimulusAUC_zscore(i) = (correct_stimulusAUC(i)-mean(shuf_stim_auc))/std(shuf_stim_auc);

    % incorrect stimulus AUC
    all_stim1 = [mean(cell_struct(i).('st_incorrect_stim_1')(avg_frame_range,:),1)];
    all_stim2 = [mean(cell_struct(i).('st_incorrect_stim_2')(avg_frame_range,:),1)];
    labels = [ones(1,length(all_stim1)),2.*ones(1,length(all_stim2))]';
    scores = [ all_stim1  all_stim2]';
    [~,~,~, incorrect_stimulusAUC(i)] = perfcurve(labels,scores,1);
    
    %% shuffle to get zscore stim auc
    shuf_stim_auc = nan(1,num_shuf);
    for s = 1:num_shuf
        shuf_labels = labels(randperm(length(labels))');
        [~,~,~, shuf_stim_auc(s)] = perfcurve(shuf_labels,scores,1);
    end
    incorrect_stimulusAUC_zscore(i) = (incorrect_stimulusAUC(i)-mean(shuf_stim_auc))/std(shuf_stim_auc);
    
    % separating stimulus types
    for stim_type = 1:2 %
        
        % using trial mean
        this_correct{stim_type} =  mean(cell_struct(i).(['st_correct_stim_'  num2str(stim_type)])(avg_frame_range,:),1);
        this_incorrect{stim_type} =  mean(cell_struct(i).(['st_incorrect_stim_'  num2str(stim_type)])(avg_frame_range,:),1);
        
        % using trial max
        this_correct{stim_type} =  max(cell_struct(i).(['st_correct_stim_'  num2str(stim_type)])(avg_frame_range,:),[],1);
        this_incorrect{stim_type} =  max(cell_struct(i).(['st_incorrect_stim_'  num2str(stim_type)])(avg_frame_range,:),[],1);
        
        this_bs{stim_type} = mean(cell_struct(i).(['st_correct_stim_'  num2str(stim_type)])(bs_frame_range,:),1);
        this_active{stim_type} = mean(cell_struct(i).(['st_correct_stim_'  num2str(stim_type)])(avg_frame_range,:),1);
        labels = [ones(1,length(this_bs{stim_type})),2.*ones(1,length(this_active{stim_type}))]';
        scores = [ this_bs{stim_type}  this_active{stim_type}]';
        [~,~,~,correct_texAUC{stim_type}(i)] = perfcurve(labels,scores,2);
        
        this_bs{stim_type} = mean(cell_struct(i).(['st_incorrect_stim_'  num2str(stim_type)])(bs_frame_range,:),1);
        this_active{stim_type} = mean(cell_struct(i).(['st_incorrect_stim_'  num2str(stim_type)])(avg_frame_range,:),1);
        labels = [ones(1,length(this_bs{stim_type})),2.*ones(1,length(this_active{stim_type}))]';
        scores = [ this_bs{stim_type}  this_active{stim_type}]';
        [~,~,~,incorrect_texAUC{stim_type}(i)] = perfcurve(labels,scores,2);
        
        % AUC (correct vs incorrect)
        labels = [ones(1,length(this_correct{stim_type})),2.*ones(1,length(this_incorrect{stim_type}))]';
        scores = [ this_correct{stim_type}  this_incorrect{stim_type}]';
        [~,~,~,outcomeAUC{stim_type}(i)] = perfcurve(labels,scores,1); 
        
        % shuffle to get zscore stim auc
        shuf_outcome_auc = nan(1,num_shuf);
        for s = 1:num_shuf
            shuf_labels = labels(randperm(length(labels))');
            [~,~,~, shuf_outcome_auc(s)] = perfcurve(shuf_labels,scores,1);
        end
        outcomeAUC_zscore{stim_type}(i) = (outcomeAUC{stim_type}(i)-mean(shuf_outcome_auc))/std(shuf_outcome_auc);
        
    end
    

    pcorrect(i) = ranksum(this_correct{1},this_correct{2});
    pincorrect(i) = ranksum(this_incorrect{1},this_incorrect{2});
    pbaseline(i) = ranksum( cell2mat(this_bs),cell2mat(this_active));
    IfFlip(i) = ((mean(this_correct{1})- mean(this_correct{2}))*( mean(this_incorrect{1})- mean(this_incorrect{2}))<0);
end


% using p value
port_idx = unique(find(pcorrect<0.05  & pincorrect<0.2 & IfFlip==1))';
relavent_idx =  unique(find(pbaseline <0.01))';
stim_idx = unique(find(pcorrect<0.05  & pincorrect<0.2 & IfFlip==0))';
trialcode_idx = find(pcorrect<0.05);

% using shuf auc 20190516
trialcode_idx = find(abs(correct_stimulusAUC_zscore)>N);
stim_idx =  find(abs(correct_stimulusAUC_zscore)>N&abs(incorrect_stimulusAUC_zscore)>N...
    &correct_stimulusAUC_zscore.*incorrect_stimulusAUC_zscore>0);
port_idx =  find(abs(correct_stimulusAUC_zscore)>N&abs(incorrect_stimulusAUC_zscore)>N...
    &correct_stimulusAUC_zscore.*incorrect_stimulusAUC_zscore<0);

% cells prefering one texture in correct trials (potential photostim targets)
tex1 = find(correct_stimulusAUC_zscore>N& correct_texAUC{1}>0.6); % cells prefering texture1 in correct trials
tex2 = find(correct_stimulusAUC_zscore<-N& correct_texAUC{2}>0.6); % cells prefering texture2 in correct trials
tex = unique([tex1,tex2]);

% cells showed different response in correct and incorrect trials (for identify incorrect trials)
outcome1 = find(abs(outcomeAUC_zscore{1})>N & correct_texAUC{1}>0.6);
outcome2 = find(abs(outcomeAUC_zscore{2})>N & correct_texAUC{2}>0.6);

cell_idx_struct.relavent = relavent_idx;
cell_idx_struct.port = port_idx; % this is actually port cells
cell_idx_struct.stim = stim_idx;
% cell_idx_struct.outcome1 = intersect(find(abs(outcomeAUC{1}-0.5)>0.2),trialcode_idx);
% cell_idx_struct.outcome2 = intersect(find(abs(outcomeAUC{2}-0.5)>0.2),trialcode_idx);
cell_idx_struct.tex1 = tex1;
cell_idx_struct.tex2 = tex2;
cell_idx_struct.tex = tex;


cell_idx_struct.outcome1 = outcome1;
cell_idx_struct.outcome2 = outcome2;
cell_idx_struct.outcome = unique([cell_idx_struct.outcome1 cell_idx_struct.outcome2]);
cell_idx_struct.trialcode = trialcode_idx;


end

