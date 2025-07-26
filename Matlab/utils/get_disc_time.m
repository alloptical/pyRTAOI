function [ decod_struct ] =  get_disc_time( proj_struct, decod_struct,opt,varargin )
% adapted from get_classfier_dc_time in dropbox 
Nstd = getOr(opt,'Nstd',1.5);
min_frames = getOr(opt,'min_frames',10);
try
    num_frames = length(opt.frame_range);
catch
    num_frames = opt.trial_length;

end
num_cv_shuf = 100;
IF_USE_DB = false; % use decision boundary as thresh
IF_FRAMEWISE = getOr(opt,'IF_FRAMEWISE',false);
IF_REVERSE = false;
IF_USE_CROSSVAL_RESULT = false;
threshold = [];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'db')
        this_thresh = varargin{v+1};
        IF_USE_DB = true;
    end
    
    if strcmpi(varargin{v},'IF_FRAMEWISE')
        IF_FRAMEWISE = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'threshold')
        threshold = varargin{v+1};
    end

   if strcmpi(varargin{v},'IF_REVERSE')
        IF_REVERSE = varargin{v+1};
   end
    
   if strcmpi(varargin{v},'IF_USE_CROSSVAL_RESULT')
       IF_USE_CROSSVAL_RESULT = varargin{v+1};
   end


end


this_struct = proj_struct;
fds = fields(this_struct);

% initiation
correct_fd ={fds{contains(fds, '_correct_')}};
incorrect_fd = {fds{contains(fds, '_incorrect_')}};

for stim_type = 1:2
    this_fd =  ['stim' num2str(stim_type)];
    if isempty(threshold)
        if ~IF_USE_DB
            if IF_FRAMEWISE
                this_thresh = decod_struct.thresh_framewise.(this_fd);
            else
                this_thresh = decod_struct.thresh_fix.(this_fd);
            end
        end
    else
        this_thresh= threshold.(this_fd);
        decod_struct.thresh_framewise.(this_fd) = this_thresh;
    end
    if length(this_thresh)==1
        this_thresh = repmat(this_thresh,[1,num_frames]);
    end
    
    this_correct_fd = correct_fd{contains(correct_fd, ['_' num2str(stim_type)])};
    this_incorrect_fd = incorrect_fd{contains(incorrect_fd, ['_' num2str(stim_type)])};
    data.correct = this_struct.(this_correct_fd);
    data.incorrect =  this_struct.(this_incorrect_fd);
    
    num_sample.correct = size(data.correct,1);
    num_sample.incorrect = size(data.incorrect,1);
    tot_num_samples = num_sample.correct+num_sample.incorrect;
    hit_rate = nan(1,num_frames);
    fa_rate = nan(1,num_frames);
    classif_accuracy = nan(1,num_frames);
    shuf_classif_accuracy = {};
    shuf_classif_mean = nan(1,num_frames);
    shuf_classif_sd = nan(1,num_frames);

    for f = 1:num_frames
        this_proj = [data.correct(:,f); data.incorrect(:,f)];
        this_proj = this_proj-this_thresh(f);
        if IF_REVERSE
            this_proj = (-1)^(stim_type-1).*this_proj;
        end
        
        this_yes =  find(this_proj <0);
        this_no = find(this_proj>0);
        
        this_positive =  num_sample.correct+1:num_sample.correct+num_sample.incorrect; % incorrect trial indices
        this_negative = 1:num_sample.correct; % correct trial indices
        hit_rate(f) = numel(intersect(this_positive,this_yes))/numel(this_positive);
        fa_rate(f)= numel(intersect(this_negative,this_yes))/numel(this_negative);
        
        
        %% SHUFFLE DATA TO GET CLASSIFICATION ACCURACY EXPECTED BY CHANCE
        if IF_USE_CROSSVAL_RESULT
            this_shuf_accuracy = decod_struct.all_accur_shuf.(this_fd)(:);
        else
            this_shuf_accuracy = nan(1,num_cv_shuf);
            for s = 1:num_cv_shuf
                shuf_positive = randsample(1:tot_num_samples,num_sample.incorrect); %shuf incorrect trial indices
                shuf_negtive = setdiff(1:tot_num_samples,shuf_positive);
                this_shuf_accuracy(s) =  (numel(intersect(shuf_positive,this_yes))+ numel(intersect(shuf_negtive,this_no)))/tot_num_samples;
            end
        end
        shuf_classif_accuracy{f} = this_shuf_accuracy;
        shuf_classif_mean(f) = mean(this_shuf_accuracy);
        shuf_classif_sd(f) = std(this_shuf_accuracy);
        num_correct =  numel(intersect(this_positive,this_yes))+ numel(intersect(this_negative,this_no));
        classif_accuracy(f) = num_correct/tot_num_samples; % number correctly classified/number samples

    end
    
    high_hit_frames = find(hit_rate>0.7);
    low_fa_frames = find(fa_rate<0.3);
    disc_frames = intersect(high_hit_frames,low_fa_frames);
    disc_frame_idx = pt_continuousabove(2-diff(disc_frames),0,0,min_frames,1000,0);
    if ~isempty(disc_frame_idx)
        disc_frame = disc_frames(disc_frame_idx(1));
    else
        disc_frame = nan;
    end

    
    shuf_disc_frames = find(classif_accuracy>shuf_classif_mean+Nstd.*shuf_classif_sd);
    shuf_disc_frame = pt_continuousabove(classif_accuracy-shuf_classif_mean-Nstd.*shuf_classif_sd,0,0,min_frames,1000,0);
    if ~isempty(shuf_disc_frame)
        shuf_disc_frame = shuf_disc_frame(1);
    else
        shuf_disc_frame = nan;
    end
%     shuf_disc_frames = shuf_disc_frames(shuf_disc_frames>shuf_disc_frame);
    decod_struct.(this_fd).framewise_hr = hit_rate;
    decod_struct.(this_fd).framewise_fa = fa_rate;
    decod_struct.(this_fd).disc_frames = disc_frames; % discrimination frames by hit rate and false alarm rate
    decod_struct.(this_fd).disc_frame = disc_frame;
    decod_struct.(this_fd).classif_accuracy = classif_accuracy;
    decod_struct.(this_fd).shuf_classif_accuracy = shuf_classif_accuracy;
    decod_struct.(this_fd).shuf_classif_sd = shuf_classif_sd;
    decod_struct.(this_fd).shuf_classif_mean = shuf_classif_mean;
    decod_struct.(this_fd).shuf_disc_frames = shuf_disc_frames; % discrimination frames by shuffling
    decod_struct.(this_fd).shuf_disc_frame = shuf_disc_frame;
end

end

%% temp tests
%{
incorrect_after_dc = data.incorrect(:,dc_frame:end)-sel_thresh_fix;
tresh_cr_frame = zeros(1,size(incorrect_after_dc,1));
for ff = 1:size(incorrect_after_dc,1)
    [a,b] = find(incorrect_after_dc(ff,:)<0,1);
    if ~isempty(a)
    tresh_cr_frame(ff) = b + dc_frame;
    end
end

figure
hold on
plot(data.correct','color',[0 0 0]);
plot(data.incorrect','color',[.7 .7 .7]);
plot(dc_frame:opt.trial_frames,cv_thresh(dc_frame:end),'color',[1 0 0],'linewidth',2)

figure
hold on
plot(data.incorrect','color',[.7 .7 .7]);
plot(dc_frame:opt.trial_frames,cv_thresh(dc_frame:end),'color',[0 0 0],'linewidth',2)
plot(xlim,[sel_thresh_fix sel_thresh_fix],'color',[.5 .5 .5],'linewidth',2)
scatter(tresh_cr_frame,data.incorrect(:,tresh_cr_frame),'MarkerEdgeColor',[1 0 0])

figure

this_data = data.correct;
this_diff = diff(this_data,[],2);
this_switch = this_data.*([this_data(:,2:end),this_data(:,1)])<0;
subplot(2,2,1)
hold on
imagesc(this_data);
colormap(b2r(-1, 1))
plot([dc_frame dc_frame],ylim,':','color',[0 0 0],'linewidth',2)
ylabel('Correct trials')

subplot(2,2,2)
hold on
imagesc(this_diff);
% colormap(b2r(-1, 1))

subplot(2,2,3)
this_data = data.incorrect;
this_diff = diff(this_data,[],2);
this_switch = this_data.*([this_data(:,2:end),this_data(:,1)])<0;

hold on
imagesc(data.incorrect);
colormap(b2r(-1, 1))
ylabel('Incorrect trials')
plot([dc_frame dc_frame],ylim,':','color',[0 0 0],'linewidth',2)
xlabel('Frames')

subplot(2,2,4)
hold on
imagesc(this_diff);
% colormap(b2r(-1, 1))

%}

