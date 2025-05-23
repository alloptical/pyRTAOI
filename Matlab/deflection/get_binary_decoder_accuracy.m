function [ decod_struct ] =  get_binary_decoder_accuracy( proj_struct, decod_struct,fds,frames_to_avg,varargin )
% adapted from get_classfier_dc_time in dropbox
% apply decoder threshold to frame averages
num_cv_shuf = 100;
IF_REVERSE = false;
threshold = [];
for v = 1:numel(varargin)

    if strcmpi(varargin{v},'threshold')
        threshold = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'IF_REVERSE')
        IF_REVERSE = varargin{v+1};
    end

end


this_struct = proj_struct;

% initiation
this_correct_fd =fds{1};
this_incorrect_fd = fds{2};

% get threshold
if isempty(threshold)
    this_thresh = decod_struct.thresh_fix;
else
    this_thresh= threshold;
    decod_struct.thresh_framewise = this_thresh;
end


data.correct = this_struct.(this_correct_fd);
data.incorrect =  this_struct.(this_incorrect_fd);

num_sample.correct = size(data.correct,1);
num_sample.incorrect = size(data.incorrect,1);
tot_num_samples = num_sample.correct+num_sample.incorrect;

this_proj = [mean(data.correct(:,frames_to_avg),1); mean(data.incorrect(:,frames_to_avg),1)];
this_proj = this_proj-this_thresh;
if IF_REVERSE
    this_proj = (-1).*this_proj;
end

this_yes =  find(this_proj <0);
this_no = find(this_proj>0);

this_positive =  num_sample.correct+1:num_sample.correct+num_sample.incorrect; % incorrect trial indices
this_negative = 1:num_sample.correct; % correct trial indices
hit_rate = numel(intersect(this_positive,this_yes))/numel(this_positive);
fa_rate= numel(intersect(this_negative,this_yes))/numel(this_negative);


%% SHUFFLE DATA TO GET CLASSIFICATION ACCURACY EXPECTED BY CHANCE
this_shuf_accuracy = nan(1,num_cv_shuf);
for s = 1:num_cv_shuf
    shuf_positive = randsample(1:tot_num_samples,num_sample.incorrect); %shuf incorrect trial indices
    shuf_negtive = setdiff(1:tot_num_samples,shuf_positive);
    this_shuf_accuracy(s) =  (numel(intersect(shuf_positive,this_yes))+ numel(intersect(shuf_negtive,this_no)))/tot_num_samples;
end
shuf_classif_accuracy = this_shuf_accuracy;
shuf_classif_mean = mean(this_shuf_accuracy);
shuf_classif_sd = std(this_shuf_accuracy);
num_correct =  numel(intersect(this_positive,this_yes))+ numel(intersect(this_negative,this_no));
classif_accuracy = num_correct/tot_num_samples; % number correctly classified/number samples



% shuf_disc_frames = shuf_disc_frames(shuf_disc_frames>shuf_disc_frame);
decod_struct.hit_rate = hit_rate;
decod_struct.fa_rate = fa_rate;
decod_struct.classif_accuracy = classif_accuracy;
decod_struct.shuf_classif_accuracy = shuf_classif_accuracy;
decod_struct.shuf_classif_sd = shuf_classif_sd;
decod_struct.shuf_classif_mean = shuf_classif_mean;



end

