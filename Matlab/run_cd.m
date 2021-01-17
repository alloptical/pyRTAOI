function [this_proj_struct,this_cd,this_db,this_proj_struct_framewise,this_framewise_cd] = run_cd(this_cell_struct,opt,varargin )
% see Li (Svoboda & Druckmann) motor planning Nature 2016

this_cd = [];
this_db = [];
this_framewise_cd = [];
this_proj_struct = struct();
this_proj_struct_framewise = struct();
IF_USE_DFF = false;
trace_fds = opt.fd_names; % by default, use smoothed spike trace
ref_fds = {'st_correct_stim_1_avg','st_correct_stim_2_avg'};
correct_fds = {'st_correct_stim_1','st_correct_stim_2'};
IF_NO_SCALING = false;
IF_GET_PROJ_ONLY= false;
frames_to_avg = opt.frames_to_avg;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'cd')
        this_cd = varargin{v+1};
    end
    
   if strcmpi(varargin{v},'framewise_cd')
        this_framewise_cd = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'db')
        this_db = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'IF_USE_DFF')
        IF_USE_DFF = varargin{v+1};
    end
    
   if strcmpi(varargin{v},'IF_NO_SCALING')
        IF_NO_SCALING = varargin{v+1};
   end
    
   if strcmpi(varargin{v},'ref_fds')
       ref_fds = varargin{v+1};
   end

   if strcmpi(varargin{v},'correct_fds')
       correct_fds = varargin{v+1};
   end

   if strcmpi(varargin{v},'trace_fds')
       trace_fds = varargin{v+1};
   end

   if strcmpi(varargin{v},'IF_GET_PROJ_ONLY')
       IF_GET_PROJ_ONLY = varargin{v+1};
   end

   if strcmpi(varargin{v},'this_proj_struct')
       this_proj_struct = varargin{v+1};
   end
end
if IF_USE_DFF
    trace_fds = {'st_correct_dff_1','st_incorrect_dff_1','st_correct_dff_2','st_incorrect_dff_2'};
    ref_fds = {'sta_correct_dff_1','sta_correct_dff_2'};
    correct_fds = {'st_correct_dff_1','st_correct_dff_2'};

end
% check if sta fields exist, if not get sta from st
for rf = ref_fds
    if ~any( contains(fields(this_cell_struct),rf))
        for c = 1:size(this_cell_struct,2)
         this_cell_struct(c).(rf{1}) = mean(this_cell_struct(c).(strrep(rf{1},'_avg','')),2);
        end
    end
end

% get coding direction
if ~IF_GET_PROJ_ONLY
    if IF_NO_SCALING
        % not using coding direction at all
        this_framewise_cd = ones(size(this_cell_struct,2),opt.trial_frames);
        this_cd =  ones(size(this_cell_struct,2),1);
    else
        
        if isempty(this_framewise_cd)
            traces1 = cell2mat({this_cell_struct.(ref_fds{1})})';
            traces2 = cell2mat({this_cell_struct.(ref_fds{2})})';
            
            if(IF_USE_DFF)
                traces1 = medfilt1(traces1,3,[],2);
                traces1 = movmean(traces1,[2,0],2);
                
                traces2 = medfilt1(traces2,3,[],2);
                traces2 = movmean(traces2,[2,0],2);
            end
            
            this_framewise_cd = traces1 - traces2;
        end
        
        if isempty(this_cd)
            this_cd = mean(this_framewise_cd(:,frames_to_avg),2);
        end
    end
end

% projection to coding direction
for f = 1:numel(trace_fds)
    this_fd = trace_fds{f};
    this_traces = [];
    for c = 1:size(this_cell_struct,2)
        this_trace = this_cell_struct(c).(this_fd)(:);
        if(IF_USE_DFF)
            this_trace = movmean(this_trace,[2,0]);
        end
        this_traces = [this_traces this_trace];
    end
    
    this_proj = this_traces*this_cd;
    
    this_num_trials = size(this_proj,1)/opt.trial_frames;
    this_traces_reshape = nan(this_num_trials,opt.trial_frames); %[trial,frames,factor]; in this case factor = cell
    for t = 1:this_num_trials
        this_traces_reshape(t,:) = this_proj([1:opt.trial_frames]+(t-1)*opt.trial_frames)';
    end
    
    this_proj_struct.(this_fd) = this_traces_reshape;
    
end

if IF_GET_PROJ_ONLY
    return
end
% smooth traces if using dff
proc_cell_struct = struct();
for c = 1:size(this_cell_struct,2)
    
    for f = 1:numel(trace_fds)
        this_fd = trace_fds{f};
        this_trace = this_cell_struct(c).(this_fd);
        if(IF_USE_DFF)
            this_trace = medfilt1(this_trace,3,[],2);
            this_trace = movmean(this_trace,[2,0],2);
        end      
        proc_cell_struct(c).(this_fd) = this_trace; 
    end
end
% projection to coding direction framewise
for f = 1:numel(trace_fds)
    this_fd = trace_fds{f};
    this_proj_trace = [];
    for fr = 1:size(this_framewise_cd,2)
        this_trace = cell2mat(arrayfun(@(x)x.(this_fd)(fr,:),proc_cell_struct,'UniformOutput',false)');
        this_proj = this_framewise_cd(:,fr)'*this_trace;
        this_proj_trace = [this_proj_trace;this_proj];
    end
      
    this_proj_struct_framewise.(this_fd) = this_proj_trace';
    
end

% decision boundary
if isempty(this_db)
    proj1 = this_proj_struct.(correct_fds{1});
    proj1 = nanmean(proj1(:,frames_to_avg),2);
    proj2 = this_proj_struct.(correct_fds{2});
    proj2 = nanmean(proj2(:,frames_to_avg),2);
    
    var1 = var(proj1);
    var2 = var(proj2);
    
    this_db =( mean(proj1)/var1+mean(proj2)/var2)/(1/var1+1/var2);
end
% 
% figure
% x_ticks = 1:opt.trial_frames;
% subplot(2,1,1)
% hold on
% plot_fds = {'st_correct_smooth_deconv_1','st_correct_smooth_deconv_2','st_incorrect_smooth_deconv_1','st_incorrect_smooth_deconv_2'};
% for f = 1:numel(plot_fds)
%     this_fd =plot_fds{f};
%     plot(  this_proj_struct.(this_fd)','color',color.(this_fd))
% end
% plot(xlim,[this_db this_db],'color','black')
% xlim([50,120])
% xlabel('Frames')
% ylabel('Projection to CD')
% 
% 
% 
% subplot(2,1,2)
% hold on
% plot_fds = fields(this_proj_struct);
% for f = 1:numel(plot_fds)
%     this_fd =plot_fds{f};
%     this_traces =  this_proj_struct.(this_fd)-this_db;
% %     shadedErrorBar(x_ticks,mean(this_traces,1),...
% %         std(this_traces,[],1),{'color',color.(this_fd),'linewidth',2},0.1)
%     shadedErrorBar(x_ticks,median(this_traces,1),[quantile(this_traces,0.75)-median(this_traces,1);...
%         median(this_traces,1)-quantile(this_traces,0.25)],{'color',color.(this_fd),'linewidth',2},0.5)
%     
% end
% 
% 
% plot(xlim,[0 0],'color','black')
% xlim([50,120])
% xlabel('Frames')
% ylabel('Distance to DB')





end

