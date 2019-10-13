function [proj_struct] = get_projections(this_struct,weights,proj_fds,varargin)
proj_struct = struct();
bias = 0;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'proj_struct')
        proj_struct = varargin{v+1};
    end
    if strcmpi(varargin{v},'bias')
        bias = varargin{v+1};
    end
end



for fd = 1:numel(proj_fds)
    this_fd = proj_fds{fd};
    this_data = this_struct.(this_fd);
    this_proj = nan(size(this_data,1),size(this_data,2));
    this_num_trials = size(this_data,1);this_num_bins = size(this_data,2);
    this_weights = weights;
    if size(weights,2)==1
        this_weights = repmat(weights,[1,this_num_bins]);
    end
    this_bias= bias;
    if length(bias)==1
        this_bias = repmat(bias,this_num_bins);
    end
    for t = 1:this_num_trials
        for f = 1:this_num_bins
            this_proj(t,f) = squeeze(this_data(t,f,:))'* this_weights(:,f)+this_bias(f);
        end
    end
    proj_struct.(this_fd) = this_proj;
end

end

