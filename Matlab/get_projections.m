function [proj_struct,add_proj_struct] = get_projections(this_struct,weights,proj_fds,varargin)
proj_struct = struct();
add_proj_struct = struct();
bias = 0;
trial_length = [];
IS_CELL_STRUCT = 0; % if this_struct is cell_struct,i.e. each row is a cell, then organise traces to [trial,time,cell]
IF_NOMINAL = 0; % if true then rescale for real probablilities
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'proj_struct')
        proj_struct = varargin{v+1};
    end
    if strcmpi(varargin{v},'bias')
        bias = varargin{v+1};
    end
    if strcmpi(varargin{v},'IS_CELL_STRUCT')
        IS_CELL_STRUCT = varargin{v+1};
    end
    if strcmpi(varargin{v},'trial_length')
        trial_length = varargin{v+1};
    end
    
    if strcmpi(varargin{v},'IF_NOMINAL')
        IF_NOMINAL = varargin{v+1};
    end
end



for fd = 1:numel(proj_fds)
    this_fd = proj_fds{fd};
    try
        if ~IS_CELL_STRUCT
            this_data = this_struct.(this_fd);
        else
            num_cells = size(this_struct,2);
            temp_data = {this_struct.(this_fd)};
            this_data = organise_data(temp_data,num_cells);
        end
        
        this_proj = nan(size(this_data,1),size(this_data,2));
        this_num_trials = size(this_data,1);this_num_bins = size(this_data,2);
        this_weights = weights;
        
        if size(weights,2)==1
            this_weights = repmat(weights,[1,this_num_bins]);
        end
        this_bias= bias;
        if length(bias)==1
            this_bias = repmat(bias,[1,this_num_bins]);
        end
        for t = 1:this_num_trials
            for f = 1:this_num_bins
                eta = squeeze(this_data(t,f,:))'* this_weights(:,f)+this_bias(f);
%                 [pred] = mnrval([bias; weights],squeeze(this_data(t,f,:))')
                if IF_NOMINAL               
                    this_proj(t,f) =  1./(1+exp(-eta))-0.5;
                else                  
                    this_proj(t,f) = eta;
                end
            end
            
        end
        if ~isempty(trial_length) % check if this_proj is [trial,frame]
            if size(this_proj,2)~= trial_length
                this_proj = this_proj';
            end
        end
        
        if this_num_bins == 1
            this_proj = this_proj';
        end
        
    catch
        this_proj =[];
    end
    
    proj_struct.(this_fd) = this_proj;
    add_proj_struct.(this_fd) = this_proj;
 
end

end

function this_data = organise_data(temp_data,num_cells)
this_num_trials = size(temp_data{1},2);
this_trial_length = size(temp_data{1},1);

this_data = nan(this_num_trials,this_trial_length,num_cells);
for c = 1:num_cells
    this_data(:,:,c) = temp_data{c}';
end

end
