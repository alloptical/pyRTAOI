function [p,tbl,stats] = multcmp_stats_on_struct( input_struct, varargin )
test_type = 'anova';
displayopt = 'off';
for v = 1:numel(varargin)
    if(strcmp(varargin{v},'test_type'))
        test_type = varargin{v+1};
    elseif (strcmp(varargin{v},'displayopt'))
        displayopt =  varargin{v+1};
    end
end
%
% input_struct = avg_stim;
field_names = fields(input_struct);
num_fd = numel(field_names);

tot_num_values = sum(arrayfun(@(x)numel(input_struct.(field_names{x})),1:num_fd));
group_names = cell(1,tot_num_values);
values = [];
num_proc = 1;

% prepare group names and values
for f = 1:num_fd
    this_field_name = field_names{f};
    this_num_values = numel(input_struct.(this_field_name));
    for i = num_proc:num_proc+this_num_values-1
    group_names{i} =this_field_name;
    end
    num_proc = num_proc+this_num_values;
    values = [values,input_struct.(this_field_name)(:)'];
end
if(strcmp(test_type,'anova'))
    [p,tbl,stats] = anova1(values,group_names,displayopt);
elseif (strcmp(test_type,'KW'))
    [p,tbl,stats] = kruskalwallis(values,group_names,displayopt);
elseif (strcmp(test_type,'FM'))
    [p,tbl,stats] = friedman(values,group_names,displayopt);
end

% display field names
disp(['compare fields:' ;field_names])

end

