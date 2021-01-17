function [output_struct] = parse_rtaoi_result_name(file_name)
% a very unflexible function..
temp_names = strsplit(file_name,'_');
result_type = temp_names{1};
tseries = [];
date = temp_names{2};
animal = temp_names{3};
proc_time = strrep(temp_names{end},'.mat','');

for i = 1:numel(temp_names) % add other stuff to parse here.
    if strcmpi(temp_names{i},'t')
        tseries = temp_names{i+1};
    end
end

output_struct = struct();
output_struct.animal = animal;
output_struct.date = date;
output_struct.tseries = tseries;
output_struct.proc_time = proc_time;
output_struct.result_type = result_type;
end

