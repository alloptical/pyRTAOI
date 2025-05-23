function [b2r_colors,sort_colors,norm_factor] = get_b2r_colors(values,varargin)

norm_factor = max(abs(values));
b2r_colors = ones(length(values),3);
for i = 1:length(values)
    this_value = values(i);
    if this_value >0
        b2r_colors(i,:) = tint([1,0,0],1-this_value/norm_factor);
    elseif this_value <0
        b2r_colors(i,:) = tint([0,0,1],1+this_value/norm_factor);
        
    end
    [~,sort_idx] = sort(values);
    sort_colors = b2r_colors(sort_idx,:); % for making color bar
end

