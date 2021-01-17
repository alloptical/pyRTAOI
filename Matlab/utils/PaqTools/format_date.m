function [varargout] = format_date()
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

x = clock;

if nargout < 2

    varargout{1} = [num2str(x(1)) sprintf('%2.2d',x(2)) sprintf('%2.2d',x(3))];
    
else

    varargout{1} = [num2str(x(1)) sprintf('%2.2d',x(2)) sprintf('%2.2d',x(3))];
    
    varargout{2} = [sprintf('%2.2d',x(4)) '-' sprintf('%2.2d',x(5))];
    
end

end
