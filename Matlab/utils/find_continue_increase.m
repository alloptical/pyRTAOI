function [firsts,lasts] = find_continue_increase(x,min_num)
c = find(gradient(x)>0);
d = diff([0, diff(c)==1, 0]);
firsts = c(d>0);
lasts = c(d<0);
nums = lasts-firsts;
firsts = firsts(nums>min_num);
lasts = lasts(nums>min_num);

end

