function [firsts,lasts] = find_continue(x,min_num,type)

if strcmp(type,'increase')
    c = find(gradient(x)>0);
    d = diff([0, diff(c)==1, 0]);
elseif strcmp(type,'decrease')
    c = find(gradient(x)<0);
    d = diff([0, diff(c)==1, 0]);
elseif strcmp(type,'const')
    c = 1:length(x);
    d = diff([0, diff(x)==0, 0]);
end
firsts = c(d>0);
lasts = c(d<0);
nums = lasts-firsts;
firsts = firsts(nums>min_num);
lasts = lasts(nums>min_num);

end

