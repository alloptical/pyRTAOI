function [dp,hr,fa] = get_dprime(num_hit,num_miss,num_fa,num_cr)
% d'=Z(Hit/(Hit+Miss))-Z(FA/(FA+CR))   Chen 2013
hr = num_hit/(num_hit+num_miss);% hit rate
fa = num_fa/(num_fa+num_cr); % false alarm rate
dp = norminv(hr)-norminv(fa);
end

