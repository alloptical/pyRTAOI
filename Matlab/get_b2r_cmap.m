function [cmap] = get_b2r_cmap(cell_avg)
fdd = fields(cell_avg);
fd1 = fdd(contains(fdd,'_1'));
all_val = [];
for f = 1:numel(fd1)
    all_val = [all_val;cell_avg.(fd1{f})];
end
% cmap1 = b2r(min(all_val),max(all_val));
% cmap1 = imresize(cmap1,[length(all_val),3]);
% [~,idx] = sort(all_val);
cmap1 =1- ones(length(all_val),3).*all_val./max(all_val);
all_val = abs(all_val);


fd1 = fdd(contains(fdd,'_2'));
all_val = [];
for f = 1:numel(fd1)
    all_val = [all_val;cell_avg.(fd1{f})];
end
all_val = abs(all_val);
cmap2 =1- ones(length(all_val),3).*all_val./max(all_val);

cmap = [cmap1;cmap2];

end

