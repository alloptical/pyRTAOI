% load data
s2p_data =load('C:\Users\Zihui\Dropbox\pyRTAOI-matlab\proc data\suite2p results\ex1\F_OG245_20171229_plane1_proc.mat');
s2p_data = s2p_data.dat;

%% get cell idx
iscell = arrayfun(@(x)double(cell2mat(x)),extractfield(s2p_data.stat,'iscell'));
iscell_idx = find(iscell == 1);
num_cells = numel(iscell_idx);
%% get traces
Fcell = cell2mat(s2p_data.Fcell); % cell trace
Fnp = cell2mat(s2p_data.FcellNeu); % neuropil trace
coeff = extractfield(s2p_data.stat,'neuropilCoefficient'); % neuropil coeff

Fnp = Fnp(iscell_idx,:);
Fcell = Fcell(iscell_idx,:);
Fcoeff = coeff(iscell_idx);

%% make structure
s2p_struct = struct();
for i = 1:num_cells
    % cell trace
    s2p_struct(i).Fcell = Fcell(i,:);
    s2p_struct(i).Fnp = Fnp(i,:);
    s2p_struct(i).Fsub = Fcell(i,:)-Fnp(i,:).*Fcoeff(i); % neuropil subtracted traces - could have negative values!

    % cell shape
    s2p_struct(i).xpix = s2p_data.stat(iscell_idx(i)).xpix;
    s2p_struct(i).ypix = s2p_data.stat(iscell_idx(i)).ypix;
    
    % spike trace
    spike_frames = s2p_data.stat(i).st;
    spike_amplitude = s2p_data.stat(i).c;
    s2p_struct(i).s2p_spike_trace(spike_frames) = spike_amplitude;
    s2p_struct(i).s2p_rate = sum(spike_amplitude)/length(s2p_struct(i).Fcell);
end

%% clear data - to speed up
clear s2p_data

%% show cells
roi_color = colormap(lines);
roi_color = repmat(roi_color,[5,1]);
mean_img = mean(s2p_data.mimg,3); % fov (average image)
temp_img = zeros(size(mean_img));
figure
hold on
imagesc(mean_img)
colormap(gray)
for i = 1:num_cells
    this_num_pix = numel(s2p_struct(i).xpix);
    for j = 1:this_num_pix
        temp_img(s2p_struct(i).ypix(j),s2p_struct(i).xpix(j)) = 1;
    end
    contour(temp_img,'Linecolor',roi_color(i,:))
    temp_img = zeros(size(mean_img));
end
title('Suite2p detected ROIs (after manual curation)')
axis square
%% plot traces
figure; hold on
plot_offset = mean(Fcell(:));
for i = 1:num_cells
    plot(s2p_struct(i).Fcell+i*plot_offset,'color',roi_color(i,:))
    plot(s2p_struct(i).Fsub+i*plot_offset,'color','black','linewidth',1.5)
end
xlabel('Frames')
ylabel('ROI index')