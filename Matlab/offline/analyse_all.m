%% close and clear
close all
clear all

%% choose example
example = 'ex1';
save_figs = true;
lenient_onacid = 0;
use_cnmf = 0;

%% Suite 2p results: load and process 
%% load data
folder = fullfile('C:\Users\intrinsic\Desktop\pyRTAOI-matlab\proc data\suite2p results', example);
cd(folder)
proc_file = dir('*proc.mat');
filename = proc_file.name;

s2p_data = load(filename);
s2p_data = s2p_data.dat;

cd('C:\Users\intrinsic\Desktop\pyRTAOI-matlab')
addpath(genpath(cd))

%% get cell idx
iscell_ = arrayfun(@(x)double(cell2mat(x)),extractfield(s2p_data.stat,'iscell'));
iscell_idx = find(iscell_ == 1);
num_cells_s2p = numel(iscell_idx);
%% get traces
Fcell = cell2mat(s2p_data.Fcell); % cell trace
Fnp = cell2mat(s2p_data.FcellNeu); % neuropil trace
coeff = extractfield(s2p_data.stat,'neuropilCoefficient'); % neuropil coeff

Fnp = Fnp(iscell_idx,:);
Fcell = Fcell(iscell_idx,:);
Fcoeff = coeff(iscell_idx);

%% make structure
s2p_struct = struct();
for i = 1:num_cells_s2p
    % cell trace
    s2p_struct(i).Fcell = Fcell(i,:);
    s2p_struct(i).Fnp = Fnp(i,:);
    s2p_struct(i).Fsub = Fcell(i,:)-Fnp(i,:).*Fcoeff(i); % neuropil subtracted traces - could have negative values!

    % cell shape
    s2p_struct(i).xpix = s2p_data.stat(iscell_idx(i)).xpix;
    s2p_struct(i).ypix = s2p_data.stat(iscell_idx(i)).ypix;    
end

%% extract cell contours
mean_img = mean(s2p_data.mimg,3); % fov (average image)
temp_img = zeros(size(mean_img));

% % CC_suite2p = cell(num_cells,1);
CC_suite2p = cell(num_cells_s2p,1);

for i = 1:num_cells_s2p
    this_num_pix = numel(s2p_struct(i).xpix);
    for j = 1:this_num_pix
        temp_img(s2p_struct(i).ypix(j),s2p_struct(i).xpix(j)) = 1;
    end
    
    CC_suite2p{i} = contour(temp_img,'Linecolor', 'r');
    
% %     [x,y] = find(temp_img == 1);
% %     CC_suite2p{i}(1,:) = x';
% %     CC_suite2p{i}(2,:) = y';
    temp_img = zeros(size(mean_img));
end
close gcf

%% Offset suite2p contours for original pixel locations - make sure first is x and second y!
x_pix = s2p_data.ops.xrange;
y_pix = s2p_data.ops.yrange;

% offsets
x_lower = x_pix(1);
x_upper = x_pix(end);

y_lower = y_pix(1);
y_upper = y_pix(end);

for i=1:num_cells_s2p
    CC_suite2p{i}(1,:) = CC_suite2p{i}(1,:)+x_lower;
    CC_suite2p{i}(2,:) = CC_suite2p{i}(2,:)+y_lower;
end

%% create binary mask for suite2p
separate_cells = 1;
mean_img = mean(s2p_data.mimg,3); % fov (average image)
if separate_cells
    mask_suite2p = zeros([num_cells_s2p,512,512]); %size(mean_img)]);
else
    mask_suite2p = zeros([512,512]);%(size(mean_img));
end

for i = 1:num_cells_s2p
    this_num_pix = numel(s2p_struct(i).xpix);
    for j = 1:this_num_pix
        if separate_cells
            mask_suite2p(i,s2p_struct(i).ypix(j)+y_lower,s2p_struct(i).xpix(j)+x_lower) = 1;
        else
            mask_suite2p(s2p_struct(i).ypix(j)+y_lower,s2p_struct(i).xpix(j)+x_lower) = 1;
        end
    end
end

%% Cnmf results: load and process
%% load data
if use_cnmf
    folder = fullfile('C:\Users\intrinsic\Desktop\pyRTAOI-matlab\proc data\cnmf results', example);
    cd(folder)

    proc_file = dir('*.mat');
    filename = proc_file(1).name;

    cnmf_data = load(filename);

    cd('C:\Users\intrinsic\Desktop\pyRTAOI-matlab')

    %% make data structure
    cnmf_struct = struct();
    cnmf_dims = cnmf_data.cnm_dims;
    cnmf_image = reshape(cnmf_data.cnm_b,cnmf_dims);
    cnmf_A = full(cnmf_data.cnm_A);

    if iscell(cnmf_A)
        cnmf_A = cell2mat(cnmf_A);
    end

    try
        num_frames_cnmf = cnmf_data.t_cnm;
        full_cnmf = false;
    catch
        f_size = size(cnmf_data.cnm_f);
        num_frames_cnmf = f_size(2);
        full_cnmf = true;
    end

    num_comp_cnmf = size(cnmf_A,2);
    comp_shape_cnmf = [cnmf_dims(1),cnmf_dims(2),num_comp_cnmf];

    try
        cm_cnmf = cnmf_data.coms;
    catch
        cm_cnmf = com(sparse(double(cnmf_A)),cnmf_dims(1),cnmf_dims(2)); % center of mass
    end


    for i = 1:num_comp_cnmf
        cnmf_struct(i).shape = reshape(cnmf_A(:,i),cnmf_dims);
        cnmf_struct(i).noisyC = cnmf_data.noisyC(i,1:num_frames_cnmf);
        cnmf_struct(i).deconvC = cnmf_data.cnm_C(i,1:num_frames_cnmf);  % deconvC
        cnmf_struct(i).centroid = cm_cnmf(i,:);
    end

    %% extract contours
    cnm_plot_options = CNMFSetParms;
    cnm_plot_options.roi_color = colormap(lines(num_comp_cnmf));

    CC_cnmf = plot_contours(sparse(double(cnmf_A)),cnmf_image,cnm_plot_options,1,[],[],[1 1 1]);
    close gcf

    %% create (binary) cell mask
    thresh = 2e-4;
    binary = false;

    if separate_cells
        mask_cnmf = zeros([num_comp_cnmf,512,512]);
        for i=1:num_comp_cnmf
            mask = reshape(cnmf_A(:,i), cnmf_dims);
            if binary
                mask(mask>thresh) = 1;
            end
            mask_cnmf(i,:,:) = mask;
        end
    else
        mask_cnmf = reshape(mean(cnmf_A,2), cnmf_dims);
        if binary
            mask_cnmf(mask_cnmf>thresh) = 1;
        end
    end
    
    thresh = 2e-4;
    binary_mask_cnmf = mask_cnmf;
    binary_mask_cnmf(binary_mask_cnmf>thresh) = 1;
end

%% Onacid results: load and process
%% load data
folder = fullfile('C:\Users\intrinsic\Desktop\pyRTAOI-matlab\proc data\onacid results', example);
cd(folder)

if lenient_onacid
    proc_file = dir('lenient*.mat');
else
    proc_file = dir('ex*onacid*.mat');
    if length(proc_file)==0
        proc_file = dir('onacid*.mat');
    end
end

if strcmp(example, 'ex3')
    filename = proc_file(2).name;
else
    filename = proc_file(1).name;
end

% proc_file = dir('*.mat');

% if strcmp(example, 'ex3')
%     filtered = true;
% else
%     filtered = false;  % use raw for most; use filtered for ex3
% end
% ds_factor = 1.5;
% 
% if filtered
%     filename = proc_file(2).name;
% else
%     filename = proc_file(1).name;
% end

onacid_data = load(filename);

cd('C:\Users\intrinsic\Desktop\pyRTAOI-matlab')

%% make data structure
onacid_struct = struct();
onacid_dims = onacid_data.cnm_dims;
onacid_image = reshape(onacid_data.cnm_b,onacid_dims);
onacid_A = full(onacid_data.cnm_A);

if iscell(onacid_A)
    onacid_A = cell2mat(onacid_A);
end

try
    num_frames_onacid = onacid_data.t_cnm;
    full_cnmf = false;
catch
    f_size = size(onacid.cnm_f);
    num_frames_onacid = f_size(2);
    full_cnmf = true;
end

num_comp_onacid = size(onacid_A,2);
comp_shape_onacid = [onacid_dims(1),onacid_dims(2),num_comp_onacid];

try
    ds_factor = onacid_data.ds_factor;
catch
    ds_factor = 1.5;
end

% try
%     cm_onacid = onacid_data.coms*ds_factor;  % center of mass
% catch
%     cm_onacid = com(sparse(double(onacid_A)),onacid_dims(1),onacid_dims(2))*ds_factor;
% end
cm_onacid = com(sparse(double(onacid_A)),onacid_dims(1),onacid_dims(2))*ds_factor;

init = onacid_data.num_frames_init;

for i = 1:num_comp_onacid
    onacid_struct(i).shape = reshape(onacid_A(:,i),onacid_dims);
    onacid_struct(i).noisyC = onacid_data.noisyC(i,1:num_frames_onacid);
    onacid_struct(i).deconvC = onacid_data.cnm_C(i,1:num_frames_onacid);  % deconvC
    onacid_struct(i).centroid = cm_onacid(i,:);
end

%% extract contours
cnm_plot_options = CNMFSetParms;
cnm_plot_options.roi_color = colormap(lines(num_comp_onacid));

CC_onacid = plot_contours(sparse(double(onacid_A)),onacid_image,cnm_plot_options,1,[],[],[1 1 1]);

% scale up contours back to original size
for i=1:num_comp_onacid
    CC_onacid{i} = CC_onacid{i}*ds_factor;
end
close gcf

%% create (binary) cell mask
binary = 0;
separate_cells = 1;

if separate_cells
    mask_onacid = zeros([num_comp_onacid,512,512]);
    for i=1:num_comp_onacid
        mask = reshape(full(onacid_A(:,i)), onacid_dims);
        mask = imresize(mask,ds_factor);
        if binary
             mask(mask>0) = 1;
        end
        mask_onacid(i,:,:) = mask;
    end
else
    mask_onacid = reshape(mean(onacid_A,2), onacid_dims);
    mask_onacid = imresize(mask_onacid,ds_factor);
    if binary
        mask_onacid(mask_onacid>0) = 1;
    end
end

thresh = 2e-4;
binary_mask_onacid = mask_onacid;
binary_mask_onacid(binary_mask_onacid>thresh) = 1;
% binary_mask_onacid(binary_mask_onacid<=thresh) = 0;

%% display mean cell mask
figure,imagesc(squeeze(mean(mask_onacid)))
axis square

%%
%%%%%%%%%%%%%% ANALYSIS OF ALL %%%%%%%%%%%%%%%%%%%%

%% Common contour plot
close gcf
save_folder = strcat(cd, '\results\', example, '\');

% figure;
figH = figure('Position', get(0, 'Screensize'));
set(figH,'color','w');
try
    imagesc(cnmf_image)
catch
    imagesc(imresize(onacid_image,ds_factor))
%     imagesc(mean_img)
end
colormap(gray)
hold on
axis square
text_on = 0;

% suite2p
for i = 1:num_cells_s2p
    off = 15;
    plot(CC_suite2p{i}(1,:),CC_suite2p{i}(2,:),'.r','LineWidth',1)
    
    if text_on
        x = median(CC_suite2p{i}(1,:));
        y = median(CC_suite2p{i}(2,:));
        text(round(x-off),round(y-off),num2str(i),...
        'color','r','fontsize',16,'fontname','helvetica','fontweight','bold')
    end
end

% cnmf
if use_cnmf
    for i=1:num_comp_cnmf
        off = 15;
        plot(CC_cnmf{i}(1,:),CC_cnmf{i}(2,:),'.y','LineWidth',1)

        if text_on
            text(round(cm_cnmf(i,2)+off),round(cm_cnmf(i,1)+off),num2str(i),...
                'color','y','fontsize',16,'fontname','helvetica','fontweight','bold')
        end
    end
end

% onacid
for i=1:num_comp_onacid
    off = 15;
    plot(CC_onacid{i}(1,:),CC_onacid{i}(2,:),'.g','LineWidth',1)
    if text_on
        text(round(cm_onacid(i,2)-off),round(cm_onacid(i,1)+off),num2str(i),...
            'color','g','fontsize',16,'fontname','helvetica','fontweight','bold')
    end
end

if use_cnmf
    title(sprintf('Cells detected\n red: suite2p, yellow: cnmf, green: onacid'))
    img_title = strcat(save_folder, 'plot_all_cells', '.png');
else
    title(sprintf('Cells detected\n red: suite2p, green: onacid'))
    img_title = strcat(save_folder, 'plot_suite2p_onacid', '.png');
end

if save_figs
    % saveas(gcf, strcat(save_folder, 'plot_all_cells', '.png'))
    F = getframe(figH);
    imwrite(F.cdata, img_title, 'png')
end

%% Compare binary masks to determine which cells overlap
c = 1:3;
if ~use_cnmf
    c = 1;
end

%%
for comparison=2 %c

    % suite2p and onacid
    if comparison==1
        compare1 = mask_suite2p;
        compare2 = mask_onacid;
        compare2b = binary_mask_onacid;
        set1 = 'suite2p';
        set2 = 'onacid';
        
    % suite2p and cnmf
    elseif comparison==2
        compare1 = mask_suite2p;
        compare2 = mask_cnmf;
        compare2b = mask_cnmf;
        set1 = 'suite2p';
        set2 = 'cnmf';
        
    % cnmf and onacid
    elseif comparison==3   % code doesn't work well for this case...
        compare1 = mask_cnmf;
        compare1b = binary_mask_cnmf;
        compare2 = mask_onacid;
        compare2b = binary_mask_onacid;
        set1 = 'cnmf';
        set2 = 'onacid';
    end

    k = 1;
    same_cell = {};

    for i=1:size(compare1,1)  % using binary onacid mask for now
        if comparison == 1 | comparison == 2
            cell1 = squeeze(compare1(i,:,:));
        else
            cell1 = squeeze(compare1b(i,:,:));
%             cell1b = squeeze(compare1b(i,:,:));
        end
       

        for j=1:size(compare2b,1)
            cell2 = squeeze(compare2b(j,:,:));
            pixel_overlap = sum(sum(cell1 & cell2));
            
%             pixel_overlap/sum(sum(cell1))
            if pixel_overlap/sum(sum(cell1))>0.7 %min(sum(sum(cell1)),sum(sum(cell2)))>0.8 %pixel_overlap>250
                same_cell{k} = [i,j];
                k = k+1;
            end
        end
    end

    k = size(same_cell,2);
    same_cell_orig = same_cell; % keep orig copy

    %% reset (for testing)
    same_cell = same_cell_orig;
    k = size(same_cell,2);
    
    %% Same cells and others
    idx1 = [];
    idx2 = [];

    for i = 1:k
        idx1(i) = same_cell{i}(1);
        idx2(i) = same_cell{i}(2);
    end

    [u,i] = unique(idx1);
    repeated1 = setdiff(1:length(idx1), i);
    repeated1 = unique(idx1(repeated1));
    
    repeated1_ix = [];
    for i=1:length(repeated1)
        repeated1_ix = [repeated1_ix find(idx1==repeated1(i))];
    end
    
%     sorted_repeats1 = sort(repeats1);
    
    [u,i] = unique(idx2);
    repeated2 = setdiff(1:length(idx2),i);
    repeated2 = unique(idx2(repeated2));
    
    repeated2_ix = [];
    for i=1:length(repeated2)
        repeated2_ix = [repeated2_ix find(idx2==repeated2(i))];
    end
    
    %% repeats: find the one with biggest overlap area and remove rest
    
    for r=1:2
        % reset previous values
        max_ix = [];
        ixs = [];
        curr = 0;
        
        if r==1
            repeated = repeated1;
            repeated_ix = repeated1_ix;
            idx = idx1;
            idx_other = idx2;
        elseif r==2
            repeated = repeated2;
            repeated_ix = repeated2_ix;
            idx = idx2;
            idx_other = idx1;
        end
    
        if length(repeated) > 0
            repeated_grouped = cell(1,length(repeated));
%             others_grouped = cell(1,length(repeated));
            
            j = 1;
            prev = idx(repeated_ix(1));
            for i = 1:length(repeated_ix)
                curr = idx(repeated_ix(i));
                if curr ~= prev
                    j = j+1;
                end
                repeated_grouped{j} = [repeated_grouped{j} repeated_ix(i)];
%                 others_grouped{j} = [others_grouped{j} ];
                prev = curr;
            end
%             for i = repeated_ix
%                 curr = idx(i);
%                 if curr ~= prev
%                     j = j+1;
%                 end
%                 repeated_grouped{j} = [repeated_grouped{j} i];
%                 prev = curr;
%             end

            l = 1;
            for i=1:length(repeated)
                ixs = repeated_grouped{i};
                areas = [];
                for j = ixs %1:length(ixs)
                    cell1 = squeeze(compare1(idx1(j),:,:));
                    cell2 = squeeze(compare2b(idx2(j),:,:));  % use finer non-binary mask?
                    overlap = cell1&cell2;
                    areas(j) = sum(sum(overlap));
                end

                [val,max_ix(l)] = max(areas);
                l = l+1;
            end
        end

        remove_ix = setdiff(repeated_ix,max_ix);
        if r==1
            remove1_ix = remove_ix;
        elseif r==2
            remove2_ix = remove_ix;
        end
    end
    
    %% Quick visual inspection of same cells - individual cells
    inspect = 0;

    if inspect
        figure
        for i=repeated1_ix %1:k
            i
        %     same_cell{i}
        %     imagesc(sum
            cell1 = squeeze(compare1(same_cell{i}(1),:,:));
            cell2 = squeeze(compare2b(same_cell{i}(2),:,:));
            overlap = cell1+cell2;
            imagesc(overlap)
            pause(2)
        end
        close
    end
    
     %% remove: manually or automatically
     manual = 0;
     
     if manual
        remove1_ix = [];
        remove2_ix = [];
     else
         remove_ix = sort([remove1_ix remove2_ix]);
     end
    
    idx1(remove_ix) = [];
    idx2(remove_ix) = [];
    
    ok_size = size(idx1,2) == size(idx2,2);
    if ok_size
        k = size(idx1,2);
    else
        display('Check size of idx1 and idx2!')
    end
    
%     k =  k - length(remove_ix);
    repeated1 = [];
    repeated2 = [];


    %%
    overlap_mask = zeros([512,512]);

    for i = 1:k
        cell1 = squeeze(compare1(idx1(i),:,:)); %(same_cell{i}(1),:,:));
        cell2 = squeeze(compare2(idx2(i),:,:)); %(same_cell{i}(2),:,:));
        if strcmp(set1,'onacid') || strcmp(set1, 'cnmf')
            cell1 = 10*cell1;
        end
        if strcmp(set2,'onacid') || strcmp(set2, 'cnmf')
            cell2 = 10*cell2;
        end
        overlap = cell1+cell2;

        overlap_mask = overlap_mask + overlap;
    end

    all1 = 1:size(compare1,1);
    missing1 = setdiff(all1,idx1);

    all2 = 1:size(compare2,1);
    missing2 = setdiff(all2,idx2);

    thresh=2e-4;
    
    for i = missing1
        miss_cell = squeeze(compare1(i,:,:));
        if strcmp(set1, 'cnmf')
            miss_cell(miss_cell>thresh) = -1-5*miss_cell(miss_cell>thresh);
        else
            miss_cell(miss_cell==1) = -1;
        end
        overlap_mask = overlap_mask + miss_cell;
    end

    for i = missing2
        miss_cell = squeeze(compare2(i,:,:));
        miss_cell(miss_cell>thresh) = -3-1*miss_cell(miss_cell>thresh);  % == 1

        overlap_mask = overlap_mask + miss_cell;
    end
    
%     figure,
    figH = figure('Position', get(0, 'Screensize'));
    set(figH,'color','w');
    imagesc(overlap_mask)
    axis square
    colorbar
    title(sprintf('Same cells detected for %s and %s (yellow)\nlight blue: %s only; dark blue: %s only', ...
        set1,set2,set1,set2))
    single_plot = strcat(save_folder,set1,'_',set2,'_same_cells.png');
    
    if save_figs
        F = getframe(figH);
        imwrite(F.cdata, single_plot, 'png')
    %     saveas(gcf, plot_title, 'png')
    end

    %% Compare traces of same cells
    single_plot = 1;
%     figH = figure('Position', get(0, 'Screensize'));

    for x=1:2
        if single_plot
            if x==1
                figH = figure('Position', get(0, 'Screensize'));
            end
            subplot(num2str(strcat(num2str(12),num2str(x))))
        else
            figH = figure('Position', get(0, 'Screensize'));
        end
%         subplot(num2str(strcat(num2str(12),num2str(x))))
        hold on
        set(figH,'color','w');

        compare = [];

        if x==1
            compare = compare1;
            idx = idx1;
        elseif x==2
            compare = compare2;
            idx = idx2;
        end
        

        if size(compare) == size(mask_suite2p)
            if compare == mask_suite2p
                plot_offset = 80000;

                j = 1;
                for i = idx
                    if find(repeated1==i)
                        col = 'red';
                    else
                        col = 'black';
                    end
                %     plot(s2p_struct(i).Fcell+i*plot_offset,'color',roi_color(i,:))
                    plot(s2p_struct(i).Fsub+j*plot_offset,'color',col,'linewidth',1.5)
                    j = j+1;
                end
                title('suite2p')
                xlim([0, size(s2p_struct(i).Fsub,2)])
            end

        elseif use_cnmf & size(compare) == size(mask_cnmf)
            if compare == mask_cnmf
                plot_offset = 80000;

                j = 1;
                for i = idx
                %     plot(cnmf_struct(i).noisyC+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
                    plot(cnmf_struct(i).deconvC+j*plot_offset,'color','black','linewidth',1.5)
                    j = j+1;
                end
                title('cnmf')
                xlim([0, size(cnmf_struct(i).deconvC,2)])
            end

        elseif size(compare) == size(mask_onacid)
            if compare == mask_onacid
                plot_offset = 40; % 40 or 45

                j = 1;
                for i = idx
                    if find(repeated2==i)
                        col = 'red';
                    else
                        col = 'black';
                    end
                %     plot(cnm_struct(i).noisyC+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
                    plot(onacid_struct(i).deconvC+j*plot_offset,'color',col,'linewidth',1.5)
                    j = j+1;
                end
                title('onacid')
                xlim([0, size(onacid_struct(i).deconvC,2)])
            end
        end

        xlabel('Frames')
        ylabel('ROI')

        yt = plot_offset*(size(idx,2)+2); %yticks;
        ylim([0,yt])
        
        if x==1
            y = yticks;
            common_yticks = int64(y/plot_offset);
            yticks(common_yticks*plot_offset);
            yticklabels({common_yticks});
        elseif x==2
            yticks(common_yticks*plot_offset);
            yticklabels({yticks/plot_offset});
        end
        
        yticklabels(common_yticks)
        
    end
    
    suptitle('Matched traces of overlapping cells')
    
    if save_figs
        F = getframe(figH);
        if single_plot
            same_plots = strcat(save_folder,set1,'_',set2,'_same_cells_plots.png');
        else
            if x==1
                set = set1;
                set_other = set2;
            elseif x==2
                set = set2;
                set_other = set1;
            end
            same_plots = strcat(save_folder,set,'_same_cells_plots',set_other,'.png');
        end
        imwrite(F.cdata, same_plots, 'png')
    end

    
    %% Compare traces of non-overlapping cells
    single_plot = 1;
    
%     figH = figure('Position', get(0, 'Screensize'));
    
    for x=1:2
        if single_plot
            if x==1
                figH = figure('Position', get(0, 'Screensize'));
            end
            subplot(num2str(strcat(num2str(12),num2str(x))))
        else
            figH = figure('Position', get(0, 'Screensize'));
        end
        
        set(figH,'color','w');
        hold on
        compare = [];

        if x==1
            compare = compare1;
            idx = missing1;
        elseif x==2
            compare = compare2;
            idx = missing2;
        end

        if size(compare) == size(mask_suite2p)
            if compare == mask_suite2p
                plot_offset = 80000;

                j = 1;
                for i = idx
                %     plot(s2p_struct(i).Fcell+i*plot_offset,'color',roi_color(i,:))
                    plot(s2p_struct(i).Fsub+j*plot_offset,'color','black','linewidth',1.5)
                    j = j+1;
                end
                title('suite2p')
                xlim([0, size(s2p_struct(i).Fsub,2)])
            end
            
        elseif use_cnmf & size(compare) == size(mask_cnmf)
            if compare == mask_cnmf
                plot_offset = 80000;

                j = 1;
                for i = idx
                %     plot(cnmf_struct(i).noisyC+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
                    plot(cnmf_struct(i).deconvC+j*plot_offset,'color','black','linewidth',1.5)
                    j = j+1;
                end
                title('cnmf')
                xlim([0, size(cnmf_struct(i).deconvC,2)])
            end

        elseif size(compare) == size(mask_onacid)
            if compare == mask_onacid
                plot_offset = 40; % 40 or 45

                j = 1;
                for i = idx
                %     plot(cnm_struct(i).noisyC+i*plot_offset,'color',cnm_plot_options.roi_color(i,:))
                    plot(onacid_struct(i).deconvC+j*plot_offset,'color','black','linewidth',1.5)
                    j = j+1;
                end
                title('onacid')
                xlim([0, size(onacid_struct(i).deconvC,2)])
            end
        end

        xlabel('Frames')
        ylabel('ROI')
        suptitle('Traces of non-overlapping cells')

        yt = plot_offset*(size(idx,2)+2);
        ylim([0,yt])
        
        yticks_now = yticks;
        yticks_new = unique(int64(yticks/plot_offset))*plot_offset;
        yticks(yticks_new);
        yticklabels({yticks_new/plot_offset})
    end
    
    if save_figs
        F = getframe(figH);
        if single_plot
            diff_plots = strcat(save_folder,set1,'_',set2,'_different_cells_plots.png');
        else
            if x==1
                set_ = set1;
            elseif x==2
                set_ = set2;
            end
            diff_plots = strcat(save_folder,set_,'_only_plots.png');
        end
        imwrite(F.cdata, diff_plots, 'png')
    end
    
    %% Summary
%     clc

    % save summary as .txt file
    if save_figs
        file = strcat(save_folder, strcat(set1,'_',set2,'_summary.txt'));
    else
        file = 'temp.txt';
    end
    
    fid = fopen(file,'wt');
    fprintf(fid,'Total number of cells in %s set: %d', set1, size(all1,2));
    fprintf(fid,'\nTotal number of cells in %s set: %d', set2, size(all2,2));
    fprintf(fid,'\nTotal number of cells in both sets: %d', k);
    % keep just to double check removing repeats worked
    fprintf(fid,'\nNumber of unique repeats in %s set: %d', set1, size(unique(idx1),2));
    fprintf(fid,'\nNumber of unique repeats in %s set: %d', set2, size(unique(idx2),2));
    fclose(fid);
    
    % read saved file
    fid = fopen(file,'r');
    fscanf(fid,'%c')
    fclose(fid);

end

%% COMPARE ONACID RESULTS WITH 'GROUND TRUTH' OF SUITE2P
% 1. the cells that are detected by suite2p but missed by onACID are 'quiet
% cells' 
% ******* DONE **********
% --> non-overapping cells plot in separate mode (single_plot = 0)
% --> quantify total quietness

% 2. the cells detected by onACID are 'quiet' before detection
% ******* DONE **********
% --> plot of mean s2p spike rate before detection

% 3. the cells detected by onACID are actually cells
% ******* DONE? **********
% --> contour plot of cells
% --> onacid traces (separate mode available)

%% Quantify quietness
%  dat.stat.st, and dat.stat.c for spike times and their amplitude

% s2p_spike_trace = zeros(1, total_number_frames);
% s2p_spike_trace(spike_frames) = spike_amplitude;
% duration = length(frame_range)/frame_rate;
% s2p_spike_rate = sum(s2p_spike_trace(frame_range))/duration;

frame_rate = 30; % Hz
tot_frames = size(s2p_struct(1).Fcell,2);
tot_duration = tot_frames/frame_rate; % in s
    
s2p_spike_trace = zeros(num_cells_s2p, tot_frames);

% total spike rate
tot_frame_range = 1:tot_frames;
tot_s2p_spike_rate = zeros(num_cells_s2p,1);

% get deconvolved spikes
for i=1:num_cells_s2p
    spike_frames = s2p_data.stat(i).st;
    spike_amplitude = s2p_data.stat(i).c;

    s2p_spike_trace(i,spike_frames) = spike_amplitude;
    tot_s2p_spike_rate(i) = sum(s2p_spike_trace(i,:))/tot_duration; % total spike rate
end

%% Plot total quietness for all cells
comparison = 1;
compare1 = mask_suite2p;
compare2 = mask_onacid;
compare2b = binary_mask_onacid;
set1 = 'suite2p';
set2 = 'onacid';

if comparison == 1                               % 1 = suite2p; 2 = onacid

    SR_det = tot_s2p_spike_rate(unique(idx1));
    SR_undet = tot_s2p_spike_rate(unique(missing1));
    mean_SR_det = mean(SR_det);
    std_SR_det = std(SR_det);
    mean_SR_undet = mean(SR_undet);
    std_SR_undet = std(SR_undet);
    
    % sort all SR data
    [sorted_SR,orig_idx] = sort(tot_s2p_spike_rate,'ascend');
    
    new_idx1 = [];
    for i=unique(idx1)
        new_idx1 = [new_idx1 find(orig_idx==i)];
    end
    new_idx2 = [];
    for i=unique(missing1)
        new_idx2 = [new_idx2 find(orig_idx==i)];
    end

    sorted = 1;
    figH = figure('Position', [1920/2-700/2,0,700,1080-100]);
    set(figH,'color','w');
    hold on
    if sorted
        p1 = barh(new_idx1,SR_det);
        p2 = barh(new_idx2,SR_undet);
%         p0 = bar(sorted_SR);
%         p1 = bar(1:length(unique(idx1)),sort(SR_det,'descend'));
%         p2 = bar(1:length(unique(missing1)),sort(SR_undet,'descend'));
    else
        p1 = barh(unique(idx1),SR_det);
        p2 = barh(unique(missing1),SR_undet);
    end

    xlim([0,ceil(max(tot_s2p_spike_rate))])
    ylim([0,num_cells_s2p+1])
    
    set(p1,'FaceColor','blue');
    set(p2,'FaceColor','red');
    lgd = legend('OnACID detected','OnACID not detected','Location','southeast');
%     lgd.FontSize = 13;
   
    if sorted
        ylabel('Sorted cells','FontSize',13)
    else
        ylabel('Cell number','FontSize',13)
    end
    xlabel('Spikes per second (sps)','FontSize',13)
    
    title(sprintf(['Mean spike rate of %d suite2p cells\n' ...
        'Mean SR detected: %1.2f',char(177),'%1.2f sps; '...
        'mean SR not detected: %1.2f', char(177),'%1.2f sps'], ...
        num_cells_s2p, mean_SR_det,std_SR_det,mean_SR_undet,std_SR_undet),...
        'FontSize',13)

    if save_figs
        F = getframe(figH);
        tot_sr_bar = strcat(save_folder,'tot_sr_quietness.png');
        imwrite(F.cdata, tot_sr_bar, 'png')
    end
    
%% Check quitness of cells before OnACID detection
    onacid_cells = length(unique(idx2));
    frame_det = zeros(onacid_cells,1); % when onacid detects the cell
    s2p_sr_undet = zeros(size(frame_det));
 
    j = 1;
    
    for i=unique(idx2)
        frame_det(j) = find(onacid_struct(i).deconvC>0,1);
        undetected = 1:frame_det(j);
        duration = length(undetected);
        s2p_sr_undet(j) = sum(s2p_spike_trace(undetected))/duration;
        
        j = j+1;
    end
    
    mean_sr_undet = mean(s2p_sr_undet);
    std_sr_undet = std(s2p_sr_undet);
    
    % set(figH, 'Position', [x y width height])
    figH = figure('Position', [1920/2-700/2,0,700,1080-100]);
    set(figH,'color','w');
    barh(sort(s2p_sr_undet,'ascend'))
    
    xlim([0, ceil(max(s2p_sr_undet))])
    ylim([0, onacid_cells+1])
    
    xlabel('Spikes per second (sps)','FontSize',13)
    ylabel('Sorted cells','FontSize',13)
    
    title(sprintf(['Spike rate of %d common cells before OnACID detection\n',...
        'Mean spike rate: %1.2f',char(177),'%1.2f sps'], ...
        onacid_cells,mean_sr_undet,std_sr_undet),'FontSize',13)
    
    if save_figs
        F = getframe(figH);
        before_det_bar = strcat(save_folder,'before_det_sr_quietness.png');
        imwrite(F.cdata, before_det_bar, 'png')
    end

else
        display('Wrong comparison selected above!')
end
