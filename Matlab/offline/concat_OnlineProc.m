function [cell_struct,fov_struct] = concat_OnlineProc(data_path,num_init_comp)
% concatenate traces in pyROAOI result structures
% initialisation traces are cropped out
% (initialised with the same file)
% get files from folder
% they use the same initialisation file


%% find files
file_list = dir(data_path);
file_list = file_list(cellfun(@(x)contains(x,'.mat'),{file_list.name}));
file_names = {file_list.name};
num_files = numel(file_names);
disp(['Number files: ' num2str(num_files)])


%% making a suite2p-trace-like struct
Fcaiman = {};
%% loop through files
cell_struct = struct();
for f = 1:num_files
    caiman_data = load([data_path filesep file_names{f}]);
    num_frames = caiman_data.t_cnm;
    disp(['File ' num2str(f) ',number frames: ' num2str(num_frames)])

    %% get spatial footprints from the first file
    if f == 1
        % get ROIs spatial footprint in initialisation movie
        num_frames_init = caiman_data.num_frames_init;
        cnm_dims = caiman_data.cnm_dims;
        cnm_image = reshape(caiman_data.cnm_b,cnm_dims);
        cnm_A = full(caiman_data.cnm_A);
        
        num_comp = size(cnm_A,2); % all detected components
        comp_shape =[cnm_dims(1),cnm_dims(2),num_comp];
        cm = com(sparse(double(cnm_A)),cnm_dims(1),cnm_dims(2)); % center of mass
        
        num_init_comp = caiman_data.init_com_count;
        for i = 1:num_init_comp
            cell_struct(i).shape = reshape(cnm_A(:,i),cnm_dims);
            cell_struct(i).centroid = cm(i,:);
            cell_struct(i).deconvC_full =  caiman_data.cnm_C(i,num_frames_init+1:num_frames);
            cell_struct(i).noisyC_full = caiman_data.noisyC(i,num_frames_init+1:num_frames);
        end
        
        %% opsin expression
        opsin_positive = caiman_data.opsin_positive;
        accepted_idx = caiman_data.accepted_idx+1;
        opsin_positive_idx = accepted_idx(opsin_positive>0);
        
        %% plot spatial components
        com_fov = zeros(cnm_dims);
        binary_fov = zeros(cnm_dims);
        for i = 1:num_init_comp
            com_fov = com_fov+cell_struct(i).shape;
        end
        
        cnm_plot_options = CNMFSetParms;
        cnm_plot_options.roi_color = [colormap(lines);colormap(lines);colormap(lines)];
        
        figure('name','fov','position',[100 100 1200 800])
        subplot(1,2,1)
        imagesc(com_fov)
        colormap(gray)
        axis square
        
        subplot(1,2,2)
        [CC,jsf] = plot_contours(sparse(double(cnm_A)),cnm_image,cnm_plot_options,1,[],[],[1 1 1]);
        
        fov_struct.cnm_A = cnm_A;
        fov_struct.cnm_plot_options = cnm_plot_options;
        fov_struct.com_fov = com_fov;
        fov_struct.cnm_image = cnm_image;
        
        for i = 1:num_init_comp
            temp_coords = jsf(i).coordinates;
            lin_idx = zeros(size(temp_coords,1),1);
            
            for t = 1:size(temp_coords,1)
                lin_idx(t) = sub2ind(cnm_dims,temp_coords(t,1),temp_coords(t,2));
            end
            
            cell_struct(i).contour = CC{i};
            cell_struct(i).lin_coords = lin_idx;
            cell_struct(i).coordinates = jsf(i).coordinates;
            cell_struct(i).pix_values = jsf(i).values;
            cell_struct(i).centroid = jsf(i).centroid;
            cell_struct(i).opsin_positive = 0;
            if(~isempty(find(opsin_positive_idx==i)))
                cell_struct(i).opsin_positive = 1;
            end
        end
    else
        for i = 1:num_init_comp
            cell_struct(i).deconvC_full = [ cell_struct(i).deconvC_full caiman_data.cnm_C(i,num_frames_init+1:num_frames)];
            cell_struct(i).noisyC_full = [ cell_struct(i).noisyC_full caiman_data.noisyC(i,num_frames_init+1:num_frames)];
        end
        
    end
    Fcaiman{f} =  caiman_data.cnm_C(1:num_init_comp,num_frames_init+1:num_frames);

end
%% save suite2p suite2p-trace-like struct
% time = datestr(now,'yyyymmdd_HHMM');
% save([save_path filesep  'Fcaiman_' data_name ],'Fcaiman' )

end

