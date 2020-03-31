function [ img,zlimit ] = plot_value_in_rois( cell_struct, value_field,dims,ax,varargin)
% pyrtaoi post analysis
IF_NORM_PIX = 0;
IF_CONTOUR = 1;
textcolor = [.3 .3 .3];

colorlut = [];
zlimit = []; % dummy dots to force colorlut matching preset max and min range
show_cell_idx = [];
target_cell_idx = [];
trigger_cell_idx = [];
max_abs_value = [];
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'IF_NORM_PIX')
        IF_NORM_PIX = varargin{v+1};
    elseif strcmpi(varargin(v),'colorlut')
        colorlut = varargin{v+1};
    elseif  strcmpi(varargin(v),'IF_CONTOUR') % THIS IS SLOW
        IF_CONTOUR = varargin{v+1};
    elseif  strcmpi(varargin(v),'zlimit')
        zlimit = varargin{v+1};
    elseif  strcmpi(varargin(v),'show_cell_idx')
        show_cell_idx = varargin{v+1};
    elseif  strcmpi(varargin(v),'target_cell_idx')
        target_cell_idx = varargin{v+1};
    elseif  strcmpi(varargin(v),'trigger_cell_idx')
        trigger_cell_idx = varargin{v+1};
    elseif  strcmpi(varargin(v),'max_abs_value')
        max_abs_value = varargin{v+1};
    end
end

img = zeros(dims);
all_values = extractfield(cell_struct,value_field);

for i = 1:size(cell_struct,2)
    xy_coords = sub2ind(dims,cell_struct(i).coordinates(:,1),cell_struct(i).coordinates(:,2));
    plot_value = cell_struct(i).(value_field);
    if isnan(plot_value)
        continue
    end
    
    if IF_NORM_PIX
        pix_values = cell_struct(i).pix_values;
        norm_pix_values = pix_values;%./sum(pix_values);
        img(xy_coords) = norm_pix_values.*plot_value;
    else
        img(xy_coords) =plot_value;
    end
%     text(round(cell_struct(i).centroid(:,2)),...
%         round(cell_struct(i).centroid(:,1)),num2str(i),'color',textcolor,'fontweight','bold');
    
end
all_values = img(:);
if isempty(zlimit)
    zlimit = [min(all_values),max(all_values)];
    if zlimit(1)>0
       zlimit(1) = 0; % forcing min color to be white if all values are positive 
    end
    if ~isempty(max_abs_value)
        zlimit(end) = min(zlimit(end),max_abs_value);
        zlimit(1) = max(zlimit(1),-max_abs_value);
    end
else
    img(1) = zlimit(1);
    img(2) = zlimit(2);
end

hold on
imagesc(img);
if isempty(colorlut)
    colormap(ax,b2r(zlimit(1) ,zlimit(2)))
else
    colormap(ax,colorlut)
end


if IF_CONTOUR
    try
        cc_pix = cell2mat({cell_struct(target_cell_idx).contour});
        scatter(cc_pix(1,:),cc_pix(2,:),1,'r')
        cc_pix = cell2mat({cell_struct(setdiff( 1:size(cell_struct,2),target_cell_idx)).contour});
        scatter(cc_pix(1,:),cc_pix(2,:),1,[.5 .5 .5])
        IF_CONTOUR = 0;
    catch
        warning('getting contours')
    end
else
%     scatter(round(cell_struct(i).centroid(:,2)),round(cell_struct(i).centroid(:,1)),...
%         10,textcolor)
end

% mark cell indices
if ~isempty(show_cell_idx)
    plot_tex_idx = show_cell_idx;
else
    plot_tex_idx = 1:size(cell_struct,2);
end
for i = plot_tex_idx
    text(round(cell_struct(i).centroid(:,2)),...
        round(cell_struct(i).centroid(:,1)),num2str(cell_struct(i).cnm_idx),'color',textcolor,'fontweight','bold');
end

% mark target cells
if ~isempty(target_cell_idx)
    for i = target_cell_idx
     text(round(cell_struct(i).centroid(:,2)),...
        round(cell_struct(i).centroid(:,1)),num2str(cell_struct(i).cnm_idx),'color','r','fontweight','bold');
    end
    cc_pix = cell2mat({cell_struct(target_cell_idx).contour});
    scatter(cc_pix(1,:),cc_pix(2,:),1,'r')
end

% mark trigger cells
if ~isempty(trigger_cell_idx)
    for i = trigger_cell_idx
     text(round(cell_struct(i).centroid(:,2)),...
        round(cell_struct(i).centroid(:,1)),num2str(cell_struct(i).cnm_idx),'color','blue','fontweight','bold');
    end
    cc_pix = cell2mat({cell_struct(trigger_cell_idx).contour});
    scatter(cc_pix(1,:),cc_pix(2,:),1,[.5 .5 .5])

end

xlim([1,size(img,1)])
ylim([1,size(img,2)])

cb = colorbar('location','southoutside');
% cb.TickLabels = {num2str(zlimit(1),'%.2f'), num2str(zlimit(2),'%.2f')};

box off
axis off
axis square
set(gca,'YDir','reverse')


end

