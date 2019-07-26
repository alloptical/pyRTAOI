function [] = mark_cells_on_fov(cell_struct, mark_idx,color, varargin)
MarkerSize = 140;
Linewidth = 1.5;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'MarkerSize')
        MarkerSize = varargin{v+1};
    end
    if strcmpi(varargin{v},'Linewidth')
        Linewidth = varargin{v+1};
    end
end

for i = 1:length(mark_idx)
    scatter(round(cell_struct(mark_idx(i)).centroid(:,2)),round(cell_struct(mark_idx(i)).centroid(:,1)),...
        MarkerSize,'MarkerEdgeColor',color,'Linewidth',Linewidth)
end

