function [norm_weights,norm_thresh] = get_norm_weights(weights,thresh,cell_mean,cell_std)
% use weights and thresh to normalise data
norm_weights = weights./cell_std';
norm_thresh = thresh + cell_mean./cell_std*weights;
end

