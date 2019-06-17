function [output_img] = make_centroid_image(centroids,fov_size,ds_factor)
% hard coded to be 512x512 
% if movie were downsampled then need to transform back
output_img = zeros(512,512);
center_pix = fov_size(1)/ds_factor;
ROIy = round(ds_factor.*(centroids(:,1)-center_pix)+256);
ROIx = round(ds_factor.*(centroids(:,2)-center_pix)+256);

for c = 1:size(centroids,1)
    output_img(ROIy(c),ROIx(c)) = 1;
end

% figure;imshow(output_img)

end

