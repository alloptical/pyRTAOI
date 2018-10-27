Thiszoom     = 2;
Refzoom      = 1.14;


[file_name,path_name] = uigetfile('*.bmp','Select centroid ROI image(s)','MultiSelect','on');
[a, image_name] = fileparts(file_name);
if path_name == 0
    warning('Cancelled')
    return
end
file_path = [path_name filesep file_name];
cd(path_name);
AllROIimage = imread(file_name);

[ROIy, ROIx] = find(AllROIimage);
numROI = length(ROIy);
% make zoomed-in ROIs

ROIy = round(Thiszoom/Refzoom.*(ROIy-256)+256);
ROIx = round(Thiszoom/Refzoom.*(ROIx-256)+256);
zoom_centroid_img = zeros (size(AllROIimage));
for i = 1:numROI
    if(ROIy(i)>0 && ROIx(i)>0 && ROIy(i)<=512 && ROIx(i)<=512)
        zoom_centroid_img(ROIy(i),ROIx(i)) = 1;
    end
    
end
figure;imshow(zoom_centroid_img); title(['zoomed targets ' num2str(Thiszoom) 'X']);
imwrite(zoom_centroid_img,[path_name 'zoom' num2str(Thiszoom) file_name]);