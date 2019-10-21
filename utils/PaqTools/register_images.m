function varargout=register_images(basename,unregisteredname,maxpixelval,useCpcorr,varargin)
% Adam Packer  July 18th, 2008
% Quick image registration wrapper/script
% Inputs
% 1) base image
% 2) unregistered image
% 3) max pixel value, higher=brighter, 255 is default even histogram...
% 4) switch to use cpcorr to refine points
% 5) varargin use initial points if desired
% Outputs
% 1) handle to base image to change transparency
% 2) tform
% 3) registered image

% read images
base=imread(basename);
unregistered=imread(unregisteredname);

base8=base-abs(min(min(base)));
base8=uint8(base8*(maxpixelval/double(max(max(base8)))));
unregistered8=unregistered+abs(min(min(unregistered)));
unregistered8=uint8(unregistered8*(maxpixelval/double(max(max(unregistered8)))));

% run cpselect tool to match up image points
% make sure to export points to workspace!
if ~isempty(varargin)
    moving_in=varargin{1};
    fixed_in=varargin{2};
    [moving_out,fixed_out]=cpselect(unregistered8,base8,moving_in,fixed_in,'wait',true);
else
    [moving_out,fixed_out]=cpselect(unregistered8,base8,'wait',true);
end
% input_points = evalin('base','input_points');
% base_points = evalin('base','base_points');

if useCpcorr
    moving_out_corr = cpcorr(moving_out,fixed_out,unregistered,base);
end

% calculate affine image transformation using matched image points
if useCpcorr
    tform = fitgeotrans(moving_out_corr,fixed_out,'affine'); % AP20140120
else
    tform = fitgeotrans(moving_out,fixed_out,'affine'); % AP20140120
end
% tform = cp2tform(moving_out, fixed_out, 'affine');

% create registered image that fits in the size of the base image
r=imref2d([size(base,1) size(base,2)]);
registered = imwarp(unregistered8,tform,'FillValues',0,'outputview',r);
% registered = imtransform(unregistered,tform,...
% 'FillValues', 255,...
% 'XData', [1 size(base,2)],...
% 'YData', [1 size(base,1)]);

% set up the figure
figure;
imshow(base8,[]);
hold on
h = imshow(registered,[]);
set(h,'AlphaData',0.5);

% set outputs
varargout{1}=h;
varargout{2}=tform;
varargout{3}=registered;
varargout{4}=moving_out;
varargout{5}=fixed_out;
if useCpcorr
    varargout{6}=moving_out_corr;
end