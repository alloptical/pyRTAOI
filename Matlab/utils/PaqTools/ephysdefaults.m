function [rangequest,qmquestion,chunkquest,channels,filter,stackquest,laserquest]=ephysdefaults
% Adam Packer
% April 10th, 2007
% Added filter capability December 18th, 2009
% Added laserquest capability June 29th, 2010

% Set defaults for EphysViewer
% rangequest=1 lets you decide what range of data you want to import
% qmquestion=1 lets you import the data in a slower but more memory efficient manner
% chunkquest=1 lets you split the data into chunks based on some trigger
% channels=1 lets you chose only certain channels to import
% filter=x runs medfilt1 on the data with order x
% stackquest=1 asks you if you want to split the pulses to different axes (deprecated)
% laserquest=x displays a laser power rounded to tens of percentage values
% where x volts is 100%
rangequest=0;
qmquestion=0;
chunkquest=0;
channels=1;
filter=0;
stackquest=0;
laserquest=0.79;