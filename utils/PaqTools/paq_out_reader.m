%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose and open file
function [data]=paq_out_reader;
[filename,pathname,FilterIndex]=uigetfile({'*.txt;*.bin'},'Choose a data file');
if ~FilterIndex
    return
end
fullpath=[pathname filename];
fid=fopen(fullpath);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in rate, number of channels, and channel names
DeviceChars=fread(fid,1,'float32','b');
Device{1}=[];
for i=1:DeviceChars;
        Device{i}=[Device{i-1}, strrep(fread(fid,1,'float32=>char', 'b'),' ','')];
end
data=Device;
fclose(fid);
