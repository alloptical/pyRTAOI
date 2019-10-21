%Assumes that the current directory is a folder where all .fig files are
%figures created by the triggered averaging functions of EphysViewer.  This
%function collects data from each .fig file and stores it to combine all of
%it into a new figure with the same format as the triggered averaging
%figures.  Individual triggered traces (not the averages) are collected and
%averaged to make averages of larger numbers of figures than in the
%originals.  
% All of the original averages are assumed to have the same number of
% points before and after the triggers.



d = dir;
d = d(3:end);
trigtraces = [];
rawtraces = [];
zeroedtraces = [];
outputfig = [];
for fidx = 1:length(d);
    if ~(d(fidx).isdir)
%copy figure.  Change appdata, plots put info in titles
        if strcmp(d(fidx).name(end-3:end),'.fig')
            h = open(d(fidx).name);
            dat = getappdata(h,'data');
            if ~isstruct(dat);%making up for an earlier error
                x = dat;
                clear dat
%                 dat.data = x;
                dat.data{1} = x{1};
                dat.data{3} = x{2};
                dat.data{5} = x{4};
            end
            trigtraces = cat(1,trigtraces,dat.data{1});
            rawtraces = cat(1,rawtraces,dat.data{3});
            zeroedtraces = cat(1,zeroedtraces,dat.data{5});
            close(h);
        end
    end
end

outputfig = TriggeredAverageKernelFcn(trigtraces, rawtraces, zeroedtraces);

foldername = cd;
slashes = strfind(foldername,'\');
set(outputfig,'name',['Averages from ',foldername(slashes(end)+1:end)]);