function EphysViewer(varargin)
%this function will read a .daq file and will plot the data in a way
%similar to RecordDaq, but will plot each channel on a separate axes.  It
%will eventually hopefully have some nice data output options.

warning off

handles.filetypes = {'*.daq;*.bin;*.paq;*.fig'};
handles.position = [0.2266    0.3359    0.5469    0.5469];

if nargin==1;%assume input is the path to the file with the data to view
    handles.path=varargin{1};%save pathname into handles
elseif nargin==2;
    handles.path=varargin{1};
    handles.rechunk=varargin{2};
elseif nargin==0;%if no inputs
    [filename,pathname,FilterIndex]=uigetfile(handles.filetypes,'Choose a data file');
    if ~FilterIndex
        return
    end
    handles.path=[pathname,filename];
end

ViewerMainFcn(handles);%open figure with all necessary axes etc.

function varargout = ViewerMainFcn(handles)
%function ViewerMainFcn(handles)
%handles must contain .path, containing the path to the file to be opened

warning off
[pathstr, name, ext] = fileparts(handles.path);
handles.filename = name;
if strcmp(ext,'.daq')
    info=bdaqread(handles.path,'info');
elseif strcmp(ext,'.bin');
    info=paq2lab(handles.path,'info');
elseif strcmp(ext,'.paq');
    info=paq2lab(handles.path,'info');
elseif strcmp(ext,'.fig');%for opening Figures representing Triggered Averages from EphysViewer
    ReadFig = hgload(handles.path);
    set(ReadFig,'visible','off');
    try
        ReadFigDat = getappdata(ReadFig,'data');
    catch
        errdlg ('Invalid Figure.  Must have "data" in appdata containing a structure of line data')
        return
    end
    close(ReadFig)
    clear ReadFig
    info.ObjInfo.SampleRate = 1;%assume nothing
    info.HwInfo = 'From Triggered Average Figure';
    info.ObjInfo.SamplesAquired = size(ReadFigDat.data{1},2);
    AveragedChannels = [2 4 6];%these are only channels with single datastrings
    for a = 1:length(AveragedChannels);%for each averaged channel
        thischannel = AveragedChannels(a);
        info.ObjInfo.Channel(a).HwChannel = thischannel;
        info.ObjInfo.Channel(a).Units = 'Units';
        info.ObjInfo.Channel(a).ChannelName = ReadFigDat.names{thischannel};
    end
end
%Viewer expects the following information in info
% info.ObjInfo.SampleRate = Data Samples (Points) Per Second;
% info.HwInfo = Name of Acquistion system?;
% info.ObjInfo.SamplesAcquired = number of datapoints per channel;
%
% for each channel (channel = a)
% info.ObjInfo.Channel(a).HwChannel = Channel Number in the acquisition hardware;
% info.ObjInfo.Channel(a).Units = ie mV, nA, Volts, etc;
% info.ObjInfo.Channel(a).ChannelName = user-chosen name for each channel;
for a=1:length(info.ObjInfo.Channel)%for each channel
    if strcmp(ext,'.paq');
        filechans{a}=info.ObjInfo.Channel(a).HwChannel;
        nam=info.ObjInfo.Channel(a).ChannelName;
        charcell{a}=strcat(filechans{a},'-',nam);
    else
        filechans(a)=info.ObjInfo.Channel(a).HwChannel;
        nam=info.ObjInfo.Channel(a).ChannelName;
        charcell{a}=[num2str(filechans(a)),' - ',nam];
    end
end
% filechans=sort(filechans);

% Get ephys defaults
try
    [rangequest,qmquestion,chunkquest,channels,filter,stackquest,laserquest]=ephysdefaults;
catch
    chunkquest=0;
    channels=1;
    laserquest=0;
end
handles.chunkquest=chunkquest;
% Find out which channels to plot
if channels==1;
    [handles.includedchannels,ok]=listdlg('ListString',charcell,'PromptString','Plot which channels... All?',...
        'Name','Channel Selection','InitialValue',1:length(charcell));
    if isempty(ok) || ~ok;
        return
    end
elseif channels==0;
    numberofchannels=length(charcell);
    handles.includedchannels=1:numberofchannels;
end

% Open the figure
handles.displayfig=figure('Visible','off',...
    'Units','Normalized',...
    'Position',handles.position,...
    'CloseRequestFcn',@CloseViewerFcn,...
    'numbertitle','off',...
    'tag','EphysViewer',...
    'name',handles.filename);%name after the file
%,...%     'menubar','none',...
%     'ResizeFcn',@ViewerResizing);%all may help faster rendering
%     'DoubleBuffer','on',...
%     'BackingStore','off',...
%     'DithermapMode','manual',...
%     'RendererMode','manual',...
%     'Renderer','OpenGL',...

% Set up the figure menus
handles.MainFigDataMenu=uimenu('Label','&Data');
uimenu(handles.MainFigDataMenu,'Label','Load data','Callback',@LocalMenuLoadDataFcn);
try
    if ismac
        uimenu(handles.MainFigDataMenu,'Label','Load next file in folder','Callback',@LocalMenuLoadNextFileFcn,'Accelerator','g');
    else
        uimenu(handles.MainFigDataMenu,'Label','Load next file in folder','Callback',@LocalMenuLoadNextFileFcn,'Accelerator','m');
    end
catch
    uimenu(handles.MainFigDataMenu,'Label','Load next file in folder','Callback',@LocalMenuLoadNextFileFcn,'Accelerator','m');
end
uimenu(handles.MainFigDataMenu,'Label','Load previous file in folder','Callback',@LocalMenuLoadPrevFileFcn,'Accelerator','f');
uimenu(handles.MainFigDataMenu,'Label','Load fluorescent traces','Callback',@LocalMenuLoadFluorescentTracesFcn);
uimenu(handles.MainFigDataMenu,'Label','Load synchronous trace','Callback',@LocalMenuLoadSynchronousTraceFcn);
uimenu(handles.MainFigDataMenu,'Label','ReChunkerize data','Callback',@LocalMenuReChunkerizeDataFcn);
uimenu(handles.MainFigDataMenu,'Label','Export visible data','Callback',@LocalMenuExportVisibleDataFcn);

handles.MainFigDisplayMenu=uimenu('Label','Display');
uimenu(handles.MainFigDisplayMenu,'Label','Grid','Callback',@LocalMenuGridFcn);
uimenu(handles.MainFigDisplayMenu,'Label','Load associated image in separate window','Callback',@LoadImageFcn);
uimenu(handles.MainFigDisplayMenu,'Label','Plot entire VNT on image','Callback',@PlotEntireVNTOnImage);
uimenu(handles.MainFigDisplayMenu,'Label','Load contours','Callback',@LoadContoursFcn);

handles.MainFigDisplayMenu=uimenu('Label','&Analyze');
uimenu(handles.MainFigDisplayMenu,'Label','Peak-Triggered Average','Callback',@PeakTrigAvgFcn);
uimenu(handles.MainFigDisplayMenu,'Label','Dip-Triggered Average','Callback',@DipTrigAvgFcn);
uimenu(handles.MainFigDisplayMenu,'Label','Template Matching','Callback',@TemplateMatchingFcn);
uimenu(handles.MainFigDisplayMenu,'Label','Save EPSCs to file','Callback',@TemplateMatchingSaveFcn);
uimenu(handles.MainFigDisplayMenu,'Label','Load EPSCs from file','Callback',@TemplateMatchingLoadFcn);

% Horizontal scaling controls
handles.HorizScaleLookup=[.5 1 2 5 10];
% these are horizontal scales to be chosen from: referring to fractions of a second
% (but will be used at different powers of 10)
handles.HorizScaleInButton=uicontrol('Style','PushButton',...
    'parent',handles.displayfig,...
    'Units','Normalized',...
    'Position',[.12 .065 .045 .045],...
    'String','In',...
    'Callback',@HorizScaleInButtonFcn);
handles.HorizScaleOutButton=uicontrol('Style','PushButton',...
    'parent',handles.displayfig,...
    'Units','Normalized',...
    'Position',[.166 .065 .045 .045],...
    'String','Out',...
    'Callback',@HorizScaleOutButtonFcn);
handles.HorizMinEditBox=uicontrol('Style','Edit',...
    'parent',handles.displayfig,...
    'Units','Normalized',...
    'Position',[.12 .0225 .045 .039],...
    'Callback',@HorizMinSetFcn);
handles.HorizMaxEditBox=uicontrol('Style','Edit',...
    'parent',handles.displayfig,...
    'Units','Normalized',...
    'Position',[.166 .0225 .045 .039],...
    'Callback',@HorizMaxSetFcn);
handles.TimeUnitsLabel = uicontrol('style','text',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.12 .005 .1 .015],...
    'string','Seconds',...
    'backgroundcolor',[.8 .8 .8]);
handles.MoveLeftButton=uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.24 .01 .03 .1],...
    'string','<',...
    'callback',@MoveLeftFcn);
handles.MoveRightButton=uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.271 .01 .03 .1],...
    'string','>',...
    'callback',@MoveRightFcn);
handles.ResetAxesButton=uicontrol('style','pushbutton',...
    'units','normalized',...
    'position',[.33 .01 .05 .1],...
    'string',{'Reset';'Axes'},...
    'Callback',@ResetAxesButtonFcn);

% Create chunk GUI objects (buttons, textboxes, etc)
handles.PrevChunk = uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.5-.06 .065 .06 .045],...
    'string','Prev Chunk',...
    'enable','off',...
    'callback',@PrevChunkFcn);
handles.NextChunk = uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.5 .065 .06 .045],...
    'string','Next Chunk',...
    'enable','off',...
    'callback',@NextChunkFcn);
handles.CurrChunkEditBox=uicontrol('Style','Edit',...
    'parent',handles.displayfig,...
    'Units','Normalized',...
    'Position',[0.5+0.06 0.065 0.07 0.022],...
    'enable','off',...
    'Callback',@CurrChunkEditFcn);
handles.CurrChunkLabel = uicontrol('style','text',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[0.5+0.06 .065+.022 0.07 0.022],...
    'string','Current Chunk',...
    'enable','off',...
    'backgroundcolor',[.8 .8 .8]);
handles.PrevIter = uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.5-.06 .01 .06 .045],...
    'string','Prev Iter.',...
    'enable','off',...
    'callback',@PrevIterFcn);
handles.NextIter = uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.5 .01 .06 .045],...
    'string','Next Iter.',...
    'enable','off',...
    'callback',@NextIterFcn);
handles.CurrIterEditBox=uicontrol('Style','Edit',...
    'parent',handles.displayfig,...
    'Units','Normalized',...
    'Position',[0.5+0.06 0.01 0.07 0.022],...
    'enable','off',...
    'Callback',@CurrIterEditFcn);
handles.CurrIterLabel = uicontrol('style','text',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[0.5+0.06 .01+.022 0.07 0.022],...
    'string','Current Iter.',...
    'enable','off',...
    'backgroundcolor',[.8 .8 .8]);

% Create measure buttons
handles.SingleMeasureButton = uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.65 .06 .12 .05],...
    'string','Single Measure',...
    'callback',@SingleMeasureFcn);%see comments below
handles.PairedMeasureButton = uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.65 .01 .12 .05],...
    'string','Paired Measure',...
    'callback',@PairedMeasureFcn);%see comments below
handles.SingleMeasureModeButton = uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.78 .06 .05 .05],...
    'string','Mode:',...
    'callback',@SingleMeasureModeFcn);
handles.PairedMeasureModeButton = uicontrol('style','pushbutton',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.78 .01 .05 .05],...
    'string','Mode:',...
    'callback',@PairedMeasureModeFcn);
handles.SingleMeasureMode = 'Click';%default
handles.PairedMeasureMode = 'Click & Click';%default
handles.SingleMeasureModeIndic = uicontrol('style','text',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.84 .0725 .14 .025],...
    'backgroundcolor',[.9 .9 .9],...
    'string',handles.SingleMeasureMode);
handles.PairedMeasureModeIndic = uicontrol('style','text',...
    'parent',handles.displayfig,...
    'units','normalized',...
    'position',[.84 .0225 .14 .025],...
    'backgroundcolor',[.9 .9 .9],...
    'string',handles.PairedMeasureMode);

% Single Measures Initialization
handles.SingleMeasuresMatrix.MeasureNumber = [];
handles.SingleMeasuresMatrix.FileName = {};
handles.SingleMeasuresMatrix.ChannelName = {};
handles.SingleMeasuresMatrix.xs = [];
handles.SingleMeasuresMatrix.ys = [];
handles.SingleMeasuresMatrix.MeasureMode = {};
handles.SingleMeasuresMatrix.Comments = {};

% Paired Measures Initialization
handles.PairedMeasuresMatrix.MeasureNumber = [];
handles.PairedMeasuresMatrix.FileName = {};
handles.PairedMeasuresMatrix.ChannelName = {};
handles.PairedMeasuresMatrix.x1s = [];
handles.PairedMeasuresMatrix.y1s = [];
handles.PairedMeasuresMatrix.x2s = [];
handles.PairedMeasuresMatrix.y2s = [];
handles.PairedMeasuresMatrix.xdiffs = [];
handles.PairedMeasuresMatrix.ydiffs = [];
handles.PairedMeasuresMatrix.MeasureMode = {};
handles.PairedMeasuresMatrix.Comments = {};

% Load the data
% objuserdata.channelaxes=1:length(handles.includedchannels);
% handles.obj.userdata=objuserdata;
handles.obj=info.ObjInfo;
handles.obj.HwInfo=info.HwInfo;
LoadDecimateStoreData(handles.displayfig,handles.path,handles.includedchannels,handles.obj.SampleRate);%load in data when the figure creation is finished
handles.sizedata = getappdata(handles.displayfig,'sizedata');

% Chunkerize or not
% THINGS TO DO IF CHUNK MODE ON
% PlotMinWithinXlims will plot two lines chunk + overlay
% Double up axes for two plots on one axes
% Scoot axes right side in so room for ticks
% what's gonna happen to measure button?!?!
% Handle pulse,target,iteration advancement correctly
if chunkquest==1 || isfield(handles,'rechunk')
    button=questdlg('Do you want to chunkerize your data?','Chunkin'' time');
    if strcmp(button,'Yes')
        handles.ChunkMode=1;
        handles=ChunkerizeData(handles);
        if laserquest~=0
            handles.ShowLaserPower=laserquest;
        end
        if stackquest==1
            button=questdlg('Do you want to split out the pulses to different axes?','Stack pulses?','No');
            if strcmp(button,'Yes')
                handles.StackPulses=1;
                totalaxes=handles.PulsesPerTarget;
                handles.ChannelAxesMatrix=zeros(totalaxes,handles.includedchannels);
                handles.ChannelAxesMatrix(:,handles.ChunkerizerParams.RecordChannel)=1;
                handles.ChannelAxesMatrix(:,handles.ChunkerizerParams.StimChannel)=1;
                handles.PostPulseStartSecondsToView=1;
                handles.xlims=[0 handles.PostPulseStartSecondsToView];
            end
        else
            handles.StackPulses=0;
            handles.ChannelAxesMatrix=eye(length(handles.includedchannels));
            totalaxes=length(handles.includedchannels);
            handles.xlims=ConvertChunkToXlims(handles);
        end
        handles.displaylength=handles.xlims(2)-handles.xlims(1);
        set(handles.HorizMinEditBox,'String',num2str(handles.xlims(1),3));
        set(handles.HorizMaxEditBox,'String',num2str(handles.xlims(2),3));
        set(handles.PrevChunk,'enable','on');
        set(handles.NextChunk,'enable','on');
        set(handles.CurrChunkEditBox,'enable','on');
        set(handles.CurrChunkLabel,'enable','on');
        set(handles.CurrChunkEditBox,'string','1');
        set(handles.PrevIter,'enable','on');
        set(handles.NextIter,'enable','on');
        set(handles.CurrIterEditBox,'enable','on');
        set(handles.CurrIterLabel,'enable','on');
        set(handles.CurrIterEditBox,'string','1');
    else
        handles.ChunkMode=0;
        handles.StackPulses=0;
        handles.ChannelAxesMatrix=eye(length(handles.includedchannels));
        handles.displaylength=handles.sizedata(1)/handles.obj.SampleRate;%total seconds
        handles.xlims=[1/handles.obj.SampleRate handles.displaylength];
        set(handles.HorizMinEditBox,'string','0');
        set(handles.HorizMinEditBox,'string',num2str(handles.displaylength));
        % Set up the axes
        totalaxes=length(handles.includedchannels);%taking data from info
    end
elseif chunkquest==0
    handles.ChunkMode=0;
    handles.StackPulses=0;
    handles.ChannelAxesMatrix=eye(length(handles.includedchannels));
    handles.displaylength=handles.sizedata(1)/handles.obj.SampleRate;%total seconds
    handles.xlims=[1/handles.obj.SampleRate handles.displaylength];
    set(handles.HorizMinEditBox,'string','0');
    set(handles.HorizMinEditBox,'string',num2str(handles.displaylength));
    % Set up the axes
    totalaxes=length(handles.includedchannels);%taking data from info
end

% GUI layout constants and initialization
handles.vertzoomfactor=2;
vertgap=.015;
bottombottom=.145;
toptop=.98;
totalvertspace=toptop-bottombottom;
eachvertspace=totalvertspace/totalaxes;
axesheight=eachvertspace-vertgap;
for a = 1:totalaxes;%create each axes, equally spaced between certain set points
    ainv = totalaxes - a + 1;%last, last-1, last-2... 2, 1... just to make POSITION of axes right
    handles.axeshandles(a)=axes('parent',handles.displayfig,...
        'Units','Normalized',...
        'position',[.12 bottombottom+(eachvertspace*(ainv-1))+vertgap .78 axesheight],...
        'DrawMode','fast',...
        'XLimMode','manual',...%if not, then auto rescales sometimes on one axis
        'YLimMode','manual',...
        'XTick',[],...
        'ButtonDownFcn',@ZoomFcn);%immediately activate zooming
    if a == totalaxes(end)
        set(handles.axeshandles(a),'XTickMode','auto')
    end
    
    axesposition=get(handles.axeshandles(a),'position');%do relative to these
    middlevert=axesposition(2)+(.5*axesposition(4));
    widthslides=.047;
    heightslides=axesposition(4)/3;
    
    try
        handles.VertSlideUp(a)=uicontrol('Style','PushButton',...
            'parent',handles.displayfig,...
            'Units','Normalized',...
            'Position',[axesposition(1)+axesposition(3)+.005 axesposition(2)+2*heightslides widthslides heightslides],...
            'String','^',...
            'Callback',@VertSlideUpFcn);%slides up by 1/4 field of view
        handles.VertSlideDown(a)=uicontrol('Style','PushButton',...
            'parent',handles.displayfig,...
            'Units','Normalized',...
            'Position',[axesposition(1)+axesposition(3)+.005 axesposition(2)+heightslides widthslides heightslides],...
            'String','v',...
            'Callback',@VertSlideDownFcn);%slides down by 1/4 field of view
        
        handles.VertScaleIn(a)=uicontrol('Style','PushButton',...
            'parent',handles.displayfig,...
            'Units','Normalized',...
            'Position',[axesposition(1)+axesposition(3)+.005+widthslides axesposition(2)+2*heightslides widthslides heightslides],...
            'String','In',...
            'Callback',@VertScaleInFcn);
        handles.VertScaleOut(a)=uicontrol('Style','PushButton',...
            'parent',handles.displayfig,...
            'Units','Normalized',...
            'Position',[axesposition(1)+axesposition(3)+.005+widthslides axesposition(2)+heightslides widthslides heightslides],...
            'String','Out',...
            'Callback',@VertScaleOutFcn);
    catch
        error ('too many channels selected')
    end
    
    handles.VertMinEditBox(a)=uicontrol('Style','Edit',...
        'parent',handles.displayfig,...
        'Units','Normalized',...
        'Position',[axesposition(1)+axesposition(3)+.005 axesposition(2) widthslides heightslides],...
        'Callback',@VertMinSetFcn);
    handles.VertMaxEditBox(a)=uicontrol('Style','Edit',...
        'parent',handles.displayfig,...
        'Units','Normalized',...
        'Position',[axesposition(1)+axesposition(3)+.005+widthslides axesposition(2) widthslides heightslides],...
        'Callback',@VertMaxSetFcn);
    handles.ChannelsLabelAxes(a)=axes('Units','Normalized',...%an axes for putting names of channels.  Has xylims of 0,1
        'parent',handles.displayfig,...
        'Position',[0.01 axesposition(2) .08 axesposition(4)],...%to left of the plotting axes, same vert position and height
        'Visible','Off');
    handles.ChannelMeasureCheckbox(a) = uicontrol('style','checkbox',...
        'parent',handles.displayfig,...
        'units','normalized',...
        'position',[0.01 axesposition(2)+(axesposition(4)/4) .08 .02],...
        'backgroundcolor',[.8 .8 .8],...
        'string','Measure',...
        'FontSize',10,...
        'value',0);
end

% Create laser power GUI edit box
if laserquest~=0;
    if isfield(handles,'ChunkerizerParams')
        LaserPowerTextPos=get(handles.axeshandles(handles.ChunkerizerParams.StimChannel),...
            'Position');
        LaserPowerTextPos(1)=LaserPowerTextPos(1)+LaserPowerTextPos(3)/2;
        LaserPowerTextPos(2)=LaserPowerTextPos(2)+LaserPowerTextPos(4)*.8;
        LaserPowerTextPos(3)=0.1;
        LaserPowerTextPos(4)=0.04;
        handles.LaserEditBox = uicontrol('style','Edit',...
            'parent',handles.displayfig,...
            'units','normalized',...
            'position',LaserPowerTextPos,...
            'backgroundcolor','w',...
            'string','init');
    end
end

set(handles.axeshandles,'xlim',handles.xlims);%set all xlims the same
PlotMinWithinXlims(handles,get(handles.axeshandles(1),'XLim'));

decdata=getappdata(handles.displayfig,'decdata');
for a=1:totalaxes
    channelsthisaxes=find(handles.ChannelAxesMatrix(a,:));
    % 	channelsthisaxes=find(objuserdata.channelaxes==a);%find which channels are in that axes
    if ~isempty(channelsthisaxes)
        maxy=0;
        miny=0;
        for b=1:length(channelsthisaxes);%for each channel this axes
            maxy=max([maxy,max(decdata.data{1}(:,channelsthisaxes(b)))]);
            miny=min([miny,min(decdata.data{1}(:,channelsthisaxes(b)))]);
            thischannelindex = channelsthisaxes(b);
            if handles.StackPulses==1
                ThisUnits=handles.obj.Channel(handles.ChunkerizerParams.RecordChannel).Units;
                ThisPulseString=['Pulse' num2str(a)];
                handles.ChannelLabels(a,b)=text(0.01, .75,...
                    [ThisPulseString,sprintf('\n'),'(',ThisUnits,')'],...
                    'Parent',handles.ChannelsLabelAxes(a),...
                    'FontSize',10);
                %             elseif handles.ChunkMode==1
                %                 handles.ChannelLabels(a,b)=text(0.01, .75,...
                %                     [handles.obj.Channel(thischannelindex).ChannelName,sprintf('\n'),'(',handles.obj.Channel(thischannelindex).Units,')'],...
                %                     'Parent',handles.ChannelsLabelAxes(a),...
                %                     'FontSize',10);
            else % rock it old school, i.e. the way brendon did it, i.e. correctly... (AP 20081205)
                objuserdata.channelaxes=1:length(handles.includedchannels);
                channelsthisaxes=find(objuserdata.channelaxes==a);%find which channels are in that axes
                thischannelindex = handles.includedchannels(channelsthisaxes(b));
                handles.ChannelLabels(a,b)=text(0.01, .75,...
                    [handles.obj.Channel(thischannelindex).ChannelName,sprintf('\n'),'(',handles.obj.Channel(thischannelindex).Units,')'],...
                    'Parent',handles.ChannelsLabelAxes(a),...
                    'FontSize',10);
            end
        end
        if maxy==0  % Added by Adam 01/18/07 to account for a channel that is ALWAYS zero!
            if miny==0
                maxy=1;
            end
        end
        handles.maxy(a)=maxy+abs(.05*maxy);%go up 5%
        handles.miny(a)=miny-abs(.05*miny);%go up 5%
        set(handles.axeshandles(a),'ylim',[handles.miny(a) handles.maxy(a)]);
        set(handles.VertMinEditBox(a),'String',num2str(handles.miny(a),3));
        set(handles.VertMaxEditBox(a),'String',num2str(handles.maxy(a),3));
    end
end


% for a=1:totalaxes;%for every axes
%     channelsthisaxes=find(handles.ChannelAxesMatrix(a,:));
% % 	channelsthisaxes=find(objuserdata.channelaxes==a);%find which channels are in that axes
%     if ~isempty(channelsthisaxes)
%         for b=1:length(channelsthisaxes);%for each channel this axes
%             thischannelindex = handles.includedchannels(channelsthisaxes(b));
%             handles.ChannelLabels(a,b)=text(0.01, .75,...
%                 [handles.obj.Channel(thischannelindex).ChannelName,sprintf('\n'),'(',handles.obj.Channel(thischannelindex).Units,')'],...
%                 'Parent',handles.ChannelsLabelAxes(a),...
%                 'FontSize',10);%,'color',thiscolor);
%         end
%     end
% end
guidata(handles.displayfig,handles);

set(handles.displayfig,'Visible','on');
varargout{1} = 1;

function LoadDecimateStoreData(fighandle,filename,includedchannels,samplingrate)
%uses handles to store decimated data into the appdata of the main figure


[pathstr, name, ext] = fileparts(filename);
%data must be output as matrix having dims = pointnumber, channelnumber
if strcmp(ext,'.daq')
    alldata=bdaqread(filename,'Channels',includedchannels);
    %read data back to this workspace from the written file...only read specified channels
elseif strcmp(ext,'.bin')
    alldata=paq2lab(filename,'Channels',includedchannels);
elseif strcmp(ext,'.paq')
    alldata=paq2lab(filename,'Channels',includedchannels);
elseif strcmp(ext,'.fig');%again... assuming the fig is like a triggered average fig.
    ReadFig = hgload(filename);
    set(ReadFig,'visible','off');
    try
        ReadFigDat = getappdata(ReadFig,'data');
    catch
        errdlg ('Invalid Figure.  Must have "data" in appdata containing a structure of line data')
        return
    end
    close(ReadFig)
    clear ReadFig
    averagedchannels = [2 4 6];
    includedchannels = averagedchannels(includedchannels);
    for cidx = 1:length(includedchannels);
        alldata(:,cidx) = ReadFigDat.data{includedchannels(cidx)}';
    end
end

ordmag=floor(log10(length(alldata)));%get number of orders of magnitude occupied by the length
try
    [rangequest,qmquestion,chunkquest,channels,filter,stackquest,laserquest]=ephysdefaults;
catch
    filter=0;
end
if filter~=0;
    alldata=medfilt1(double(alldata),filter);
end
decdata.data{1}=alldata;%store data for plotting
clear alldata;%save space
% decdata.lengths(1)=length(alldata);
decdata.decfactor(1)=1;
decdata.timepoints{1}=[];%this is just to take up space... don't want another vector as long as alldata
%instead, for the special case of data{1}, time points will be the raw
%index numbers of the points.

waithandle = waitbar(0,'Decimating and Storing Data');%for user
for a=1:ordmag-1;%for each order of magnitude
    [decdata.data{a+1},decdata.timepoints{a+1}]=Decimatebymaxmin(decdata.data{1},10^a);
    decdata.timepoints{a+1}=decdata.timepoints{a+1}/samplingrate;
    %store a version of the data decimated by that order of magnitude
    %     decdata.lengths(a+1)=length(decdata.data{a+1});
    decdata.decfactor(a+1)=10^a;%store the decimation factor
    waithandle = waitbar(a/(ordmag-1),waithandle);
end
setappdata(fighandle,'decdata',decdata);
setappdata(fighandle,'sizedata',size(decdata.data{1}));
close(waithandle)

function varargout=Decimatebymaxmin(vect,factor)

flipindicator=0;
if length(size(vect))>2;
    error('input must be 2D or less')
elseif size(vect,1)<size(vect,2);%want size to be big x small (data stream by stream number)
    vect=vect';
    flipindicator=1;
end
numchans = size(vect,2);

factor=factor*2;
remain=rem(length(vect),factor);
if remain>0;%if number of points is not divisible by factor
    newvect=vect(1:end-remain,:);
    if remain>2;%if enough points to get a max and a min
        temp=vect(end-remain+1:end,:);
        temp(end+1:factor,:)=repmat(mean(temp,1),[factor-size(temp,1) 1]).*ones(factor-size(temp,1),size(temp,2));
        newvect(end+1:end+factor,:)=temp;
    end
else%if divisible by factor
    newvect=vect;
end
clear vect

newvect=reshape(newvect,[factor length(newvect)/factor  size(newvect,2)]);%set up for factor-fold decimation of data

[nvmax,maxpt]=max(newvect,[],1);
[nvmin,minpt]=min(newvect,[],1);

newvect=zeros(2*size(nvmax,2),numchans);
newvect(1:2:size(newvect,1)-1,:)=squeeze(nvmax);%odd numbered points are maxs
newvect(2:2:size(newvect,1),:)=squeeze(nvmin);%even points are local mins

if nargout>1;
    temp=((0:2:(size(newvect,1)-1))*(factor/2));
    timepoints=sort(cat(2,temp+round(factor/2/3),temp+round(factor/3)))';
end

if flipindicator==1;
    newvect=newvect';
    timepoints=timepoints';
end

varargout{1}=newvect;
if nargout>1;
    varargout{2}=timepoints;
end

function PlotMinWithinXlims(handles,xlims)
%this function takes the xlims it's given and takes handles, which has
%decimated data in its appdata (decdata structure) and plots the
%minimum-necessary-resolution data between those xlims.  It also deletes
%previous data and it plots specific lines in specific axes, as indicated
%by channelaxes in the daqobject userdata (handles.obj.userdata).

decdata=getappdata(handles.displayfig,'decdata');
if handles.StackPulses==1;
    for i=1:handles.PulsesPerTarget
        ThisPulseStartTime(i)=handles.PulseStarts(handles.CurrentTargetNumber,handles.CurrentIteration,i)/handles.obj.SampleRate;
        TwoDxlims(i,:)=[xlims(1)+ThisPulseStartTime(i) xlims(2)+ThisPulseStartTime(i)];
    end
end

set(handles.axeshandles(1),'units','pixels');%set units of an axes to pixels so can measure in pixels
pixelwidth=get(handles.axeshandles(1),'position');%record the pixels occupied by that axes
set(handles.axeshandles(1),'units','normalized');%restore the "units" property to its original value
pixelwidth=pixelwidth(3);%save just the width (in pixels) of the axes (not all about its location)

% Find shortest version of decimated data that has at least 2x as many points as the axes has pixels.
% If the needed resolution is greater than exists, just use the highest available.
rationow=(((xlims(2)-xlims(1))*handles.obj.SampleRate)/pixelwidth);
index=max(find(decdata.decfactor<(rationow/2)));
if isempty(index);
    index=1;
end

% objuserdata=handles.obj.userdata;

for a=1:length(handles.axeshandles);%for each axes
    channelsthisaxes=find(handles.ChannelAxesMatrix(a,:));
    %     channelsthisaxes=find(objuserdata.channelaxes==a);%find which channels are in that axes
    if ~isempty(channelsthisaxes);
        delete(get(handles.axeshandles(a),'children'));%delete all lines in the window
        DataToPlot=decdata.data{index}(:,channelsthisaxes);%get the data for just the channels in this axes
        thisaxes=handles.axeshandles(a);
        if index==1;
            if handles.StackPulses==1
                yindices=round(TwoDxlims(a,1)*handles.obj.SampleRate):round(TwoDxlims(a,2)*handles.obj.SampleRate);
                yindices(yindices>length(DataToPlot))=[];
                yindices(yindices<1)=[];
                %                 xpoints=xlims(1):1/handles.obj.SampleRate:xlims(2);
                %                 yindices=yindices(1:size(xpoints));
                
                xpoints=(yindices/handles.obj.SampleRate)-ThisPulseStartTime(a);
                ydata=DataToPlot(yindices,:);
            else
                yindices=round(xlims(1)*handles.obj.SampleRate):round(xlims(2)*handles.obj.SampleRate);
                yindices(yindices>length(DataToPlot))=[];
                yindices(yindices<1)=[];
                xpoints=yindices/handles.obj.SampleRate;
                ydata=DataToPlot(yindices,:);
            end
            %             if handles.overlay
            %                 xpoints = [xpoints;xpoints];
            %                 ydata = [DataToPlot(yindices,:); stimvector];
            %             else
            %                 ydata = DataToPlot(yindices,:);
            %             end
            try
                if handles.ShowLaserPower~=0;
                    if a==handles.ChunkerizerParams.StimChannel
                        set(handles.LaserEditBox,'String',num2str(10*round(max(ydata)/(handles.ShowLaserPower/10))));
                    end
                end
            catch
            end
            maxy=max(ydata);
            miny=min(ydata);
            handles.maxy(a)=maxy+abs(.05*maxy);%go up 5%
            handles.miny(a)=miny-abs(.05*miny);%go up 5%
            set(handles.axeshandles(a),'ylim',[handles.miny(a) handles.maxy(a)]);
            line(xpoints,ydata,'parent',thisaxes,'HitTest','Off')%hit test off for easier zooming
            set(handles.VertMinEditBox(a),'String',num2str(handles.miny(a),3));
            set(handles.VertMaxEditBox(a),'String',num2str(handles.maxy(a),3));
        else
            if handles.StackPulses==1
                bottime=max(find(decdata.timepoints{index}<TwoDxlims(a,1)));%find last point less than the bottom xlim
                if isempty(bottime); bottime=1; end
                toptime=min(find(decdata.timepoints{index}>TwoDxlims(a,2)));%find first point greater than the top xlim
                if isempty(toptime); toptime=length(decdata.timepoints{index}); end
            else
                bottime=max(find(decdata.timepoints{index}<xlims(1)));%find last point less than the bottom xlim
                if isempty(bottime); bottime=1; end
                toptime=min(find(decdata.timepoints{index}>xlims(2)));%find first point greater than the top xlim
                if isempty(toptime); toptime=length(decdata.timepoints{index}); end
            end
            %             if handles.overlay
            %                 xpoints = [xpoints;xpoints];
            %                 ydata = [DataToPlot(yindices,:); stimvector];
            %             else
            %                 ydata = DataToPlot(yindices,:);
            %             end
            line(decdata.timepoints{index}(bottime:toptime),DataToPlot(bottime:toptime,:),'parent',thisaxes,'HitTest','Off')%hit test off for easier zooming
        end
    end
end

function LocalMenuLoadDataFcn(obj,event)

handles=guidata(obj);%get general info from the figure

[pathstr, name, ext]=fileparts(handles.path);
[filename,pathname,FilterIndex]=uigetfile(handles.filetypes,'Choose a data file',[pathstr,'\']);
if ~FilterIndex
    return
end
%reinitialize handles for next version of Viewer
handles2.filetypes = {'*.daq;*.bin;*.paq;*.fig'};%because ViewerMainFcn expects a handles with very little (ie from top fcn)
handles2.path=[pathname,filename];
handles2.position = get(handles.displayfig,'position');
try
    vieweroutput = ViewerMainFcn(handles2);%vieweroutput it set to 1 at end of function
    % ie if all went well... if not all goes well the expectation of
    % an output creates an error... so go to the catch
catch
    vieweroutput = 0;
end
if vieweroutput
    close(handles.displayfig);
end

function LocalMenuLoadNextFileFcn(obj,event)

handles=guidata(obj);
[pathstr, name, ext]=fileparts(handles.path);
structlist=dir(fullfile(pathstr,['*' ext]));
celllist=struct2cell(structlist);
goodfilelist=sort_nat(celllist(1,:));
fileidx=strmatch([name ext],goodfilelist);
try
    filenametoget=goodfilelist(fileidx+1);
    handles2.filetypes = {'*.daq;*.bin;*.paq;*.fig'};%because ViewerMainFcn expects a handles with very little (ie from top fcn)
    try
        if ismac
            handles2.path=fullfile(pathstr,filenametoget{1});
        else
            handles2.path=fullfile(pathstr,'\',filenametoget{1});
        end
    catch
        handles2.path=fullfile(pathstr,'\',filenametoget{1});
    end
    handles2.position = get(handles.displayfig,'position');
    try
        vieweroutput = ViewerMainFcn(handles2);%vieweroutput it set to 1 at end of function
        % ie if all went well... if not all goes well the expectation of
        % an output creates an error... so go to the catch
    catch
        vieweroutput = 0;
    end
    if vieweroutput
        close(handles.displayfig);
    end
catch
    msgbox('There is no next file in this folder.')
end


function LocalMenuLoadPrevFileFcn(obj,event)

handles=guidata(obj);
[pathstr, name, ext]=fileparts(handles.path);
structlist=dir(fullfile(pathstr,['*' ext]));
celllist=struct2cell(structlist);
goodfilelist=sort_nat(celllist(1,:));
fileidx=strmatch([name ext],goodfilelist);
try
    filenametoget=goodfilelist(fileidx-1);
catch
    msgbox('There is no previous file in this folder.')
end
handles2.filetypes = {'*.daq;*.bin;*.paq;*.fig'};%because ViewerMainFcn expects a handles with very little (ie from top fcn)
try
    if ismac
        handles2.path=fullfile(pathstr,filenametoget{1});
    else
        handles2.path=fullfile(pathstr,'\',filenametoget{1});
    end
catch
    handles2.path=fullfile(pathstr,'\',filenametoget{1});
end
handles2.position = get(handles.displayfig,'position');
try
    vieweroutput = ViewerMainFcn(handles2);%vieweroutput it set to 1 at end of function
    % ie if all went well... if not all goes well the expectation of
    % an output creates an error... so go to the catch
catch
    vieweroutput = 0;
end
if vieweroutput
    close(handles.displayfig);
end



function LocalMenuReChunkerizeDataFcn(obj,event)
handles=guidata(obj);
if handles.chunkquest==1
    handles=ChunkerizeData(handles);
    guidata(handles.displayfig,handles);%pass data back to the figure
else
    EphysViewer(handles.path,1);
end

function LocalMenuLoadFluorescentTracesFcn(obj,event)
% Load traces, signals, and movie timing
handles=guidata(obj);

% Ask user when the movie started
handles.MovieStartTime=APinputdlg('Movie start time (sec)','Movie start');
if isempty(handles.MovieStartTime)
    return
end
handles.MovieStartTime=str2double(handles.MovieStartTime(1));

% Load traces
[filename,pathname,FilterIndex]=uigetfile('*.mat','Choose a MAT file with traces');
if FilterIndex
    tracefile=fullfile(pathname,filename);
else
    return
end
load(tracefile);
handles.tracefile=tracefile;
handles.FluorTraces=traces;

% Load signal onsets
[filename,pathname,FilterIndex]=uigetfile('*.mat','Choose a MAT file with signal onsets');
if FilterIndex
    signalfile=fullfile(pathname,filename);
else
    return
end
load(signalfile);
handles.signalfile=signalfile;
handles.FluorSignals=onsets;

% Load movie timing file (txt, tab delimited, with one header row)
[filename,pathname,FilterIndex]=uigetfile('*.txt','Choose a TXT file with movie timing');
if FilterIndex
    timingfile=fullfile(pathname,filename);
else
    return
end
handles.timingfile=timingfile;
movietiming=dlmread(timingfile,'\t',1,0);
handles.FluorTiming=movietiming(:,3:4);
handles.FrameStartTimes=handles.FluorTiming(:,1);

% Calculate the maximum number of signals from cell handles.FluorSignals
MaxNumSignals=0;
for i=1:size(handles.FluorSignals,2);
    if size(handles.FluorSignals{i},1) > MaxNumSignals
        MaxNumSignals=size(handles.FluorSignals{i},1);
    end
end

% Load handles.FluorSignalTimes (referenced to electrophysiology time)
% Reference handles.FluorTraces to electrophysiology time
handles.FluorSignalTimes=NaN(size(handles.FluorSignals,2),MaxNumSignals);
for i=1:size(handles.FluorSignals,2);
    if ~isempty(handles.FluorSignals{i})
        for k=1:length(handles.FluorSignals{i}')
            ThisSigStartMovieTime=handles.FrameStartTimes(handles.FluorSignals{i}(k));
            handles.FluorSignalTimes(i,k)=handles.MovieStartTime + ThisSigStartMovieTime;
        end
    end
end

% Plot and return data to figure
handles=PlotFiredCell(handles);
guidata(handles.displayfig,handles);%pass data back to the figure

function handles=PlotFiredCell(handles)
% Find the cell that was fired by the pulse and plot its calcium transient.
if isfield(handles,'FluorSignalTimes')
    ThisPulseStart=handles.PulseStarts(handles.CurrentTargetNumber,handles.CurrentIteration)/handles.obj.SampleRate;
    ThisPulseStop=handles.PulseStops(handles.CurrentTargetNumber,handles.CurrentIteration)/handles.obj.SampleRate;
    [CellWithSig,SigNum]=find((handles.FluorSignalTimes > ThisPulseStart-0.05 & (handles.FluorSignalTimes < ThisPulseStop+0.05)));
    NumCellsFired=length(CellWithSig);
    % if there is not a FiredCellFig already, set it up
    FiredCellFigName=['Fired ' handles.filename];
    if ~isfield(handles,'FiredCellFig')
        handles.FiredCellFig=figure('name',FiredCellFigName,'numbertitle','off');
    else
        try
            PossiblyOverwrittenFiredCellFig=get(handles.FiredCellFig);
            if ~strcmp(PossiblyOverwrittenFiredCellFig.Name,FiredCellFigName);
                handles.FiredCellFig=figure('name',FiredCellFigName,'numbertitle','off');
            else
                figure(handles.FiredCellFig)
                set(gcf,'name',FiredCellFigName,'numbertitle','off');
            end
        catch
            handles.FiredCellFig=figure('name',FiredCellFigName,'numbertitle','off');
        end
    end
    xlims=handles.xlims;
    xlims(1)=xlims(1);
    xlims(2)=xlims(2);
    if NumCellsFired==0
        clf  % if no cells are fired, clear the fired cell fig
    end
    for n=1:NumCellsFired
        subplot(NumCellsFired,1,n);
        bar(handles.FrameStartTimes + handles.MovieStartTime,handles.FluorTraces(CellWithSig(n),:));
        hold on;
        plot(handles.FluorSignalTimes(CellWithSig(n),:),0,'o','Color','r');
        hold off;
        xlim(xlims)
        %         ylim([-200 1000])
    end
end

function LocalMenuLoadSynchronousTraceFcn(obj,event)
% Load traces and movie timing
handles=guidata(obj);

% Load traces
[filename,pathname,FilterIndex]=uigetfile('*.mat','Choose a MAT file with trace');
if FilterIndex
    tracefile=fullfile(pathname,filename);
else
    return
end
handles.tracefile=tracefile;
temp=load(tracefile);
tempnames=fieldnames(temp);
tempfirstname=strcat('temp.',tempnames{1});
handles.SynchTrace=eval(tempfirstname);

% WORKED WITH TXT FILES WITH TIMING DATA FROM HCIMAGE CIRCA 2008
% Load trace timing file (txt, tab delimited, with one header row)
% [filename,pathname,FilterIndex]=uigetfile('*.txt','Choose a TXT file with trace timing');
% if FilterIndex
%     timingfile=fullfile(pathname,filename);
% else
%     return
% end
% handles.timingfile=timingfile;
% tracetiming=dlmread(timingfile,'\t',1,0);
% handles.TraceTiming=tracetiming(:,3:4);
% handles.TraceStartTimes=handles.TraceTiming(:,1);
% handles.MeanExposureTime=mean(handles.TraceTiming(:,2));
% longframes = find(handles.TraceTiming(:,2) > handles.MeanExposureTime + (handles.MeanExposureTime / 2));
% if longframes
%     answer=inputdlg(sprintf('There are long frames at %d\nAdd how many frames?',longframes),'Long frames');
%     FrameAdjust=str2double(answer);
%     handles.TraceStartTimes(longframes:end) = handles.TraceStartTimes(longframes:end) + FrameAdjust*handles.MeanExposureTime;
% end

% Load trace timing file 
[filename,pathname,FilterIndex]=uigetfile('*.mat','Choose a MAT file with trace timing');
if FilterIndex
    timingfile=fullfile(pathname,filename);
else
    return
end
handles.timingfile=timingfile;
temp=load(timingfile);
tempnames=fieldnames(temp);
tempfirstname=strcat('temp.',tempnames{1});
handles.TraceTiming=eval(tempfirstname);

% Plot and return data to figure
handles=PlotSynchronousTrace(handles);
guidata(handles.displayfig,handles);%pass data back to the figure

function handles=PlotSynchronousTrace(handles)
% if there is not a SynchTraceFig already, set it up
SynchTraceFigName=['Trace ' handles.filename];
if ~isfield(handles,'SynchTraceFig')
    handles.SynchTraceFig=figure('name',SynchTraceFigName,'numbertitle','off');
else
    try
        PossiblyOverwrittenSynchTraceFig=get(handles.SynchTraceFig);
        if ~strcmp(PossiblyOverwrittenSynchTraceFig.Name,SynchTraceFigName);
            handles.SynchTraceFig=figure('name',SynchTraceFigName,'numbertitle','off');
        else
            figure(handles.SynchTraceFig)
            set(gcf,'name',SynchTraceFigName,'numbertitle','off');
        end
    catch
        handles.SynchTraceFig=figure('name',SynchTraceFigName,'numbertitle','off');
    end
end
% plot regular
plot(handles.TraceTiming,handles.SynchTrace);
xlim(get(handles.axeshandles(1),'xlim'));

% plot unconnected dots
% plot(handles.TraceStartTimes + handles.TraceStartTime,handles.SynchTrace,'LineStyle','none','Marker','.');

% do a bar/stem/plot plot (works well with baseline subtracted traces)
% bar(handles.TraceTiming,handles.SynchTrace);
% xlim(get(handles.axeshandles(1),'XLim'))
% ylim([1750 2100])


function handles=ChunkerizeData(handles)

decdata=getappdata(handles.displayfig,'decdata');

% Load VNT
button=questdlg('Do you want to load a VNT file associated with this data?','Load target coordinates');
if strcmp(button,'Yes')
    [filename,pathname,FilterIndex]=uigetfile('*.vnt;*.csv','Choose a VNT file');
    if FilterIndex
        handles.VNTfile=[pathname,filename];
        handles.VNTlist=load(handles.VNTfile);
        [path,name,ext]=fileparts(filename);
        if strcmp(ext,'.vnt')
            handles.VNTlist(end,:)=[]; % cut the end of the VNT file which should have -1 power
        end
        handles.VNTlist(handles.VNTlist(:,5)==1,:)=[]; % keep only "imaging" target
        [VNTlistXsorted,VNTlistXindices]=sortrows(handles.VNTlist,1); % sort the list by first column
        VNTlistXsortedNotDiff=~logical(diff(VNTlistXsorted)); % find places where X is not different
        VNTlistXsortedNotDiffXY=VNTlistXsortedNotDiff(:,1).*VNTlistXsortedNotDiff(:,2); % find places where X and Y are not different
        NumPulsesTimesIterGuess=find(VNTlistXsortedNotDiffXY==0,1,'first'); % find first duplicate row
        DiffIndices=diff(VNTlistXindices); % calculate how far apart rows were in original file
        FirstNonPulseIdx=find(DiffIndices~=1,1,'first'); % find the first instance where the rows were together
        if FirstNonPulseIdx==1
            NumPulsesGuess=1;
        elseif isempty(FirstNonPulseIdx)
            NumPulsesGuess=1;
        else
            NumPulsesGuess=FirstNonPulseIdx;
        end
        handles.NumIterFromVNT=NumPulsesTimesIterGuess/NumPulsesGuess;
        handles.NumTargetsFromVNT=size(handles.VNTlist,1)/(NumPulsesTimesIterGuess);
        handles.VNTXsortedIndices=VNTlistXindices;
        handles.VNTIdxsOfMatches=reshape(VNTlistXindices,handles.NumIterFromVNT*NumPulsesGuess,handles.NumTargetsFromVNT);
        handles.ChunkerizerParamsGuess={'1','2','0','0.1','10','10000', ...
            num2str(NumPulsesGuess), ...
            num2str(handles.NumTargetsFromVNT), ...
            num2str(handles.NumIterFromVNT),'10','2000','0'};
    end
elseif strcmp(button,'No')
    handles.ChunkerizerParamsGuess=[];
else
    return
end

% Chunkerize the data
[handles.ChunkerizedDataIndices,handles.PulseIndices,handles.ChunkerizerParams,handles.DataMissing]= ...
    Chunkerizer(decdata.data{1},handles.ChunkerizerParamsGuess);
if isempty(handles.ChunkerizedDataIndices)
    msgbox('Chunkerizing failed!')
    return
elseif handles.ChunkerizedDataIndices==Inf;
    return
end

% Initialize variables
handles.OverlayPlots=0;
handles.CurrentTargetNumber=1;
handles.CurrentIteration=1;
handles.NumTargets=size(handles.ChunkerizedDataIndices,1);
handles.NumIterations=size(handles.ChunkerizedDataIndices,3);
handles.PulsesPerTarget=handles.ChunkerizerParams.PulsesPerTarget;
handles.PulseStarts=zeros(handles.NumTargets,handles.NumIterations,handles.PulsesPerTarget);
handles.PulseStops=zeros(handles.NumTargets,handles.NumIterations,handles.PulsesPerTarget);
for i=1:handles.NumTargets
    for j=1:handles.NumIterations
        FirstPulse=handles.PulsesPerTarget*(handles.NumTargets*(j-1)+i-1)+1;
        LastPulse=handles.PulsesPerTarget*(handles.NumTargets*(j-1)+i);
        handles.PulseStarts(i,j,:)=handles.PulseIndices(FirstPulse:LastPulse,1);
        handles.PulseStops(i,j,:)=handles.PulseIndices(FirstPulse:LastPulse,2);
    end
end
handles.TotalNumberOfChunks=size(handles.ChunkerizedDataIndices,1)*handles.NumIterations;
handles.SortResults=zeros(handles.NumTargets,handles.NumIterations);
handles.SignalsFromTargets=zeros(handles.NumTargets,handles.NumIterations);
handles.DisplayScalingFactor=7;

% Load extensible indexing matrix
% Load correct target coordinates now that data has been properly chunkerized
if isfield(handles,'VNTIdxsOfMatches')
    if handles.DataMissing==1;
        % This is not going to work if there are iterations in the VNT! AP 20080718
        handles.Reindexer=1:handles.NumTargets;
    else
        handles.Reindexer=mod(handles.VNTIdxsOfMatches,handles.NumTargets);
    end
    if handles.NumIterFromVNT ~= handles.NumIterations
        handles.Reindexer=repmat(handles.Reindexer,handles.NumIterations,1);
    end
    handles.targetXidxs=round(handles.VNTlist(:,1));
    handles.targetYidxs=round(handles.VNTlist(:,2));
    
    % This code block handles the case when the number of pulses is
    % not equal to the VNT target list.
    % It should be unnecessary if parsing of the VNT target list is
    % properly handled using the reindexer.
    %     if size(OriginalChunks,1)~=length(handles.VNTlist)
    %         str1=sprintf('Number of pulses does not match number of targets from VNT.\n');
    %         str2=sprintf('Did you fire %d times?',handles.ChunkerizerParams.NumIterations);
    %         button=questdlg([str1 str2],'VNT/pulse mismatch');
    %         if strcmp(button,'Yes')
    %             handles.VNTrepeatTimes=handles.ChunkerizerParams.NumIterations;
    %         elseif strcmp(button,'No')
    %             handles.NumIterationsQuestion=APinputdlg('How many times were the targets fired?','Iterations');
    %             if isempty(handles.NumIterationsQuestion)
    %                 return
    %             end
    %             handles.VNTrepeatTimes=str2double(handles.NumIterationsQuestion);
    %         else
    %             return
    %         end
    %         handles.VNTlist=repmat(handles.VNTlist,handles.VNTrepeatTimes,1);
    %     end
    
    % OriginalChunks=ones(handles.NumTargets*handles.NumIterations,1);
    % OriginalChunks=find(OriginalChunks);
    % handles.PulsesPerTargetInVNT=size(handles.VNTIdxsOfMatches,1);
else
    FakeTargetList=1:handles.NumTargets*handles.NumIterations;
    FakeTargetList=reshape(FakeTargetList,handles.NumTargets,handles.NumIterations);
    handles.Reindexer=mod(FakeTargetList',handles.NumTargets);
    handles.Reindexer=handles.Reindexer;
end
handles.Reindexer(handles.Reindexer==0)=handles.NumTargets;
handles.Reindexer=sortrows(handles.Reindexer');

function newlims=ConvertChunkToXlims(handles)
CurrTarget=handles.Reindexer(handles.CurrentTargetNumber,handles.CurrentIteration);
CurrIter=handles.CurrentIteration;
ChunkStartIdx=handles.ChunkerizedDataIndices(CurrTarget,1,CurrIter)/handles.obj.SampleRate;
ChunkStopIdx=handles.ChunkerizedDataIndices(CurrTarget,2,CurrIter)/handles.obj.SampleRate;
if handles.StackPulses==1
    newlims=[0 handles.PostPulseStartSecondsToView];
else
    newlims=[ChunkStartIdx ChunkStopIdx];
end

function LoadImageFcn(obj,evt)
handles=guidata(obj);

% read in image file
[filename,pathname,FilterIndex]=uigetfile('*.*','Choose an image file');
if ~FilterIndex
    return
end
handles.Image=imread([pathname,filename]);

% if there is not an ImageFig already, set it up
if ~isfield(handles,'ImageFig')
    handles.ImageFig=figure('name',filename);
else
    figure(handles.ImageFig)
    set(gcf,'name',filename);
end
handles.ContrastSlider = ...
    uicontrol('Style','slider', ...
    'Units','normalized', ...
    'Position',[0.0  0.0 .4 .05], ...
    'Min',0, ...
    'Max',1, ...
    'Sliderstep',[.01 .05], ...
    'Value',1/3, ...
    'Enable','on', ...
    'Callback', @ContrastBrightnessFcn);
handles.BrightnessSlider = ...
    uicontrol('Style','slider', ...
    'Units','normalized', ...
    'Position',[0.5 0.0 0.4 .05], ...
    'Min',0, ...
    'Max',1, ...
    'Sliderstep',[.01 .05], ...
    'Value',1/3, ...
    'Enable','on', ...
    'Callback', @ContrastBrightnessFcn);

% display the image in the imageaxes
figure(handles.ImageFig);
imagesc(handles.Image);
% set(handles.ImageFig,'YTickLabel',[])
% set(handles.ImageFig,'XTickLabel',[])
colormap(gray);

% Send displayfig to ImageFig
guidata(handles.ImageFig,handles.displayfig);

guidata(handles.displayfig,handles);%pass data back to the figure

function LoadContoursFcn(obj,event)
handles=guidata(obj);
[filename,pathname,FilterIndex]=uigetfile('*.mat','Choose a contour file created by caltracer');
if FilterIndex
    handles.contourfile=[pathname,filename];
    loadedcontours=load(handles.contourfile);
    if isfield(loadedcontours,'CONTS');
        handles.contours=loadedcontours.CONTS;
    elseif isfield(loadedcontours,'outputcontours');
        handles.contours=loadedcontours.outputcontours;
    else
        msgbox('file does not have contours in expected format, see ephysviewer loadcontoursfcn');
    end
end


guidata(handles.displayfig,handles);

function ContrastBrightnessFcn(obj,evnt)
% Get displayfig from this fig, then get handles from displayfig
displayfig=guidata(obj);
handles=guidata(displayfig);
brightness = get(handles.BrightnessSlider,'value');
contrast = get(handles.ContrastSlider,'value');
contrast = contrast^0.25;
if contrast < 0.5
    m = 2*contrast/255;
else
    m = 1/(510-510*contrast+eps);
end
if brightness < 0.5
    b = -255*m*(1-(2*brightness)^0.25);
else
    b = 2*brightness-1;
end

cc = (0:255)*m+b;
cc(cc<0)=0;
cc(cc>1)=1;
climg = repmat(cc',1,3);
colormap(climg);

function CurrChunkEditFcn(obj,event)
handles=guidata(obj);
handles.CurrentTargetNumber=str2double(get(handles.CurrChunkEditBox,'String'));
if handles.CurrentTargetNumber>handles.NumTargets
    handles.CurrentTargetNumber=handles.NumTargets;
end
set(handles.CurrChunkEditBox,'string',num2str(handles.CurrentTargetNumber));
handles.xlims=ConvertChunkToXlims(handles);
handles.displaylength=handles.xlims(2)-handles.xlims(1);
set(handles.axeshandles,'xlim',handles.xlims);%set all xlims the same
set(handles.HorizMinEditBox,'String',num2str(handles.xlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(handles.xlims(2),3));
PlotMinWithinXlims(handles,get(handles.axeshandles(1),'XLim'));
handles=PlotTargetOnImage(handles);
handles=PlotFiredCell(handles);
if isfield(handles,'SynchTraceFig')
    handles=PlotSynchronousTrace(handles);
end
guidata(handles.displayfig,handles);%pass data back to the figure

function NextChunkFcn(obj,event)
handles=guidata(obj);
handles.CurrentTargetNumber=handles.CurrentTargetNumber+1;
if handles.CurrentTargetNumber>handles.NumTargets
    handles.CurrentTargetNumber=handles.NumTargets;
end
set(handles.CurrChunkEditBox,'string',num2str(handles.CurrentTargetNumber));
handles.xlims=ConvertChunkToXlims(handles);
handles.displaylength=handles.xlims(2)-handles.xlims(1);
set(handles.axeshandles,'xlim',handles.xlims);%set all xlims the same
set(handles.HorizMinEditBox,'String',num2str(handles.xlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(handles.xlims(2),3));
PlotMinWithinXlims(handles,get(handles.axeshandles(1),'XLim'));
handles=PlotTargetOnImage(handles);
handles=PlotFiredCell(handles);
if isfield(handles,'SynchTraceFig')
    handles=PlotSynchronousTrace(handles);
end
guidata(handles.displayfig,handles);%pass data back to the figure

function PrevChunkFcn(obj,event)
handles=guidata(obj);
handles.CurrentTargetNumber=handles.CurrentTargetNumber-1;
if handles.CurrentTargetNumber==0
    handles.CurrentTargetNumber=1;
end
set(handles.CurrChunkEditBox,'string',num2str(handles.CurrentTargetNumber));
handles.xlims=ConvertChunkToXlims(handles);
handles.displaylength=handles.xlims(2)-handles.xlims(1);
set(handles.axeshandles,'xlim',handles.xlims);%set all xlims the same
set(handles.HorizMinEditBox,'String',num2str(handles.xlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(handles.xlims(2),3));
PlotMinWithinXlims(handles,get(handles.axeshandles(1),'XLim'));
handles=PlotTargetOnImage(handles);
handles=PlotFiredCell(handles);
if isfield(handles,'SynchTraceFig')
    handles=PlotSynchronousTrace(handles);
end
guidata(handles.displayfig,handles);%pass data back to the figure

function CurrIterEditFcn(obj,event)
handles=guidata(obj);
handles.CurrentIteration=str2double(get(handles.CurrIterEditBox,'String'));
if handles.CurrentIteration>handles.NumIterations
    handles.CurrentIteration=handles.NumIterations;
end
set(handles.CurrIterEditBox,'string',num2str(handles.CurrentIteration));
handles.xlims=ConvertChunkToXlims(handles);
handles.displaylength=handles.xlims(2)-handles.xlims(1);
set(handles.axeshandles,'xlim',handles.xlims);%set all xlims the same
set(handles.HorizMinEditBox,'String',num2str(handles.xlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(handles.xlims(2),3));
PlotMinWithinXlims(handles,get(handles.axeshandles(1),'XLim'));
handles=PlotTargetOnImage(handles);
handles=PlotFiredCell(handles);
if isfield(handles,'SynchTraceFig')
    handles=PlotSynchronousTrace(handles);
end
guidata(handles.displayfig,handles);%pass data back to the figure

function NextIterFcn(obj,event)
handles=guidata(obj);
handles.CurrentIteration=handles.CurrentIteration+1;
if handles.CurrentIteration>handles.NumIterations
    handles.CurrentIteration=handles.NumIterations;
end
set(handles.CurrIterEditBox,'string',num2str(handles.CurrentIteration));
handles.xlims=ConvertChunkToXlims(handles);
handles.displaylength=handles.xlims(2)-handles.xlims(1);
set(handles.axeshandles,'xlim',handles.xlims);%set all xlims the same
set(handles.HorizMinEditBox,'String',num2str(handles.xlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(handles.xlims(2),3));
PlotMinWithinXlims(handles,get(handles.axeshandles(1),'XLim'));
handles=PlotTargetOnImage(handles);
handles=PlotFiredCell(handles);
if isfield(handles,'SynchTraceFig')
    handles=PlotSynchronousTrace(handles);
end
guidata(handles.displayfig,handles);%pass data back to the figure

function PrevIterFcn(obj,event)
handles=guidata(obj);
handles.CurrentIteration=handles.CurrentIteration-1;
if handles.CurrentIteration==0
    handles.CurrentIteration=1;
end
set(handles.CurrIterEditBox,'string',num2str(handles.CurrentIteration));
handles.xlims=ConvertChunkToXlims(handles);
handles.displaylength=handles.xlims(2)-handles.xlims(1);
set(handles.axeshandles,'xlim',handles.xlims);%set all xlims the same
set(handles.HorizMinEditBox,'String',num2str(handles.xlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(handles.xlims(2),3));
PlotMinWithinXlims(handles,get(handles.axeshandles(1),'XLim'));
handles=PlotTargetOnImage(handles);
handles=PlotFiredCell(handles);
if isfield(handles,'SynchTraceFig')
    handles=PlotSynchronousTrace(handles);
end
guidata(handles.displayfig,handles);%pass data back to the figure

function handles=PlotTargetOnImage(handles)
if isfield(handles,'ImageFig')
    CurrTarget=handles.Reindexer(handles.CurrentTargetNumber,handles.CurrentIteration);
    CurrIter=handles.CurrentIteration;
    % Deal with VNTs that rely on 'Fire Cycles'
    % i.e. those without iterations
    if length(handles.targetXidxs)==handles.NumTargets
        OrigChunk=CurrTarget;
    else
        OrigChunk=handles.NumTargets*(CurrIter-1)+CurrTarget;
    end
    if isfield(handles,'VNTIdxsOfMatches')
        TargetX=handles.targetXidxs(OrigChunk);
        TargetY=handles.targetYidxs(OrigChunk);
        title4=['Target coordinates are ' num2str(TargetX) ', ' num2str(TargetY)];
        if isfield(handles,'TargetHandle');
            try
                % comment out next line to accrue targets on figure
                delete(handles.TargetHandle);
            catch
            end
        end
        figure(handles.ImageFig);
        hold on;
        if isfield(handles,'contours');
            CurrContour=handles.contours{CurrTarget};
            handles.TargetHandle=patch(CurrContour(:,1),CurrContour(:,2),[0 0 1]);
        else
            handles.TargetHandle=plot(TargetX,TargetY,'Marker','o','MarkerSize',12,'Color','b');
        end
        figure(handles.displayfig);
    end
else
    title4=[];
end

function PlotEntireVNTOnImage(obj,evt)
handles=guidata(obj);

if isfield(handles,'ImageFig')
    if isfield(handles,'VNTIdxsOfMatches')
        TargetX=handles.targetXidxs;
        TargetY=handles.targetYidxs;
        if isfield(handles,'TargetHandle');
            try
                delete(handles.TargetHandle);
            catch
            end
        end
        figure(handles.ImageFig);
        hold on;
        handles.TargetHandle=plot(TargetX,TargetY,'d','Marker','o','MarkerSize',12,'Color','b');
        figure(handles.displayfig);
    end
end
guidata(handles.displayfig,handles)

function LocalMenuExportVisibleDataFcn(obj,event)
%function LocalMenuExportDataFcn(obj,event);
%This function allows the user to export the data (actual numbers)
%contained in the viewed data (data within the current xlimits).  Data can
%either be exported to the base workspace, to a .mat file or to a binary
%file, which has the same precision as the data was originally encoded in.

handles=guidata(obj);%get general info from the figure
data=getappdata(handles.displayfig,'decdata');%get saved loaded data
if isempty(data)
    helpdlg('Data must be loaded before exporting')
    return
end

charcell={'Base Workspace','Base Workspace (mem efficient)','.Mat File','Binary File'};%set up a list of options for export
[selection,ok]=listdlg('ListString',charcell,'SelectionMode','single','Name','Export List','PromptString','Export to where?');
if isempty(ok) || ~ok;%if user hit cancel or closed the window
    return%exit the export function
end

data=data.data{1};%keep just the full data (not the decimated versions)
xlims=get(handles.axeshandles(1),'xlim');%get the width of the current viewing window

yindices=round(xlims(1)*handles.obj.SampleRate):round(xlims(2)*handles.obj.SampleRate);
if yindices(1)<1;
    yindices(1)=[];
end
if yindices(end)>handles.obj.SamplesAcquired/handles.obj.SampleRate;
    yindices(end)=[];
end
yindices=yindices';
data=data(yindices,:);
% pause(1);
% xpoints=yindices./handles.obj.SampleRate;
yindices=yindices./handles.obj.SampleRate;

switch selection%depending on what the user chose to do
    case 1 %'Workspace'
        if isempty(str2num(handles.filename(1)))
            structname=handles.filename;
        else
            structname=strcat('x',strrep(handles.filename,'-','_'));
        end
        for a=1:size(data,2);%also corresponds to handles.includedchannels
            eval([structname,...
                '.chan',num2str(handles.obj.Channel(handles.includedchannels(a)).HwChannel),...
                '_',handles.obj.Channel(handles.includedchannels(a)).ChannelName,...
                '=data(:,a);']);
        end
        eval([structname,'.time=yindices;']);
        assignin('base',structname,eval(structname));%export to base workspace
    case 2 %'Base Workspace (mem efficient)'
        %         data=data(yindices,:);
        assignin('base','Time',yindices);
        assignin('base','ChannelNames',{handles.obj.Channel(handles.includedchannels).ChannelName});
        assignin('base','DataFromViewer',data);
    case 3 %'.Mat File'
        % 		data=data(yindices,:);
        for a=1:size(data,2);%also corresponds to handles.includedchannels
            eval(['ExpDat.',...
                'chan',num2str(handles.obj.Channel(handles.includedchannels(a)).HwChannel),...
                '_',handles.obj.Channel(handles.includedchannels(a)).ChannelName,...
                '=data(:,a);']);
        end
        %         ExpDat.time=xpoints;
        ExpDat.time=yindices;
        try
            [pathstr,name,ext] = fileparts(handles.obj.LogFileName);%get name of the .daq file from which data was loaded
            [FileName,PathName,FilterIndex] = uiputfile('.mat','Save Data to Disk',[pathstr,'\']);%ask user for name... above.mat is default
        catch
            [FileName,PathName,FilterIndex] = uiputfile('.mat','Save Data to Disk');%ask user for name... above.mat is default
        end
        if FilterIndex%if user did not hit cancel
            savefile=[PathName,FileName];
            save(savefile,'ExpDat');%save to mat file
        else
            return
        end
    case 4 %'Binary File'
        msgbox('Not yet implemented','Binary output','Help')
        
        %         if strcmp(class(handles.obj),'analoginput');%depending on whether this is being called from RecordDaq or DaqViewer
        %     		datatype=daqhwinfo(handles.obj);%get hwinfo
        % 			datatype=datatype.NativeDataType;%ie int16
        %         else
        %             datatype=handles.obj.HwInfo.NativeDataType;%ie int16,uint16, etc
        %         end
        %         clear data%dump data... we don't need it any more since we'll read it from the file in native format
        %         yindices=[yindices(1) yindices(end)];%this is all we'll need for daqread, this way is better for ram if big stuff later
        %         [pathstring,name,ext,versn] = fileparts(handles.path);
        %         for a=1:length(handles.obj.Channel);%for every hwchannel
        %             hwchans(a)=handles.obj.Channel(a).HwChannel;%record it's number
        %         end
        %         for a=1:length(handles.includedchannels);%for each channel
        %             thishwchannel=hwchans(handles.includedchannels(a));%record the number of the specified hwchannel
        %             totalpathstr=[pathstring,'\',name,'_chan',num2str(thishwchannel),'_',...%make a path for a new file to be saved for each channel
        %                     handles.obj.Channel(handles.includedchannels(a)).ChannelName,'_',datatype,'.bin'];%including it's HW channel number and name
        % 			[FileName,PathName,FilterIndex] = uiputfile('.bin','Save This Channel to Disk?',totalpathstr);
        %             if FilterIndex%if user did not hit cancel
        %                 data=bdaqread(handles.path,'DataFormat','Native',...%read data only for specified channels for specified samples
        %                     'Channels',handles.includedchannels(a),'Samples',yindices);
        %                 savefile=[PathName,FileName];
        %                 [trash1,trash2,ext,versn] = fileparts(savefile);
        %                 if isempty(ext);%if no extension;
        %                     savefile=[savefile,'.bin'];%add one
        %                 end
        %                 fid=fopen(savefile,'w');%create a file with the user-specified name
        %                 eval(['fwrite(fid,data,''',datatype,''');']);%save data into it
        %                 fclose(fid);
        %             else
        %                 continue
        %             end
        %         end
        %         totalpathstr=[pathstring,'\',name,'_times_double.bin'];%including it's HW channel number and name
        % 		[FileName,PathName,FilterIndex] = uiputfile('.bin','Save This Channel to Disk?',totalpathstr);
        %         if FilterIndex%if user did not hit cancel
        %             savefile=[PathName,FileName];
        %             [trash1,trash2,ext,versn] = fileparts(savefile);
        %             if isempty(ext);%if no extension;
        %                 savefile=[savefile,'.bin'];%add one
        %             end
        %             fid=fopen(savefile,'w');%create a file with the user-specified name
        %             fwrite(fid,xpoints,'double');%save data into it
        %             fclose(fid);
        %         end
end

function LocalMenuGridFcn(obj,event)

handles=guidata(obj);%get general info from the figure

if strcmp(get(handles.axeshandles(1),'XGrid'),'off');
    set(handles.axeshandles,'XGrid','on','YGrid','on')
else
    set(handles.axeshandles,'XGrid','off','YGrid','off')
end

function CloseViewerFcn(obj,event)

delete(obj)%deleting the figure

function HorizScaleInButtonFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
xlimits=get(handles.axeshandles(1),'XLim');%in seconds
power10=10^floor(log10(handles.displaylength));
m=handles.displaylength/power10;
value=max(find(handles.HorizScaleLookup<m));%find indices of the preset largest value that is less than the current display length
value=handles.HorizScaleLookup(value)*power10;%seconds

if value>handles.sizedata(1)/handles.obj.SampleRate%if wants to zoom out to wider than total num of points
    value = handles.sizedata(1)/handles.obj.SampleRate;
end

if value<5/handles.obj.SampleRate;%don't allow to zoom too much
    value=5/handles.obj.SampleRate;
end
meanpos=mean(xlimits);
newlims=[meanpos-.5*value meanpos+.5*value];
set(handles.HorizMinEditBox,'String',num2str(newlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(newlims(2),3));
set(handles.axeshandles,'XLim',newlims);

PlotMinWithinXlims(handles,newlims)
if isfield(handles,'SynchTraceFig')
    handles=PlotSynchronousTrace(handles);
end

handles.displaylength=value;
% set(handles.HorizScaleTextBox,'String',num2str(handles.displaylength));%set value of text box in seconds, and set it late
%     %in the function so it represents the most recent information possible

guidata(handles.displayfig,handles);%pass data back to the figure

function HorizScaleOutButtonFcn(obj,ev);

handles=guidata(obj);%get general info from the figure

xlimits=get(handles.axeshandles(1),'XLim');
power10=10^floor(log10(handles.displaylength));
m=handles.displaylength/power10;
value=min(find(handles.HorizScaleLookup>m));%find indices of the preset largest value that is less than the current display length
value=handles.HorizScaleLookup(value)*power10;
meanpos=mean(xlimits);
newlims=[meanpos-.5*value meanpos+.5*value];

defaultstart=1/handles.obj.SampleRate;
% defaultend=handles.obj.SamplesAcquired/handles.obj.SampleRate;
defaultend=handles.sizedata(1)/handles.obj.SampleRate;%total seconds
if value>defaultend%both the left and right ends will be over the limits
    newlims=[defaultstart defaultend];%fix that
    value=handles.obj.SamplesAcquired/handles.obj.SampleRate;%save value correctly so later handles.displaylength will be correct
elseif newlims(1)<defaultstart;%if only the left side is too low (but right is fine)
    newlims=[defaultstart value];%plot from zero to the required length
elseif newlims(2)>defaultend;%if only the right side is too high
    newlims=[defaultend-value defaultend];%plot required length back from the end of the data
end
set(handles.HorizMinEditBox,'String',num2str(newlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(newlims(2),3));
set(handles.axeshandles,'XLim',newlims);

PlotMinWithinXlims(handles,newlims)
if isfield(handles,'SynchTraceFig')
    handles=PlotSynchronousTrace(handles);
end

handles.displaylength=value;
% set(handles.HorizScaleTextBox,'String',num2str(handles.displaylength));%set value of text box in seconds, and set it late
%in the function so it represents the most recent information possible

guidata(handles.displayfig,handles);%pass data back to the figure

function MoveLeftFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
% set(handles.obj,'TimerFcn','');%pause plotting (resume later)
% pause(.1)
xlimits=get(handles.axeshandles(1),'XLim');
dl=xlimits(2)-xlimits(1);

newlims=[xlimits(1)-.5*dl xlimits(2)-.5*dl];
if newlims(1)<0;
    newlims=[1/handles.obj.SampleRate 1/handles.obj.SampleRate+dl];
end
set(handles.axeshandles,'XLim',newlims);
set(handles.HorizMinEditBox,'String',num2str(newlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(newlims(2),3));
PlotMinWithinXlims(handles,newlims)
if isfield(handles,'SynchTraceFig')
    handles=PlotSynchronousTrace(handles);
end

function MoveRightFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
% set(handles.obj,'TimerFcn','');%pause plotting (resume later)
% pause(.1)
xlimits=get(handles.axeshandles(1),'XLim');
dl=xlimits(2)-xlimits(1);

newlims=[xlimits(1)+.5*dl xlimits(2)+.5*dl];
defaultend=handles.sizedata(1)/handles.obj.SampleRate;%total seconds
if newlims(2)>defaultend;%if right side is going beyond the end of the data
    newlims=[defaultend-dl defaultend];
end
set(handles.axeshandles,'XLim',newlims);
set(handles.HorizMinEditBox,'String',num2str(newlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(newlims(2),3));
PlotMinWithinXlims(handles,newlims)
if isfield(handles,'SynchTraceFig')
    handles=PlotSynchronousTrace(handles);
end

function VertSlideUpFcn(obj,ev)

handles=guidata(obj);%get general info from the figure

ind=find(handles.VertSlideUp==obj);%find number of axis to be changed

ylims=get(handles.axeshandles(ind),'Ylim');
range=ylims(2)-ylims(1);
newylims=ylims+.3*range;
maxy=handles.maxy(ind);
if max(newylims)>maxy;
    newylims=[maxy-range maxy];
end

set(handles.axeshandles(ind),'Ylim',newylims);%slide up 1/4 field of view
set(handles.VertMinEditBox(ind),'String',num2str(newylims(1),3));
set(handles.VertMaxEditBox(ind),'String',num2str(newylims(2),3));

function VertSlideDownFcn(obj,ev)

handles=guidata(obj);%get general info from the figure

ind=find(handles.VertSlideDown==obj);%find number of axis to be changed
ylims=get(handles.axeshandles(ind),'Ylim');
range=ylims(2)-ylims(1);
newylims=ylims-.3*range;

miny=handles.miny(ind);
if min(newylims)<miny;
    newylims=[miny miny+range];
end

set(handles.axeshandles(ind),'Ylim',newylims);%slide up 1/4 field of view
set(handles.VertMinEditBox(ind),'String',num2str(newylims(1),3));
set(handles.VertMaxEditBox(ind),'String',num2str(newylims(2),3));

function VertScaleInFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
ind=find(handles.VertScaleIn==obj);%find number of axis to be changed

ylims=get(handles.axeshandles(ind),'Ylim');
m=mean(ylims);
siderange=m-ylims(1);
newylims=[m-(1/handles.vertzoomfactor)*siderange m+(1/handles.vertzoomfactor)*siderange];
set(handles.axeshandles(ind),'Ylim',newylims);%slide up 1/4 field of view
set(handles.VertMinEditBox(ind),'String',num2str(newylims(1),3));
set(handles.VertMaxEditBox(ind),'String',num2str(newylims(2),3));

function VertScaleOutFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
if isempty(ev)
    ind=find(handles.VertScaleOut==obj);%find number of axis to be changed
else
    ind=find(handles.axeshandles==obj);%find the axis to be changed
end

ylims=get(handles.axeshandles(ind),'Ylim');
m=mean(ylims);
siderange=m-ylims(1);
newylims=[m-(handles.vertzoomfactor)*siderange m+(handles.vertzoomfactor)*siderange];

miny=handles.miny(ind);%use the index number found above
maxy=handles.maxy(ind);
if siderange*2*handles.vertzoomfactor>maxy-miny;
    newylims=[miny maxy];
elseif newylims(1)<miny && newylims(2)<=maxy;%if top side is greater than allowable, but bottom is ok
    newylims=[miny miny+siderange*2*handles.vertzoomfactor];%keep same range, but set top to maxy
elseif newylims(2)>maxy && newylims(1)>=miny;%if bottom side is less than allowable, but top is ok
    newylims=[maxy-siderange*2*handles.vertzoomfactor maxy];%keep same range, but set bottom to miny
end

set(handles.axeshandles(ind),'Ylim',newylims);%slide up 1/4 field of view
set(handles.VertMinEditBox(ind),'String',num2str(newylims(1),3));
set(handles.VertMaxEditBox(ind),'String',num2str(newylims(2),3));

function VertMinSetFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
ind=find(handles.VertMinEditBox==obj);%find number of axis to be changed
newmin=str2double(get(handles.VertMinEditBox(ind),'String'));%get user input
if newmin<handles.miny(ind);%if value was beyond possible input values
    newmin=handles.miny(ind);%set to the limit input value
    set(handles.VertMinEditBox(ind),'String',num2str(newmin,3));
end
newylims=get(handles.axeshandles(ind),'Ylim');%get old limits
newylims(1)=newmin;%sub in new value
set(handles.axeshandles(ind),'Ylim',newylims);%set axes

function VertMaxSetFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
ind=find(handles.VertMaxEditBox==obj);%find number of axis to be changed
newmax=str2double(get(handles.VertMaxEditBox(ind),'String'));%get user input
if newmax>handles.maxy(ind);%if value was beyond possible input values
    newmax=handles.maxy(ind);%set to the limit input value
    set(handles.VertMaxEditBox(ind),'String',num2str(newmax,3));
end
newylims=get(handles.axeshandles(ind),'Ylim');%get old limits
newylims(2)=newmax;%sub in new value
set(handles.axeshandles(ind),'Ylim',newylims);%set axes

function HorizMinSetFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
xlimits=get(handles.axeshandles(1),'XLim');%in seconds

newlims = xlimits;
newlims(1) = str2double(get(obj,'string'));
if newlims(1)<0;
    newlims(1)=1/handles.obj.SampleRate;
end
if newlims(2)<=newlims(1)
    return
end
set(handles.axeshandles,'XLim',newlims);
PlotMinWithinXlims(handles,newlims)
if isfield(handles,'SynchTraceFig')
    handles=PlotSynchronousTrace(handles);
end

function HorizMaxSetFcn(obj,ev)

handles=guidata(obj);%get general info from the figure
xlimits=get(handles.axeshandles(1),'XLim');%in seconds

newlims = xlimits;
newlims(2) = str2double(get(obj,'string'));
defaultend=handles.sizedata(1)/handles.obj.SampleRate;%total seconds
if newlims(2)>defaultend;%if right side is going beyond the end of the data
    newlims(2) = defaultend;
end
if newlims(2)<=newlims(1)
    return
end
set(handles.axeshandles,'XLim',newlims);
PlotMinWithinXlims(handles,newlims)
if isfield(handles,'SynchTraceFig')
    handles=PlotSynchronousTrace(handles);
end

function ZoomFcn(obj,ev)

handles=guidata(obj);

point1 = get(handles.displayfig,'CurrentPoint'); % button down detected
if strcmp(get(handles.displayfig,'SelectionType'),'normal');%if a left click
    A=ver;
    LaterThanR2014b=any(strcmp(A(1).Release,{'(R2014b)' '(R2015a Prerelease)'}));
    if ismac
     if LaterThanR2014b
         useimrect=1;
     else
         useimrect=0;
     end
    else
        useimrect=0;
    end
    useimrect=0;
    if useimrect;
        imr=imrect;
        r2=getPosition(imr);
        xlims(1)=r2(1);
        xlims(2)=r2(1)+r2(3);
        ylims(1)=r2(2);
        ylims(2)=r2(2)+r2(4);
    else
        if ~LaterThanR2014b
            rect = [0 0 point1(1,1) point1(1,2)];
            r2 = rbbox(rect,[point1(1,1) point1(1,2)]); %in figure units... need to convert to x and y coords... mult by axes relative position and by xlim and ylim of ax
        else
            r2=rbbox;
        end
        axpos=get(obj,'position');
        oldx=get(obj,'xlim');
        oldy=get(obj,'ylim');
        xlims(1)=(r2(1)-axpos(1))/axpos(3);%convert first point to x start in axes normalized
        xlims(2)=r2(3)/axpos(3)+xlims(1);%convert third point to x stop in axes normalized
        xlims(1)=oldx(1)+xlims(1)*(oldx(2)-oldx(1));
        xlims(2)=oldx(1)+xlims(2)*(oldx(2)-oldx(1));
        ylims(1)=(r2(2)-axpos(2))/axpos(4);
        ylims(2)=r2(4)/axpos(4)+ylims(1);
        ylims(1)=oldy(1)+ylims(1)*(oldy(2)-oldy(1));
        ylims(2)=oldy(1)+ylims(2)*(oldy(2)-oldy(1));
    end
    if (r2(3)+r2(4))>=.01;%if at least some minimum box was made
        handles.displaylength=xlims(2)-xlims(1);
        set(handles.axeshandles,'xlim',xlims);%set all xlims the same
        set(obj,'ylim',ylims);
        PlotMinWithinXlims(handles,xlims)
        if isfield(handles,'SynchTraceFig')
            handles=PlotSynchronousTrace(handles);
        end
        
        handles.displaylength=(xlims(2)-xlims(1));%
        
        %set x and y text boxes
        set(handles.HorizMinEditBox,'String',num2str(xlims(1),3));
        set(handles.HorizMaxEditBox,'String',num2str(xlims(2),3));
        ind = find(handles.axeshandles==obj);%find number of axis to be changed
        set(handles.VertMinEditBox(ind),'String',num2str(ylims(1),3));%set text boxes
        set(handles.VertMaxEditBox(ind),'String',num2str(ylims(2),3));
        
        % 		set(handles.HorizScaleTextBox,'String',num2str(handles.displaylength));%set value of text box in seconds, and set it late
        %in the function so it represents the most recent information possible
        guidata(handles.displayfig,handles);%pass data back to the figure
        
    end
elseif strcmp(get(handles.displayfig,'SelectionType'),'alt');%if a left click
    HorizScaleOutButtonFcn(obj,'notempty');
    VertScaleOutFcn(obj,'notempty');
end

function ResetAxesButtonFcn(obj,ev)
handles=guidata(obj);
if handles.ChunkMode==1;
    if handles.StackPulses==1;
        handles.xlims=[0 handles.PostPulseStartSecondsToView];
    else
        handles.xlims=ConvertChunkToXlims(handles);
    end
    handles.displaylength=handles.xlims(2)-handles.xlims(1);
else
    handles.displaylength=handles.sizedata(1)/handles.obj.SampleRate;%total seconds
    handles.xlims=[1/handles.obj.SampleRate handles.displaylength];
end

set(handles.HorizMinEditBox,'String',num2str(handles.xlims(1),3));
set(handles.HorizMaxEditBox,'String',num2str(handles.xlims(2),3));
set(handles.axeshandles,'xlim',[handles.xlims(1) handles.xlims(2)]);

for a=1:length(handles.axeshandles)
    set(handles.axeshandles(a),'ylim',[handles.miny(a) handles.maxy(a)])
    set(handles.VertMinEditBox(a),'String',num2str(handles.miny(a),3));
    set(handles.VertMaxEditBox(a),'String',num2str(handles.maxy(a),3));
end

PlotMinWithinXlims(handles,[1/handles.obj.SampleRate handles.displaylength])
if isfield(handles,'SynchTraceFig')
    handles=PlotSynchronousTrace(handles);
end
guidata(handles.displayfig,handles)

function SingleMeasureFcn(obj,ev)

handles = guidata(obj);
measuremode = handles.SingleMeasureMode;
[x,y,channelstomeasure] = MeasurementKernelFcn(handles,measuremode);
if ~isempty(handles.SingleMeasuresMatrix.MeasureNumber)
    measnum = handles.SingleMeasuresMatrix.MeasureNumber(end)+1;
else
    measnum = 1;
end
for chidx = 1:size(y,2);
    handles.SingleMeasuresMatrix.MeasureNumber(end+1) = measnum;
    handles.SingleMeasuresMatrix.FileName{end+1} = handles.filename;
    cnn = handles.includedchannels(channelstomeasure(chidx));
    chname = handles.obj.Channel(cnn).ChannelName;
    handles.SingleMeasuresMatrix.ChannelName{end+1} = chname;
    handles.SingleMeasuresMatrix.xs(end+1) = x;
    handles.SingleMeasuresMatrix.ys(end+1) = y(chidx);
    handles.SingleMeasuresMatrix.MeasureMode{end+1} = measuremode;
    handles.SingleMeasuresMatrix.Comments{end+1} = '';
end

handles = SingleMeasuresBox(handles);

guidata(handles.displayfig,handles);

function [x,y,channelstomeasure] = MeasurementKernelFcn(handles,measuremode);
if ~strcmp('Click',measuremode)
    [tempx,tempy] = ginput(1);
    point1 = get(handles.displayfig,'CurrentPoint'); % button down detected
    if strcmp(get(handles.displayfig,'SelectionType'),'normal');%if a left click
        rect = [point1(1,1) point1(1,2) 0 0];
        r2 = rbbox(rect);%in figure units = seconds
        thisaxishand = gca;
        axpos = get(thisaxishand,'position');
        axxlims = get(thisaxishand,'xlim');
        xlims(1)=(r2(1)-axpos(1))/axpos(3);%convert first point to x start in axes normalized
        xlims(2)=r2(3)/axpos(3)+xlims(1);%convert third point to x stop in axes normalized
        
        
        xstartsec = axxlims(1) + xlims(1)*(axxlims(2)-axxlims(1));
        xstopsec = axxlims(1) + xlims(2)*(axxlims(2)-axxlims(1));
        %         if (r2(3)+r2(4))>=.001;%if at least some minimum box was made
        %
        %get the number of the channel in the axes clicked.  Will
        %defnitely get measures from this channel
        thisaxisidx = find(handles.axeshandles == thisaxishand);%grab the number of the axes clicked...
        channelsthisaxes=find(handles.ChannelAxesMatrix(thisaxisidx,:));
        %             channelsthisaxes=find(handles.obj.userdata.channelaxes==thisaxisidx);%get the channel in the axes...
        %will be only one channel in this program
        %find which channels had boxes checked indicating to measure them
        for a = 1:length(handles.ChannelMeasureCheckbox);
            axeschecked(a) = get(handles.ChannelMeasureCheckbox(a),'value');
        end
        axeschecked=find(axeschecked);
        %     The following command needs to be debugged for the instance of
        %     plotting more than one channel per axes, which I've never fully
        %     fixed. Until then, I'm just going to assume the axis that is checked
        %     corresponds directly to the one channel that should be there.  AP 20090302
        %           channelschecked=handles.ChannelAxesMatrix(find(handles.ChannelAxesMatrix(axeschecked,:)'));
        %make final list of all channels to measure from channel
        %clicked on and those with boxes checked
        channelstomeasure = union(channelsthisaxes,axeschecked);%auto sorted in ascending order by union fcn
        
        decdata=getappdata(handles.displayfig,'decdata');
        channeldata=decdata.data{1}(:,channelstomeasure);
        clear decdata;%decdata is big, get rid of it
        
        yindices=round(xstartsec*handles.obj.SampleRate):round(xstopsec*handles.obj.SampleRate);
        yindices(yindices>length(channeldata))=[];
        yindices(yindices<1)=[];
        if isempty(yindices);
            return
        end
        xdata=yindices/handles.obj.SampleRate;
        ydata = channeldata(yindices,:);
        clear channeldata
        %for max and min: needs to be max/min of clicked channel and same time points on
        %other channels.  For mean: take mean across selected region on
        %each channel.
        %Then give just one x output
        if strcmp('Minimum',measuremode)
            selectedchidx = find(channelstomeasure == channelsthisaxes);
            [trash,x] = min(ydata(:,selectedchidx),[],1);%get x&y of min in the CLICKED ON channel
            %use x to get all the y's
            y = ydata(x,:);
            x = xdata(x);
        elseif strcmp('Maximum',measuremode)
            selectedchidx = find(channelstomeasure == channelsthisaxes);
            [trash,x] = max(ydata(:,selectedchidx),[],1);%get x&y of max in the CLICKED ON channel
            %use x to get all the y's
            y = ydata(x,:);
            x = xdata(x);
        elseif strcmp('Mean',measuremode)
            y = mean(ydata,1);
            x = mean(xdata);
        end
        %         end
    end
else
    [x,trash] = ginput(1);
    
    thisaxishand = gca;
    thisaxisidx = find(handles.axeshandles == thisaxishand);%grab the number of the axes clicked...
    channelsthisaxes=find(handles.ChannelAxesMatrix(thisaxisidx,:));
    %     channelsthisaxes=find(handles.obj.userdata.channelaxes==thisaxisidx);%get the channel in the axes...
    %will be only one channel in this program
    %find which channels had boxes checked indicating to measure them
    for a = 1:length(handles.ChannelMeasureCheckbox);
        axeschecked(a) = get(handles.ChannelMeasureCheckbox(a),'value');
    end
    axeschecked=find(axeschecked);
    %     The following command needs to be debugged for the instance of
    %     plotting more than one channel per axes, which I've never fully
    %     fixed. Until then, I'm just going to assume the axis that is checked
    %     corresponds directly to the one channel that should be there.  AP 20090302
    %   channelschecked=handles.ChannelAxesMatrix(find(handles.ChannelAxesMatrix(axeschecked,:)'));
    
    %make final list of all channels to measure from channel
    %clicked on and those with boxes checked
    channelstomeasure = union(channelsthisaxes,axeschecked);%auto sorted in ascending order by union fcn
    
    decdata=getappdata(handles.displayfig,'decdata');
    channeldata=decdata.data{1}(:,channelstomeasure);
    clear decdata;%decdata is big, get rid of it
    
    yindices=round(x*handles.obj.SampleRate);
    if yindices > length(channeldata) || yindices<1
        return
    end
    x = yindices/handles.obj.SampleRate;
    y = channeldata(yindices,:);
end

function handles = SingleMeasuresBox(handles);
%makes an auxiliary figure to display the Single Measures Data
thisfig = findobj('type','figure','name','Single Measures','userdata',handles.displayfig);
if isempty(thisfig);
    handles.SingleMeasuresWindow = figure('units','points',...
        'position',[0 40 325 1],...
        'toolbar','none',...
        'menubar','none',...
        'name','Single Measures',...
        'numbertitle','off',...
        'userdata',handles.displayfig);
else
    handles.SingleMeasuresWindow = thisfig;
end

delete(get(handles.SingleMeasuresWindow,'children'));

handles.SingleMeasuresExportMenu = uimenu('Label','Export','Parent',handles.SingleMeasuresWindow);
uimenu(handles.SingleMeasuresExportMenu,'Label','To Workspace','callback',@ExportSingleMeasuresToWS);
uimenu(handles.SingleMeasuresExportMenu,'Label','To .mat File','callback',@ExportSingleMeasuresToMAT);
uimenu(handles.SingleMeasuresExportMenu,'Label','To .xls File','callback',@ExportSingleMeasuresToXLS);
uimenu(handles.SingleMeasuresExportMenu,'Label','To current Excel Doc','callback',@ExportSingleMeasuresToCurrentExcelBook);
handles.SingleMeasuresDeleteMenu = uimenu('Label','Delete','Parent',handles.SingleMeasuresWindow);
uimenu(handles.SingleMeasuresDeleteMenu,'Label','All','callback',@DeleteSingleMeasuresAll);
uimenu(handles.SingleMeasuresDeleteMenu,'Label','Select','callback',@DeleteSingleMeasuresSelect);

precision = 8;
lineheight = 12;
fontsize = 8;
set(handles.SingleMeasuresWindow,'units','points');
nummeas = length(handles.SingleMeasuresMatrix.xs);
winsz = get(handles.SingleMeasuresWindow,'position');
winsz = [winsz(1) max([0 winsz(2)-lineheight]) winsz(3) (nummeas+2)*lineheight];
set(handles.SingleMeasuresWindow,'position',winsz);

%Column titles
uicontrol('parent',handles.SingleMeasuresWindow,'style','Text','units','points',...
    'position',[10 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','#');
uicontrol('parent',handles.SingleMeasuresWindow,'style','Text','units','points',...
    'position',[35 winsz(4)-12 55 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Channel');
uicontrol('parent',handles.SingleMeasuresWindow,'style','Text','units','points',...
    'position',[120 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','X');
uicontrol('parent',handles.SingleMeasuresWindow,'style','Text','units','points',...
    'position',[167 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Y');
uicontrol('parent',handles.SingleMeasuresWindow,'style','Text','units','points',...
    'position',[227 winsz(4)-12 55 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Comment');

nummeas = length(handles.SingleMeasuresMatrix.xs);
for midx = 1:nummeas;
    uicontrol('parent',handles.SingleMeasuresWindow,'style','frame','units','points',...
        'position',[7  winsz(4)-((lineheight-1)*(midx+1))-1 20 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.SingleMeasuresWindow,'style','text','units','points',...
        'position',[8  winsz(4)-((lineheight-1)*(midx+1)) 18 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.SingleMeasuresMatrix.MeasureNumber(midx)));
    
    uicontrol('parent',handles.SingleMeasuresWindow,'style','frame','units','points',...
        'position',[26  winsz(4)-((lineheight-1)*(midx+1))-1 75 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.SingleMeasuresWindow,'style','text','units','points',...
        'position',[27  winsz(4)-((lineheight-1)*(midx+1)) 73 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',handles.SingleMeasuresMatrix.ChannelName{midx});
    
    uicontrol('parent',handles.SingleMeasuresWindow,'style','frame','units','points',...
        'position',[101  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.SingleMeasuresWindow,'style','text','units','points',...
        'position',[102  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.SingleMeasuresMatrix.xs(midx),precision));
    
    uicontrol('parent',handles.SingleMeasuresWindow,'style','frame','units','points',...
        'position',[150  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.SingleMeasuresWindow,'style','text','units','points',...
        'position',[151  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.SingleMeasuresMatrix.ys(midx),precision));
    
    uicontrol('parent',handles.SingleMeasuresWindow,'style','frame','units','points',...
        'position',[199  winsz(4)-((lineheight-1)*(midx+1))-1 119 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.SingleMeasuresWindow,'style','edit','units','points',...
        'position',[201  winsz(4)-((lineheight-1)*(midx+1)) 117 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',handles.SingleMeasuresMatrix.Comments{midx},...
        'tag',num2str(midx),'userdata',handles.displayfig,'callback',@SingleMeasuresCommentFcn);
end

function SingleMeasuresCommentFcn(obj,ev);
EphysFig = get(obj,'userdata');
handles = guidata(EphysFig);
midx = str2double(get(obj,'tag'));
str = get(obj,'string');
handles.SingleMeasuresMatrix.Comments{midx}= str;
guidata(handles.displayfig,handles);

function PairedMeasureFcn(obj,ev)
% take handles.PairedMeasureMode and split with the _&_
handles = guidata(obj);
measuremode = handles.PairedMeasureMode;
andspot = strfind(measuremode,' & ');
measuremode1 = measuremode(1:andspot-1);
measuremode2 = measuremode(andspot+3:end);
[x1,y1,channelstomeasure] = MeasurementKernelFcn(handles, measuremode1);
[x2,y2,channelstomeasure] = MeasurementKernelFcn(handles, measuremode2);

if ~isempty(handles.PairedMeasuresMatrix.MeasureNumber)
    measnum = handles.PairedMeasuresMatrix.MeasureNumber(end)+1;
else
    measnum = 1;
end
for chidx = 1:size(y2,2);
    xdiff = x2-x1;
    ydiff = y2(chidx)-y1(chidx);
    handles.PairedMeasuresMatrix.MeasureNumber(end+1) = measnum;
    handles.PairedMeasuresMatrix.FileName{end+1} = handles.filename;
    cnn = handles.includedchannels(channelstomeasure(chidx));
    chname = handles.obj.Channel(cnn).ChannelName;
    handles.PairedMeasuresMatrix.ChannelName{end+1} = chname;
    handles.PairedMeasuresMatrix.x1s(end+1) = x1;
    handles.PairedMeasuresMatrix.y1s(end+1) = y1(chidx);
    handles.PairedMeasuresMatrix.x2s(end+1) = x2;
    handles.PairedMeasuresMatrix.y2s(end+1) = y2(chidx);
    handles.PairedMeasuresMatrix.xdiffs(end+1) = xdiff;
    handles.PairedMeasuresMatrix.ydiffs(end+1) = ydiff;
    handles.PairedMeasuresMatrix.MeasureMode{end+1} = measuremode;
    handles.PairedMeasuresMatrix.Comments{end+1} = '';
end

handles = PairedMeasuresBox(handles);

guidata(handles.displayfig,handles);

function handles = PairedMeasuresBox(handles)
%% Creates/handles auxiliary figure to display the Paired Measures Data
thisfig = findobj('type','figure','name','Paired Measures','userdata',handles.displayfig);
set(0,'units','points');
screenpix = get(0,'ScreenSize');
pixperinch = get(0,'ScreenPixelsPerInch');
inches = screenpix/pixperinch;
screenpts = inches * 72;
boxwidth = 525;

if isempty(thisfig);
    handles.PairedMeasuresWindow = figure('units','points',...
        'position',[screenpts(3)-boxwidth 40 boxwidth 1],...
        'toolbar','none',...
        'menubar','none',...
        'name','Paired Measures',...
        'numbertitle','off',...
        'userdata',handles.displayfig);
else
    handles.PairedMeasuresWindow = thisfig;
end

delete(get(handles.PairedMeasuresWindow,'children'));

precision = 8;
lineheight = 12;
fontsize = 8;
set(handles.PairedMeasuresWindow,'units','points');
nummeas = length(handles.PairedMeasuresMatrix.x1s);
winsz = get(handles.PairedMeasuresWindow,'position');
winsz = [winsz(1) max([0 winsz(2)-lineheight]) winsz(3) (nummeas+2)*lineheight];
set(handles.PairedMeasuresWindow,'position',winsz);

handles.PairedMeasuresExportMenu = uimenu('Label','Export','Parent',handles.PairedMeasuresWindow);
uimenu(handles.PairedMeasuresExportMenu,'Label','To Workspace','callback',@ExportPairedMeasuresToWS);
uimenu(handles.PairedMeasuresExportMenu,'Label','To .mat File','callback',@ExportPairedMeasuresToMAT);
uimenu(handles.PairedMeasuresExportMenu,'Label','To .xls File','callback',@ExportPairedMeasuresToXLS);
uimenu(handles.PairedMeasuresExportMenu,'Label','To current Excel Doc','callback',@ExportPairedMeasuresToCurrentExcelBook);
handles.PairedMeasuresDeleteMenu = uimenu('Label','Delete','Parent',handles.PairedMeasuresWindow);
uimenu(handles.PairedMeasuresDeleteMenu,'Label','All','callback',@DeletePairedMeasuresAll);
uimenu(handles.PairedMeasuresDeleteMenu,'Label','Select','callback',@DeletePairedMeasuresSelect);



%Column titles
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[10 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','#');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[35 winsz(4)-12 55 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Channel');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[120 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','X1');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[167 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Y1');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[220 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','X2');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[267 winsz(4)-12 12 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Y2');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[311 winsz(4)-12 30 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','X2-X1');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[358 winsz(4)-12 30 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Y2-Y1');
uicontrol('parent',handles.PairedMeasuresWindow,'style','Text','units','points',...
    'position',[430 winsz(4)-12 55 12],...
    'BackgroundColor',[.8 .8 .8],'fontweight','bold','fontsize',fontsize,'String','Comment');

nummeas = length(handles.PairedMeasuresMatrix.x1s);
for midx = 1:nummeas;
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[7  winsz(4)-((lineheight-1)*(midx+1))-1 20 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[8  winsz(4)-((lineheight-1)*(midx+1)) 18 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.MeasureNumber(midx)));
    
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[26  winsz(4)-((lineheight-1)*(midx+1))-1 75 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[27  winsz(4)-((lineheight-1)*(midx+1)) 73 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',handles.PairedMeasuresMatrix.ChannelName{midx});
    
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[101  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[102  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.x1s(midx),precision));
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[150  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[151  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.y1s(midx),precision));
    
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[201  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[202  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.x2s(midx),precision));
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[250  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[251  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.y2s(midx),precision));
    
    
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[301  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[302  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.xdiffs(midx),precision));
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[350  winsz(4)-((lineheight-1)*(midx+1))-1 50 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','text','units','points',...
        'position',[351  winsz(4)-((lineheight-1)*(midx+1)) 48 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',num2str(handles.PairedMeasuresMatrix.ydiffs(midx),precision));
    
    
    uicontrol('parent',handles.PairedMeasuresWindow,'style','frame','units','points',...
        'position',[401  winsz(4)-((lineheight-1)*(midx+1))-1 119 lineheight],'BackgroundColor',[1 1 1]);
    uicontrol('parent',handles.PairedMeasuresWindow,'style','edit','units','points',...
        'position',[403  winsz(4)-((lineheight-1)*(midx+1)) 117 lineheight-2],'BackgroundColor',[1 1 1],...
        'fontsize',fontsize,'String',handles.PairedMeasuresMatrix.Comments{midx},...
        'tag',num2str(midx),'userdata',handles.displayfig,'callback',@PairedMeasuresCommentFcn);
end

function SingleMeasureModeFcn(obj,ev)
handles = guidata(obj);
choicelist = {'Click';'Mean';'Maximum';'Minimum'};
initialidx = strmatch(handles.SingleMeasureMode, choicelist,'exact');
[Selection,ok] = listdlg('ListString',choicelist,'SelectionMode','Single',...
    'InitialValue',initialidx,'PromptString','Choose a measurement mode');
if ~ok
    return
end
handles.SingleMeasureMode = choicelist{Selection};
set(handles.SingleMeasureModeIndic,'String',choicelist{Selection});
guidata(handles.displayfig,handles);

function PairedMeasureModeFcn(obj,ev)
handles = guidata(obj);
choicelist = {'Click & Click';'Click & Mean';'Click & Maximum';'Click & Minimum';...
    'Mean & Click';'Mean & Mean';'Mean & Maximum';'Mean & Minimum';...
    'Maximum & Click';'Maximum & Mean';'Maximum & Maximum';'Maximum & Minimum';...
    'Minimum & Click';'Minimum & Mean';'Minimum & Maximum';'Minimum & Minimum'};
initialidx = strmatch(handles.PairedMeasureMode, choicelist,'exact');
[Selection,ok] = listdlg('ListString',choicelist,'SelectionMode','Single',...
    'InitialValue',initialidx,'PromptString','Choose a measurement mode');
if ~ok
    return
end
handles.PairedMeasureMode = choicelist{Selection};
set(handles.PairedMeasureModeIndic,'String',choicelist{Selection});
guidata(handles.displayfig,handles);

function DeleteSingleMeasuresAll(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));

handles.SingleMeasuresMatrix.MeasureNumber = [];
handles.SingleMeasuresMatrix.FileName = {};
handles.SingleMeasuresMatrix.ChannelName = {};
handles.SingleMeasuresMatrix.xs = [];
handles.SingleMeasuresMatrix.ys = [];
handles.SingleMeasuresMatrix.MeasureMode = {};
handles.SingleMeasuresMatrix.Comments = {};

delete(handles.SingleMeasuresExportMenu);
handles.SingleMeasuresExportMenu = [];
delete(handles.SingleMeasuresDeleteMenu);
handles.SingleMeasuresDeleteMenu = [];

set(thisfig,'userdata',[]);
set(thisfig,'name','Single Measures Deleted')

guidata(handles.displayfig,handles);

function DeleteSingleMeasuresSelect(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));

choices = {};
for midx = 1:size(handles.SingleMeasuresMatrix.MeasureNumber,2)
    choices{end+1} = [num2str(handles.SingleMeasuresMatrix.MeasureNumber(midx)),' ',...
        handles.SingleMeasuresMatrix.ChannelName{midx}];
end
[Selection,ok] = listdlg('ListString',choices,'PromptString','Choose measurements to delete');
if ~ok
    return
end

if length(Selection) == length(choices);%if user chose all measures... reinitialize all matrices
    handles.SingleMeasuresMatrix.MeasureNumber = [];
    handles.SingleMeasuresMatrix.FileName = {};
    handles.SingleMeasuresMatrix.ChannelName = {};
    handles.SingleMeasuresMatrix.xs = [];
    handles.SingleMeasuresMatrix.ys = [];
    handles.SingleMeasuresMatrix.MeasureMode = {};
    handles.SingleMeasuresMatrix.Comments = {};
    
    set(thisfig,'userdata',[]);
    set(thisfig,'name','Single Measures - Deleted')
else%else remove just selected measures
    handles.SingleMeasuresMatrix.MeasureNumber(Selection) = [];
    handles.SingleMeasuresMatrix.FileName(Selection) = [];
    handles.SingleMeasuresMatrix.ChannelName(Selection) = [];
    handles.SingleMeasuresMatrix.xs(Selection) = [];
    handles.SingleMeasuresMatrix.ys(Selection) = [];
    handles.SingleMeasuresMatrix.MeasureMode(Selection) = [];
    handles.SingleMeasuresMatrix.Comments(Selection) = [];
    
    set(thisfig,'userdata',[]);
    set(thisfig,'name','Single Measures - Partially Deleted')
end

delete(handles.SingleMeasuresExportMenu);
handles.SingleMeasuresExportMenu = [];
delete(handles.SingleMeasuresDeleteMenu);
handles.SingleMeasuresDeleteMenu = [];

guidata(handles.displayfig,handles);

function DeletePairedMeasuresAll(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%PairedMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));

handles.PairedMeasuresMatrix.MeasureNumber = [];
handles.PairedMeasuresMatrix.FileName = {};
handles.PairedMeasuresMatrix.ChannelName = {};
handles.PairedMeasuresMatrix.x1s = [];
handles.PairedMeasuresMatrix.y1s = [];
handles.PairedMeasuresMatrix.x2s = [];
handles.PairedMeasuresMatrix.y2s = [];
handles.PairedMeasuresMatrix.xdiffs = [];
handles.PairedMeasuresMatrix.ydiffs = [];
handles.PairedMeasuresMatrix.MeasureMode = {};
handles.PairedMeasuresMatrix.Comments = {};

delete(handles.PairedMeasuresExportMenu);
handles.PairedMeasuresExportMenu = [];
delete(handles.PairedMeasuresDeleteMenu);
handles.PairedMeasuresDeleteMenu = [];

set(thisfig,'userdata',[]);
set(thisfig,'name','Paired Measures Deleted')

guidata(handles.displayfig,handles);

function DeletePairedMeasuresSelect(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%PairedMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));

choices = {};
for midx = 1:size(handles.PairedMeasuresMatrix.MeasureNumber,2)
    choices{end+1} = [num2str(handles.PairedMeasuresMatrix.MeasureNumber(midx)),' ',...
        handles.PairedMeasuresMatrix.ChannelName{midx}];
end
[Selection,ok] = listdlg('ListString',choices,'PromptString','Choose measurements to delete');
if ~ok
    return
end

if length(Selection) == length(choices);%if user chose all measures... reinitialize all matrices
    handles.PairedMeasuresMatrix.MeasureNumber = [];
    handles.PairedMeasuresMatrix.FileName = {};
    handles.PairedMeasuresMatrix.ChannelName = {};
    handles.PairedMeasuresMatrix.x1s = [];
    handles.PairedMeasuresMatrix.y1s = [];
    handles.PairedMeasuresMatrix.x2s = [];
    handles.PairedMeasuresMatrix.y2s = [];
    handles.PairedMeasuresMatrix.xdiffs = [];
    handles.PairedMeasuresMatrix.ydiffs = [];
    handles.PairedMeasuresMatrix.MeasureMode = {};
    handles.PairedMeasuresMatrix.Comments = {};
    
    set(thisfig,'userdata',[]);
    set(thisfig,'name','Paired Measures - Deleted')
else%else remove just selected measures
    handles.PairedMeasuresMatrix.MeasureNumber(Selection) = [];
    handles.PairedMeasuresMatrix.FileName(Selection) = [];
    handles.PairedMeasuresMatrix.ChannelName(Selection) = [];
    handles.PairedMeasuresMatrix.x1s(Selection) = [];
    handles.PairedMeasuresMatrix.y1s(Selection) = [];
    handles.PairedMeasuresMatrix.x2s(Selection) = [];
    handles.PairedMeasuresMatrix.y2s(Selection) = [];
    handles.PairedMeasuresMatrix.xdiffs(Selection) = [];
    handles.PairedMeasuresMatrix.ydiffs(Selection) = [];
    handles.PairedMeasuresMatrix.Comments(Selection) = [];
    handles.PairedMeasuresMatrix.MeasureMode(Selection) = [];
    
    set(thisfig,'userdata',[]);
    set(thisfig,'name','Paired Measures - Partially Deleted')
end

delete(handles.PairedMeasuresExportMenu);
handles.PairedMeasuresExportMenu = [];
delete(handles.PairedMeasuresDeleteMenu);
handles.PairedMeasuresDeleteMenu = [];

guidata(handles.displayfig,handles);

function ExportSingleMeasuresToWS(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = Getsinglemeasuresoutputmat(handles);

assignin('base','SingleMeasures',outputmat);

function ExportSingleMeasuresToMAT(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = Getsinglemeasuresoutputmat(handles);

[FileName,PathName] = uiputfile('.mat','Select a File to Write',['SingleMeasures_',handles.filename]);
if FileName == 0 & PathName == 0
    return
end
save([PathName,FileName],'outputmat');

function ExportSingleMeasuresToXLS(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = Getsinglemeasuresoutputmat(handles);

[FileName,PathName] = uiputfile('.xls','Select a File to Write',['SingleMeasures_',handles.filename]);
if FileName == 0 & PathName == 0
    return
end
xlswrite([PathName,FileName],outputmat);

function ExportSingleMeasuresToCurrentExcelBook(obj,ev);
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = Getsinglemeasuresoutputmat(handles);

answer = inputdlg('Enter sheet cell number to insert at');
if isempty(answer)
    return
end
answer=lower(answer{1});
letter = isstrprop(answer,'alpha');
if sum(letter) > 1
    error('Can only handle single-letter column names');
end
if letter(1) ~= 1;
    error('Please enter a single letter followed by a number');
end
numbers = isstrprop(answer,'digit');
letter = answer(letter);
numbers = answer(numbers);
colnum = letter-'a'+1;
rownum = str2double(numbers);

chan = ddeinit('excel', 'Sheet1');
for midx = 1:size(outputmat,1);
    for iidx = 1:size(outputmat,2);
        rc = ['r',num2str(rownum+midx-1),'c',num2str(colnum+iidx-1)];%start/stop cell number
        rc = [rc,':',rc];
        ddepoke(chan,rc,outputmat{midx,iidx})
    end
end
%go thru each element and ddepoke
ddeterm(chan);

function outputmat = Getsinglemeasuresoutputmat(handles);

for midx = 1:length(handles.SingleMeasuresMatrix.MeasureNumber);
    outputmat{midx,1} = handles.SingleMeasuresMatrix.FileName{midx};
    outputmat{midx,2} = handles.SingleMeasuresMatrix.MeasureNumber(midx);
    outputmat{midx,3} = handles.SingleMeasuresMatrix.ChannelName{midx};
    outputmat{midx,4} = handles.SingleMeasuresMatrix.xs(midx);
    outputmat{midx,5} = handles.SingleMeasuresMatrix.ys(midx);
    outputmat{midx,6} = handles.SingleMeasuresMatrix.MeasureMode{midx};
    outputmat{midx,7} = handles.SingleMeasuresMatrix.Comments{midx};
end

function ExportPairedMeasuresToWS(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%PairedMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = Getpairedmeasuresoutputmat(handles);

assignin('base','PairedMeasures',outputmat);

function ExportPairedMeasuresToMAT(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%PairedMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = Getpairedmeasuresoutputmat(handles);

[FileName,PathName] = uiputfile('.mat','Select a File to Write',['PairedMeasures_',handles.filename]);
if FileName == 0 & PathName == 0
    return
end
save([PathName,FileName],'outputmat');

function ExportPairedMeasuresToXLS(obj,ev)
thisfig = get(get(obj,'parent'),'parent');%PairedMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = Getpairedmeasuresoutputmat(handles);

[FileName,PathName] = uiputfile('.xls','Select a File to Write',['PairedMeasures_',handles.filename]);
if FileName == 0 & PathName == 0
    return
end
xlswrite([PathName,FileName],outputmat);

function ExportPairedMeasuresToCurrentExcelBook(obj,ev);
thisfig = get(get(obj,'parent'),'parent');%SingleMeasures Display Fig
handles = guidata(get(thisfig,'userdata'));
outputmat = Getpairedmeasuresoutputmat(handles);

answer = inputdlg('Enter sheet cell number to insert at');
if isempty(answer)
    return
end
answer=lower(answer{1});
letter = isstrprop(answer,'alpha');
if sum(letter) > 1
    error('Can only handle single-letter column names');
end
if letter(1) ~= 1;
    error('Please enter a single letter followed by a number');
end
numbers = isstrprop(answer,'digit');
letter = answer(letter);
numbers = answer(numbers);
colnum = letter-'a'+1;
rownum = str2double(numbers);

chan = ddeinit('excel','Sheet1');
for midx = 1:size(outputmat,1);
    for iidx = 1:size(outputmat,2);
        rc = ['r',num2str(rownum+midx-1),'c',num2str(colnum+iidx-1)];%start/stop cell number
        rc = [rc,':',rc];
        ddepoke(chan,rc,outputmat{midx,iidx});
    end
end
%go thru each element and ddepoke
ddeterm(chan);

function outputmat = Getpairedmeasuresoutputmat(handles);

for midx = 1:length(handles.PairedMeasuresMatrix.MeasureNumber);
    outputmat{midx,1} = handles.PairedMeasuresMatrix.FileName{midx};
    outputmat{midx,2} = handles.PairedMeasuresMatrix.MeasureNumber(midx);
    outputmat{midx,3} = handles.PairedMeasuresMatrix.ChannelName{midx};
    outputmat{midx,4} = handles.PairedMeasuresMatrix.x1s(midx);
    outputmat{midx,5} =  handles.PairedMeasuresMatrix.y1s(midx);
    outputmat{midx,6} =  handles.PairedMeasuresMatrix.x2s(midx);
    outputmat{midx,7} =  handles.PairedMeasuresMatrix.y2s(midx);
    outputmat{midx,8} =  handles.PairedMeasuresMatrix.xdiffs(midx);
    outputmat{midx,9} =  handles.PairedMeasuresMatrix.ydiffs(midx);
    outputmat{midx,10} = handles.PairedMeasuresMatrix.MeasureMode{midx};
    outputmat{midx,11} = handles.PairedMeasuresMatrix.Comments{midx};
end

function PeakTrigAvgFcn(obj,ev)

answer = inputdlg({'Trigger Channel (Axis Number)',...
    'Trigger Level (Triggers above this level)',...
    'Trigger Minimum Duration (ms)',...
    'Trigger Maximum Duration (ms)',...
    'Channel to Average (Axis Number.  If multiple channels separate by comma.)',...
    'Milliseconds Before Trigger to Grab','Milliseconds After Trigger to Grab',...
    'Milliseconds After Trigger To Ignore New Triggers',...
    'Start of data range to use (seconds)','End of data range to use (seconds)',...
    'Ignore Trigger Deviations Less Than Following # of Data Points',...
    'Failure analysis (requires template matching) 1-yes 2-extensive(on display 2)'},...
    'Peak-Triggered Average',...% dialog title
    1,...%default number of lines per box
    {'3','0.1','0','Inf','2','10','50','0','0','Inf','0','0'});%default answers.  Default range is all of data
if length(answer)~=12
    return
end
if isempty(answer{1}) || isempty(answer{2}) ||...
        isempty(answer{3}) || isempty(answer{4}) || isempty(answer{5}) ||...
        isempty(answer{6}) || isempty(answer{7}) || isempty(answer{8}) ||...
        isempty(answer{9}) || isempty(answer{10})|| isempty(answer{11})||...
        isempty(answer{12})
    return
end
handles = guidata(obj);
decdata=getappdata(handles.displayfig,'decdata');
allrawdata = decdata.data{1};
clear decdata

trigchan = str2double(answer{1});
abovethresh = str2double(answer{2});
mindur = str2double(answer{3});
maxdur = str2double(answer{4});
beforems = str2double(answer{6});
afterms = str2double(answer{7});
trigignoredelayms = str2double(answer{8});
starttime = str2double(answer{9});
stoptime = str2double(answer{10});
ignoredeviations = str2double(answer{11});
dofailure = str2double(answer{12});

%if multiple channels entered to be averaged
commas = findstr(answer{5},',');
if isempty(commas)
    avgchan = str2double(answer{5});
else%hopefully resistant to comma+space separation
    for cidx = 1:length(commas);
        avgchan(cidx) = str2double(answer{5}(commas(cidx)-1));
    end
    avgchan(end+1) =  str2double(answer{5}(end));
end

beforepts = round(beforems * handles.obj.SampleRate / 1000);
afterpts = round(afterms * handles.obj.SampleRate / 1000);
trigignoredelaypts = round(trigignoredelayms * handles.obj.SampleRate / 1000);
mindur = round(mindur * handles.obj.SampleRate / 1000);
maxdur = round(maxdur * handles.obj.SampleRate / 1000);
if starttime > stoptime
    errordlg 'Start time must be less than stop time'
    return
end
startpt = starttime * handles.obj.SampleRate;
startpt = max([startpt 1]);%startpt is at least point number 1
stoppt = stoptime * handles.obj.SampleRate;
stoppt = min([stoppt size(allrawdata,1)]);%stoppt not greater than samples acquired

%get triggers
trigdata = allrawdata(:,trigchan);
trigtimes = pt_continuousabove(trigdata,zeros(size(trigdata)),abovethresh,mindur,maxdur,ignoredeviations);
if isempty(trigtimes)
    msgbox('No triggers found')
    return
end
trigtimes = trigtimes(:,1);
trigtimes(trigtimes > stoppt) = [];
trigtimes(trigtimes < startpt) = [];
% trigtraces = zeros(length(trigtimes),(beforepts+afterpts+1));
lastsavedtrig = -Inf;
trigtraces = [];
newtrigtimes = [];
for tidx = 1:length(trigtimes)
    if trigtimes(tidx)>(lastsavedtrig+trigignoredelaypts)
        trigtraces(end+1,:) = trigdata(trigtimes(tidx)-beforepts:trigtimes(tidx)+afterpts);
        lastsavedtrig = trigtimes(tidx);
        newtrigtimes(end+1) = lastsavedtrig;
    end
end
trigtimes = newtrigtimes;

for cidx = 1:length(avgchan);
    thischan = avgchan(cidx);
    avgdata = allrawdata(:,thischan);
    rawtraces = zeros(length(trigtimes),(beforepts+afterpts+1));
    for tidx = 1:length(trigtimes)
        rawtraces(tidx,:) = avgdata(trigtimes(tidx)-beforepts:trigtimes(tidx)+afterpts);
        %         use the sgolay filter to smooth the trace. doesnt do jack
        %         shit for 60Hz noise. AP 20091008
        %                 rawtraces(tidx,:) = sgolayfilt(rawtraces(tidx,:),3,21);
    end
    zeroedtraces = rawtraces - repmat(rawtraces(:,beforepts), [1 (beforepts+afterpts+1)]);
    
    outputfig = TriggeredAverageKernelFcn(trigtraces, rawtraces, zeroedtraces, handles.obj.SampleRate,trigtimes);
    trigname = handles.obj.Channel(handles.includedchannels(trigchan)).ChannelName;
    grabbedname = handles.obj.Channel(handles.includedchannels(avgchan(cidx))).ChannelName;
    figname = [handles.filename,': "',trigname,'" triggering "',grabbedname,'"'];
    set(outputfig,'name',figname);
end

% Adam Packer post-synaptic failure analysis
% Assumes user has completed template analysis on the post channel
% Assumes the actual event of interest starts 12.9 or 14.1 msec into the template
% (This is based on EPSCtemplate and EPSPtemplate - specific templates made by me.)
% Assumes the user only cares about post-synaptic events matched within 5 ms
% Finds the minimum based on things within 10 ms
% (This is because the template isn't perfect, nor are the triggers.)

samplesperms = handles.obj.SampleRate / 1000;
if dofailure
    for cidx=1:length(avgchan);
        numtriggers=length(trigtimes);
        if isfield(handles,'TemplateResults')
            thischan = avgchan(cidx);
            handles.TemplateResults(thischan).lookingfor=1;
            if handles.TemplateResults(thischan).lookingfor==1 %epsc
                eventdelay = 10 * samplesperms;
                postdelay = 20 * samplesperms;
                amplitudes=min(zeroedtraces(:,beforepts:(beforepts+10*samplesperms)),[],2);
            elseif handles.TemplateResults(thischan).lookingfor==0 % epsp
                eventdelay = 24.0 * samplesperms;
                postdelay = 20 * samplesperms;
                amplitudes=max(zeroedtraces(:,beforepts:(beforepts+10*samplesperms)),[],2);
            elseif handles.TemplateResults(thischan).lookingfor==2 % ipsc
                eventdelay = 12.9 * samplesperms;
                postdelay = 20 * samplesperms;
                amplitudes=max(zeroedtraces(:,beforepts:(beforepts+10*samplesperms)),[],2);
            end
            templatetimes = handles.TemplateResults(thischan).columndetected + eventdelay;
            [X,Y] = meshgrid(trigtimes,templatetimes);
            latencies = Y - X;  %find latency from AP to EPSC
            latencies(latencies<(0-eventdelay))=Inf; %set all negative latencies to Inf for next step
            [EPSCidx,successes] = find(latencies < postdelay); %find any latencies less than postdelay
            %         successlatency=latencies(EPSCidx,successes);
            %         successlatency=spdiags(successlatency,0);
            successlatency=zeros(1,length(successes));
            for j=1:length(successes)
                [peak,peakidx] = max(trigtraces(j,beforepts:(beforepts+5*samplesperms)));
                longhighdiffs=pt_continuousabove(abs(zscore(diff(zeroedtraces(j,beforepts:(beforepts+5*samplesperms))))),0,0.6,5,100,0);
                if isempty(longhighdiffs) || (longhighdiffs(1)-peakidx<=0)
                    % following min only works for inward events!
                    [onset,onsetidx] = min(zscore(diff(zeroedtraces(j,beforepts:(beforepts+5*samplesperms)))));
                else
                    onsetidx=longhighdiffs(1);
                end
                thislatency = (onsetidx - peakidx) / samplesperms;
                if thislatency <= 0
                    successlatency(j) = NaN;
                else
                    successlatency(j) = thislatency;
                end
            end
            meanlatency = mean(successlatency(~isnan(successlatency)));
            successamp = amplitudes(successes);
            meanamp = mean(successamp);
            if dofailure==2
                figure('units','normalized','position',[1.0016 0.4767 0.3089 0.7317]);
                hist(successamp);title(sprintf('%s \nHistogram success amplitudes',figname))
                figure('units','normalized','position',[1.3146 0.4767 0.3083 0.7317]);
                stem(successamp);title(sprintf('%s \nSuccess amplitudes in order',figname))
                figure('units','normalized','position',[1.0016   -0.3233    0.3089    0.7317]);
                hist(successlatency);title(sprintf('%s \nHistogram success latencies',figname))
                figure('units','normalized','position',[1.3146 -0.3233 0.3089 0.7317]);
                stem(successlatency);title(sprintf('%s \nSuccess latencies in order',figname))
            end
            numsuccesses = length(successes);
            numfailures = numtriggers - numsuccesses;
            probsuccess = (length(successes)) / numtriggers * 100;
            successes = [0 successes']; % to handle if the string of failures is at the beginning
            diffsuccess = [diff(successes) (1 + numtriggers - successes(end))]; %add the diff from the end to numtriggers to the end
            [longestrunfailures] = max(diffsuccess) - 1; %diff overcounts by one
            if isempty(longestrunfailures) % handle the case when there are no failures
                longestrunfailures = 0;
            end
            probfailrun=RunCalculator(longestrunfailures,numtriggers,1-(probsuccess/100));
            diffdiffsuccess = diff(diffsuccess);
            diffdiffsuccess = [1 diffdiffsuccess 1];
            [longestrunsuccesses,longestrunsuccessesidx] = max(diff(find(logical(diffdiffsuccess))));
            if longestrunsuccessesidx==1
                longestrunsuccesses = longestrunsuccesses;
            else
                longestrunsuccesses = longestrunsuccesses + 1;
            end
            probsuccessrun=RunCalculator(longestrunsuccesses,numtriggers,probsuccess/100);
            % Setup data if have epscs
            colnames = {'Id' 'Trigs', 'Latency', 'Amplitude', 'Successes', 'Failures', 'P Success',...
                'Max Run Failures','P max run failures',...
                'Max Run Successes','P max run successes'};
            data = {figname numtriggers meanlatency meanamp numsuccesses numfailures ...
                probsuccess longestrunfailures probfailrun ...
                longestrunsuccesses probsuccessrun};
        end
        
        % Plot all the stimulations
        f = figure('units','normalized','position',[0.3 0.05 0.58 0.6]);
        if isfield(handles,'TemplateResults')
            t = uitable(f, 'Units','normalized', ...
                'Data', data, ...
                'Position', [0 0.95 1 0.05], ...
                'RowName',[], ...
                'ColumnName', colnames, ...
                'ColumnWidth',{270 'auto' 'auto' 'auto' 'auto' 'auto' 'auto' 'auto' 'auto' 'auto' 'auto'});
        end
        
        % Deal with having e.g. 501 data points (just scale to 500)
        if mod(size(rawtraces,2),100)==1
            STAxlim = [0 100 * floor(size(rawtraces,2)/100)];
        else
            STAxlim = 'auto';
        end
        
        % Subplot each trace
        m = floor(sqrt(numtriggers));
        n = ceil(numtriggers/m);
        x=0:length(zeroedtraces(1,:))-1;
        for i=1:numtriggers
            subplot(m,n,i);
            [AX,h1,h2]=plotyy(x,zeroedtraces(i,:),x,trigtraces(i,:));
            xlim(AX(1),STAxlim);
            xlim(AX(2),STAxlim);
            if isfield(handles,'TemplateResults')
                if find(successes==i);
                    title('success')
                else
                    title('failure')
                end
            end
        end
        
        if isfield(handles,'TemplateResults')
            goodzeroedtraces=mean(zeroedtraces(successes(2:end),:));
            goodtrigtraces=mean(trigtraces(successes(2:end),:));
            x=0:length(goodzeroedtraces)-1;
            figure;plotyy(x,goodzeroedtraces,x,goodtrigtraces);xlim(STAxlim);
            %             figure;subplot(2,1,1);plot(goodzeroedtraces);subplot(2,1,2);plot(goodtrigtraces);
        end
    end
end

function DipTrigAvgFcn(obj,ev);

answer = inputdlg({'Trigger Channel (Axis Number)',...
    'Trigger Level (Triggers below this level)',...
    'Trigger Minimum Duration (ms)',...
    'Trigger Maximum Duration (ms)',...
    'Channel to Average (Axis Number.  If multiple channels separate by comma.)',...
    'Milliseconds Before Trigger to Grab','Milliseconds After Trigger to Grab',...
    'Milliseconds After Trigger To Ignore New Triggers',...
    'Start of data range to use (seconds)','End of data range to use (seconds)',...
    'Ignore Deviations Less Than Following # of Data Points'},...
    'Dip-Triggered Average',...% dialog title
    1,...%default number of lines per box
    {'1','1','0','Inf','2','10','20','0','0','Inf','0'});%default answers.  Default range is all of data
if length(answer)~=11
    return
end
if isempty(answer{1}) || isempty(answer{2}) ||...
        isempty(answer{3}) || isempty(answer{4}) || isempty(answer{5}) ||...
        isempty(answer{6}) || isempty(answer{7}) || isempty(answer{8}) ||...
        isempty(answer{9}) || isempty(answer{10}) || isempty(answer{11})
    return
end
handles = guidata(obj);
decdata=getappdata(handles.displayfig,'decdata');
allrawdata = decdata.data{1};
clear decdata

trigchan = str2double(answer{1});
belowthresh = str2double(answer{2});
mindur = str2double(answer{3});
maxdur = str2double(answer{4});
beforems = str2double(answer{6});
afterms = str2double(answer{7});
trigignoredelayms = str2double(answer{8});
starttime = str2double(answer{9});
stoptime = str2double(answer{10});
ignoredeviations = str2double(answer{11});

%if multiple channels entered to be averaged
commas = findstr(answer{5},',');
if isempty(commas)
    avgchan = str2double(answer{5});
else
    for cidx = 1:length(commas);
        avgchan(cidx) = str2double(answer{5}(commas(cidx)-1));
    end
    avgchan(end+1) =  str2double(answer{5}(end));
end

beforepts = round(beforems * handles.obj.SampleRate / 1000);
afterpts = round(afterms * handles.obj.SampleRate / 1000);
trigignoredelaypts = round(trigignoredelayms * handles.obj.SampleRate / 1000);
mindur = round(mindur * handles.obj.SampleRate / 1000);
maxdur = round(maxdur * handles.obj.SampleRate / 1000);
if starttime > stoptime
    errordlg 'Start time must be less than stop time'
    return
end
startpt = starttime * handles.obj.SampleRate;
startpt = max([startpt 1]);%startpt is at least point number 1
stoppt = stoptime * handles.obj.SampleRate;
stoppt = min([stoppt size(allrawdata,1)]);%stoppt not greater than samples acquired

%get triggers
trigdata = allrawdata(:,trigchan);
trigtimes = pt_continuousbelow(trigdata,zeros(size(trigdata)),belowthresh,mindur,maxdur,ignoredeviations);
trigtimes = trigtimes(:,1);
trigtimes(trigtimes > stoppt) = [];
trigtimes(trigtimes < startpt) = [];
lastsavedtrig = -Inf;
trigtraces = [];
newtrigtimes = [];
for tidx = 1:length(trigtimes)
    if trigtimes(tidx)>(lastsavedtrig+trigignoredelaypts)
        trigtraces(end+1,:) = trigdata(trigtimes(tidx)-beforepts:trigtimes(tidx)+afterpts);
        lastsavedtrig = trigtimes(tidx);
        newtrigtimes(end+1) = lastsavedtrig;
    end
end
trigtimes = newtrigtimes;

for cidx = 1:length(avgchan);
    thischan = avgchan(cidx);
    avgdata = allrawdata(:,thischan);
    rawtraces = zeros(length(trigtimes),(beforepts+afterpts+1));
    for tidx = 1:length(trigtimes)
        rawtraces(tidx,:) = avgdata(trigtimes(tidx)-beforepts:trigtimes(tidx)+afterpts);
    end
    zeroedtraces = rawtraces - repmat(rawtraces(:,beforepts), [1 (beforepts+afterpts+1)]);
    
    outputfig = TriggeredAverageKernelFcn(trigtraces, rawtraces, zeroedtraces, handles.obj.SampleRate,trigtimes);
    trigname = handles.obj.Channel(handles.includedchannels(trigchan)).ChannelName;
    grabbedname = handles.obj.Channel(handles.includedchannels(avgchan(cidx))).ChannelName;
    figname = [handles.filename,': "',trigname,'" triggering "',grabbedname,'"'];
    set(outputfig,'name',figname);

end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TemplateMatchingFcn(obj,evt)
handles=guidata(obj);
decdata=getappdata(handles.displayfig,'decdata');
allrawdata = decdata.data{1};
clear decdata
params=APinputdlg({'Channel','Looking for Cell_attached_Spikes=3 IPSCs=2 EPSCs=1 EPSPs=0','Coefficient of correlation threshold',...
    'Start of data range to use (seconds)',...
    'End of data range to use (seconds)',...
    'Flip the template? 1=yes'},'Parameters',1,{'1' '1' '0.6' '0' 'Inf' '0'});
if length(params)~=6
    return
end
if isempty(params{1}) || isempty(params{2}) ||...
        isempty(params{3}) || isempty(params{4}) ||...
        isempty(params{5}) || isempty(params{6})
    return
end
channel=str2double(params{1});
lookingfor=str2double(params{2});
threshold=str2double(params{3});
starttime = str2double(params{4});
stoptime = str2double(params{5});
flipit = str2double(params{6});
if lookingfor==1
    template=load('EPSCtemplate.mat');
    sense = 'inward';
elseif lookingfor==0
    template=load('EPSPtemplate.mat');
    sense = 'outward';
elseif lookingfor==2
    template=load('IPSCtemplateUP.mat');
    sense = 'outward';
elseif lookingfor==3
    template=load('CellAttachSpiketemplate.mat');
    sense = 'inward';
end
if flipit
    invert = 'yes';
elseif ~flipit
    invert = 'no';
end
startpt = starttime * handles.obj.SampleRate;
startpt = max([startpt 1]);%startpt is at least point number 1
stoppt = stoptime * handles.obj.SampleRate;
stoppt = min([stoppt size(allrawdata,1)]);%stoppt not greater than samples acquired
rateabs = handles.obj.SampleRate / 1000; % template matcher takes acq rate in kHz
templateresults = APTemplateMatching(allrawdata(startpt:stoppt,channel),invert,sense,template.template,threshold,rateabs);
if lookingfor==3
    handles.TemplateResults(channel).columndetected = templateresults.columndetected + startpt - 1 + 150;
else
    handles.TemplateResults(channel).columndetected = templateresults.columndetected + startpt - 1;
end
handles.TemplateResults(channel).lookingfor = lookingfor;
% *NB! Previous line might cause events to be one data point off if
% template matcher starts at 0!
guidata(handles.displayfig,handles);%pass data back to the figure

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TemplateMatchingSaveFcn(obj,evt)
handles=guidata(obj);
[filename,pathname,FilterIndex]=uiputfile('*.mat','Choose a data file to write EPSCs');
if ~FilterIndex
    return
end
OuterStruct=handles.TemplateResults;
save([pathname,filename],'OuterStruct');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function TemplateMatchingLoadFcn(obj,evt)
handles=guidata(obj);
[filename,pathname,FilterIndex]=uigetfile('*.mat','Choose a data file with EPSCs');
if ~FilterIndex
    return
end
OuterStruct=load([pathname,filename],'OuterStruct');
handles.TemplateResults=OuterStruct.OuterStruct;
guidata(handles.displayfig,handles);%pass data back to the figure

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ViewerResizing(obj,ev)
%
% handles=guidata(obj);
% for a=1:size(handles.ChannelLabels,1);
%     set(handles.axeshandles(a),'units','points');
%     pts=get(handles.axeshandles(a),'position');
%     set(handles.axeshandles(a),'units','normalized');
%     pts=pts(4);%get full height of the axes
%     for b=1:size(handles.ChannelLabels,2);
%         ypos=pts-((b-1)*18);
%         set(handles.ChannelLabels(a,b),'units','points','position',[1 ypos]);
%     end
% end


% function HorizScaleTextBoxFcn(obj,ev);
%
% handles=guidata(obj);%get general info from the figure
%
% value=str2double(get(handles.HorizScaleTextBox,'string'));%in seconds
% if ~isempty(value);
%     if value<5/handles.obj.SampleRate;
%         value=5/handles.obj.SampleRate;
%         disp('Entered value is too low')
%     end
% 	xlimits=get(handles.axeshandles(1),'XLim');
%     value=value;%set this value for later
%     meanpos=mean(xlimits);
%     newlims=[meanpos-.5*value meanpos+.5*value];
%
% 	defaultstart=1/handles.obj.SampleRate;
% 	defaultend=handles.obj.SamplesAcquired/handles.obj.SampleRate;
%     if value>defaultend%both the left and right ends will be over the limits
%         newlims=[defaultstart defaultend];%fix that
%         value=defaultend;%save value correctly so later handles.displaylength will be correct
%         disp('Entered value is too high')
%     elseif newlims(1)<defaultstart;%if only the left side is too low (but right is fine)
%         newlims=[defaultstart value];%plot from zero to the required length
%     elseif newlims(2)>defaultend;%if only the right side is too high
%         newlims=[defaultend-value defaultend];%plot required length back from the end of the data
%     end
% 	set(handles.axeshandles,'XLim',newlims);
% 	PlotMinWithinXlims(handles,newlims)
%
%     handles.displaylength=value;
% 	set(handles.HorizScaleTextBox,'String',num2str(handles.displaylength));%set value of text box in seconds, and set it late
%         %in the function so it represents the most recent information possible
% 	guidata(handles.displayfig,handles);%pass data back to the figure
% end
%