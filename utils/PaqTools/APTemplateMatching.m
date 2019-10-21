function [TemplateResults]=APTemplateMatching(varargin)

%Function to select events using template matching (coeficient of correlation method)
%
%Input:
%     data (a vector containing the trace to scan)
%     invert (text to indicate template inversion or not: 'yes' or 'no') {default: 'no'}
%     sense (direction of the events to be detected: 'inward' or 'outward') {default: 'outward'}
%     template (with which to scan trace)
%
%Otputs:
%     minimat (events detected)
%     minimatbase (events detected with baseline subtracted)
%     miniav (average of events detected)
%     miniamp (peak amplitudes)
%     frequency
%     columndetected (position of each event in the original trace)
%     template
%     minimati (events detected with inverted template)
%     SNR (signal to noise ratio)
%
%Syntax:
%     [minimat, minimatbase, miniav, miniamp, frequency, columndetected, template, minimati, SNR]=templatematching(data, invert, sense);
%
%Go into the function to change parameters
%Parameter list {defaults}:
%     Acquisition rate in KHz {10}
%     Time window before event peak (in ms) {15}
%     Time window after event peak (in ms) {30}
%     Baseline starting from the first point (in ms) {8}
%     Size of chunks for automatic event detector (in ms) {5000}
%     Time window not to consider after an event is detected (in ms) {20}
%     Step size to scan the trace in ms {0.5}
%     Coefficient of correlation threshold {0.9}
%     Savitzky-Golay filter parameter number 1 (sgn1) to smooth the components of the template {3}
%     Savitzky-Golay filter parameter number 2 (sgn2) to smooth the components of the template {21}
%
%Notes:
%     Input parameters invert and sense are optional.
%     Savitzky-Golay filter parameter number 2 (sgn2) must be odd.
%     Savitzky-Golay filter parameter number 1 (sgn1) must be smaller than sgn2.
%     Use sgn1=2 and sgn2=3 to bypass Savitzky-Golay filtering.
%
%Created 06/19/2007
%Modified 09/17/2007
%inward or outward event detection option added
%
%Matlab R2006a
%by Emiliano Rial Verde

% Matlab R2007b
% Modified by Adam Packer
% 09/20/2007

%Parameters
rateabs=10; %Acquisition rate in KHz
windowbefore=15; %Time window before mini peak (in ms)
windowafter=30; %Time window after mini peak (in ms)
windowbase=8; %Baseline starting from the first point (in ms)
timesplit=5000; %Size of chunks for automatic event detector (in ms)
bufferzone=20; %Time window not to consider after a mini is detected (in ms)
stepsize=0.5; %Step size to scan the trace in ms
threshold=0.8; %Coefficient of correlation threshold
sgn1=3; %Savitzky-Golay filter parameter number 1
sgn2=21; %Savitzky-Golay filter parameter number 2


if nargin==0;
    return
elseif nargin==1;
    data=varargin{1};
    invert='no';
    sense='outward';
elseif nargin==2
    data=varargin{1};
    invert=varargin{2};
    sense='outward';
elseif nargin==3
    data=varargin{1};
    invert=varargin{2};
    sense=varargin{3};
elseif nargin==4   
    data=varargin{1};
    invert=varargin{2};
    sense=varargin{3};
    template=varargin{4};
elseif nargin==5
    data=varargin{1};
    invert=varargin{2};
    sense=varargin{3};
    template=varargin{4};
    threshold=varargin{5};
elseif nargin==6
    data=varargin{1};
    invert=varargin{2};
    sense=varargin{3};
    template=varargin{4};
    threshold=varargin{5};
    rateabs=varargin{6};
elseif nargin>6;
    return
end

if nargin<4
    %Splits the trace
    fileone=reshape(data(1:floor(length(data)/timesplit)*timesplit), timesplit, floor(length(data)/timesplit));

    %Selects the maximum/minimum for each trace and shows it to select events for template
    if strcmpi(sense, 'inward')
        [a,index]=min(fileone);
    else
        [a,index]=max(fileone);
    end
    minione=[];
    for i=1:length(index)
        if index(i)-windowbefore*rateabs<=0
        elseif index(i)+windowafter*rateabs>=size(fileone,1)
        else
            minione=[minione fileone(index(i)-windowbefore*rateabs:index(i)+windowafter*rateabs, i)];
        end
    end
    minione=minione';
    minione=flipud(sortrows(minione,windowbefore*rateabs+1));
    minione=minione';

    h0=figure;
    set(h0, 'numbertitle', 'off', ...
        'name', 'S: select. N: next. P: previous. Click: finish.', ...
        'units', 'normalized', ...
        'position', [.7 .1 .26 .8]);
    minione2=[];
    selected=0;
    i=1;
    while i<=size(minione,2)
        plot(rateabs/100:rateabs/100:windowbefore+windowafter+rateabs/100, minione(:,i));
        hold on
        ylim([min(minione(:,i))-abs(min(minione(:,i))*.02) max(minione(:,i))+abs(max(minione(:,i))*.05)]);
        xlim([rateabs/100 windowbefore+windowafter+rateabs/100]);
        set(gca, 'fontsize', 8, 'XLabel', text('string', 'msec'), 'YLabel', text('string', 'mV or pA'));
        line([windowbefore-2 windowbefore-2], [min(minione(:,i))-abs(min(minione(:,i))*.02) max(minione(:,i))+abs(max(minione(:,i))*.05)], 'Color', 'k');
        line([windowbefore+2 windowbefore+2], [min(minione(:,i))-abs(min(minione(:,i))*.02) max(minione(:,i))+abs(max(minione(:,i))*.05)], 'Color', 'k');
        hold off
        title(['S: select. N: next. P: previous. Click: finish. Selected: ', num2str(selected)]);
        w=waitforbuttonpress;
        if w==0
            break
        elseif w==1
            a=get(h0, 'currentcharacter');
            if a=='s'
                minione2=[minione2 minione(:,i)];
                selected=selected+1;
                i=i+1;
            elseif a=='p'
                if i>1
                    i=i-1;
                else
                    i=1;
                end
            else
                i=i+1;
            end
        end
    end
    close(h0);
    drawnow;

    %Substract the baseline to all the traces in minione2 and filters them before
    base=median(minione2(1:windowbase*rateabs,:),1);
    base=repmat(base, size(minione2,1), 1);
    % commented out by AP 20101210
    %     minione3=sgolayfilt(minione2, sgn1, sgn2)-base;
    minione3=minione2;

    %Scales the minis in minione3 and averages them to create the template
    peak=max(minione3);
    peak=repmat(peak, size(minione3,1), 1);
    template=mean(minione3./peak,2);
end

%Scans the trace
stepsize=stepsize*rateabs;
minimat=[];
columndetected=[];
minitemplatesize=length(template);
j=1;
% h0=figure;
% set(h0, 'numbertitle', 'off', ...
%     'menubar', 'none', ...
%     'name', 'ERV template matching', ...
%     'units', 'normalized', ...
%     'position', [.4 .5 .19 .06]);
% a=get(h0, 'Color');
% uicontrol(h0, ...
%     'Style', 'text', ...
%     'String', 'Scanning trace, please wait ...', ...
%     'HorizontalAlignment', 'center', ...
%     'BackgroundColor', a, ...
%     'Units','normalized',...
%     'Position',[0.1 0.6 0.8 0.2]);
% h1=uicontrol(h0, ...
%     'Style', 'text', ...
%     'String', 'Events found: 0', ...
%     'HorizontalAlignment', 'center', ...
%     'BackgroundColor', a, ...
%     'Units','normalized',...
%     'Position',[0.1 0.2 0.8 0.2]);
drawnow;
k=0;
tic
h = waitbar(0,'Please wait...');
while j<length(data)-minitemplatesize
    b=data(j:j+minitemplatesize-1);
    c=corrcoef(template, b);
    if c(1,2)>threshold %Coefficient of correlation threshold
        minimat=[minimat b];
        columndetected=[columndetected; j];
        j=j+bufferzone*rateabs;
        k=k+1;
        %         set(h1, 'String', ['Events found: ', num2str(k)]);
        %         drawnow;
        waitbar(j/length(data),h,['Events found: ', num2str(k)]);
    else
        j=j+stepsize;
    end
end
ElapsedTimeScanning=toc;
close(h)
% close(h0)

%To calculate detection threshold use this lines
%This inverts the template to detect noise in the shape of events
if strcmpi(invert, 'yes')
    templatei=template.*(-1);
    %Scans the trace
    minimati=[];
    j=1;
    h0=figure;
    set(h0, 'numbertitle', 'off', ...
        'menubar', 'none', ...
        'name', 'ERV template matching', ...
        'units', 'normalized', ...
        'position', [.4 .5 .19 .06]);
    a=get(h0, 'Color');
    uicontrol(h0, ...
        'Style', 'text', ...
        'String', 'Scanning inverted trace, please wait ...', ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', a, ...
        'Units','normalized',...
        'Position',[0.1 0.6 0.8 0.2]);
    h1=uicontrol(h0, ...
        'Style', 'text', ...
        'String', 'Events found: 0', ...
        'HorizontalAlignment', 'center', ...
        'BackgroundColor', a, ...
        'Units','normalized',...
        'Position',[0.1 0.2 0.8 0.2]);
    drawnow;
    k=0;
    tic
    while j<length(data)-minitemplatesize
        b=data(j:j+minitemplatesize-1);
        c=corrcoef(templatei, b);
        if c(1,2)>threshold %Coefficient of correlation threshold
            minimati=[minimati b];
            j=j+bufferzone*rateabs;
            k=k+1;
            set(h1, 'String', ['Events found: ', num2str(k)]);
            drawnow;
        else
            j=j+stepsize;
        end
    end
    toc
    close(h0)
    SNR=size(minimat,2)/size(minimati,2);
else
    minimati=[];
    SNR='NA';
end

%Calculates event frequency in events/sec
seconds=(length(data)/(rateabs*1000));
frequency=size(minimat,2)/seconds;

TemplateResults.minimat=minimat;
if k % don't calculate this stuff if you didn't find any events!
    %Calculates event amplitude
    
    minimatbase=mean(minimat(1:10*rateabs,:),1);
    minimatbase=repmat(minimatbase, size(minimat,1), 1);
    minimatbase=minimat-minimatbase;
    if strcmpi(sense, 'inward')
        miniamp=min(minimatbase);
    else
        miniamp=max(minimatbase);
    end
    
    %Calculates the average event and plots it together with the template
    miniav=mean(minimatbase,2);
    % figure;
    % plot(rateabs/100:rateabs/100:windowbefore+windowafter+rateabs/100, miniav);
    % xlim([0 windowbefore+windowafter+rateabs/100]);
    % hold on
    % if strcmpi(sense, 'inward')
    %     plot(rateabs/100:rateabs/100:windowbefore+windowafter+rateabs/100, template.*abs(min(miniav)), 'r');
    % else
    %     plot(rateabs/100:rateabs/100:windowbefore+windowafter+rateabs/100, template.*max(miniav), 'r');
    % end
    % hold off
    TemplateResults.minimatbase=minimatbase;
    TemplateResults.miniav=miniav;
    TemplateResults.miniamp=miniamp;
else
    TemplateResults.minimatbase=[];
    TemplateResults.miniav=[];
    TemplateResults.miniamp=[];
end
% set the rest of the output
TemplateResults.frequency=frequency;
TemplateResults.columndetected=columndetected;
TemplateResults.template=template;
TemplateResults.minimati=minimati;
TemplateResults.SNR=SNR;
TemplateResults.ElapsedTimeScanning=ElapsedTimeScanning;


