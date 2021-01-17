function h=PlotEach(x,y,opt)
% Automatic subplotting with a title
% input: x = number of things to plot
% input: y = variable to plot
% input: opt = if 0, regular plot; if 1, plotyy
% output: h = figure handle

h=figure;
z=length(y.TriggeredTraces(1,:));
ax=zeros(1,x);
for i=1:x;
    ax(i)=subplot(x,1,i);
    if opt==1
        plotyy(1:z,y.Triggers(i,:),1:z,y.TriggeredTraces(i,:));
    else
        plot(y.TriggeredTraces(i,:));
    end
end
if opt==0
    handles=get(ax(1));
    handlekid=get(handles.Children);
    size(handlekid.XData);
    if mod(size(handlekid.XData,2),100)==1
        STAxlim = [0 100 * floor(size(handlekid.XData,2)/100)];
    else
        STAxlim = 'auto';
    end
    subplot(x,1,1)
    linkaxes(ax,'x')
    xlim(STAxlim)
end