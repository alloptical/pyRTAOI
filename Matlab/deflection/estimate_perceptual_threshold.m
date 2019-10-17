function [piezoLeft, piezoRight] = estimate_perceptual_threshold(behaviour)
figure
    rightwardchoices = behaviour ;%behaviour.fullsessionsummary.responsesRight;
    VoltageMax = 2.5;

   x_ax = linspace(0,2.5,length(rightwardchoices));
    
    [params, ~,f] = sigm_fit(1:5,rightwardchoices,[],[],[],'-r',2); hold on
    plot(1:5,rightwardchoices,'ok','MarkerSize',10,'MarkerFaceColor',[.5 .5 .5])
    
    ylabel('Probability Rightward Choice')
    xlabel('Left Piezo (V)')
    set(gcf,'color','w')
    box off
    ylim([0,1])
    axis square
    set(gca,'XTick',[1,3,5])
    set(gca,'XTickLabel',{'0','1.25','2.5'});

    xArray = 1:.0001:5;
    trueVoltage = linspace(0,VoltageMax,length(xArray));
    yVals = f(params,xArray);
    x50_min = find(abs(yVals - 0.5) == min(abs(yVals - 0.5)));
    x50_amplitude = trueVoltage(x50_min);

    t1 = text(x50_amplitude+ .1 , .4,['Piezo Left = ' num2str(x50_amplitude) ' (V)']);
    t2 = text(x50_amplitude+ .1 , .35,['Piezo Right = ' num2str(VoltageMax - x50_amplitude) ' (V)']);
    t3 = text(x50_amplitude+ .1 , .45,'Estimated Perceptual threshold:');


    title(['Piezo L = ' num2str(x50_amplitude) ', R = ' num2str(VoltageMax - x50_amplitude) '  (V)'])

    t1.FontSize = 14;
    t2.FontSize = 14;
    t3.FontSize = 14;

    piezoLeft = x50_amplitude;
    piezoRight = VoltageMax - x50_amplitude;
ylim([0,1.1])
    plot([0,5],[0.5,0.5],'--k');
    hold on
    plot([3,3],[0,1],'--k');
    hold on
    plot([xArray(x50_min),xArray(x50_min)],[0,1],'-r');
  
end