function [ P,x_positions,y_pos,stats] = scatter_cmp_conditions(values,plot_name,subplot,input_color_lut,varargin)
% scatter plot of fields in values
% inputs: values,plot_name,subplot,input_color_lut,varargin
% value.fields are row vectors
% TO DO: compare more than 2 fields
flag_color = 0;

if(nargin>2)
    ifsubplot = subplot;
    if(~ifsubplot)
        figure('Name',plot_name)
    end
    try
        temp = input_color_lut;
        if(~isempty(temp))
            flag_color = 1;
            color_lut = input_color_lut;
        else
            color_lut = [0.3 0.3 0.3;[231 73 24]./255;0.3 0.3 0.3];
        end
    catch
        color_lut = [0.3 0.3 0.3;[231 73 24]./255;0.3 0.3 0.3];
    end
else
    figure('Name',plot_name)
end

par = inputParser;
defaultIfPrePostOnly = 0;
default_test_type = 'ranksum';
addParameter(par,'connect_scatter',1);
addParameter(par,'PlotPrePostOnly',defaultIfPrePostOnly);
addParameter(par,'test_type',default_test_type);
addParameter(par,'IfAvgByAnimal',0);
addParameter(par,'IfBoxWhisker',0);
addParameter(par,'x_shift',0);
addParameter(par,'x_interval',1);
addParameter(par,'barwidth',.5);
addParameter(par,'barEdgecolor',[0 0 0]);
addParameter(par,'add_jitter',0);
addParameter(par,'animal_colors',[]);
addParameter(par,'IfColorAnimals',0);
addParameter(par,'tail','both');
addParameter(par,'IfHistogram',0);
addParameter(par,'plot_stats',0);
addParameter(par,'BinWidth',10);
addParameter(par,'displayopt','off');
addParameter(par,'Normalization','probability')
addParameter(par,'DisplayStyle','stairs')


parse(par,varargin{:})
if(~ isempty(par.Results.tail))
    tail = par.Results.tail;
end
if(~isempty(par.Results.PlotPrePostOnly))
    IfPlotPrePostOnly = par.Results.PlotPrePostOnly;
end
if(~isempty(par.Results.test_type))
    test_type = par.Results.test_type;
end
if(~isempty(par.Results.connect_scatter))
    connect_scatter = par.Results.connect_scatter;
end
if(~isempty(par.Results.IfBoxWhisker))
    IfBoxWhisker = par.Results.IfBoxWhisker;
end
if(~isempty(par.Results.x_shift))
    x_shift = par.Results.x_shift;
end
if(~isempty(par.Results.x_interval))
    x_interval = par.Results.x_interval;
end
if(~isempty(par.Results.barwidth))
    barwidth = par.Results.barwidth;
end
if(~isempty(par.Results.barEdgecolor))
    barEdgecolor = par.Results.barEdgecolor;
end
if(~isempty(par.Results.add_jitter))
    IfAddJitter = par.Results.add_jitter;
end
if(~isempty(par.Results.animal_colors))
    animal_colors = par.Results.animal_colors;
end

if(~isempty(par.Results.IfColorAnimals))
    IfColorAnimals = par.Results.IfColorAnimals;
end

if(~isempty(par.Results.IfHistogram))
    IfHistogram = par.Results.IfHistogram;
end

if(~isempty(par.Results.plot_stats))
    plot_stats = par.Results.plot_stats;
end

if(~isempty(par.Results.BinWidth))
    BinWidth = par.Results.BinWidth;
end

if(~isempty(par.Results.displayopt))
    displayopt = par.Results.displayopt;
end

if(~isempty(par.Results.Normalization))
    Normalization = par.Results.Normalization;
end

if(~isempty(par.Results.DisplayStyle))
    DisplayStyle = par.Results.DisplayStyle;
end
% plot each field in value
hold on
fields = fieldnames(values);
num_plot = numel(fields);

% initialise
max_size = 0;
for m = 1:num_plot
    temp_value = extractfield(values,char(fields{m}))';
    this_size = numel(temp_value);
    if(this_size>max_size)
        max_size = this_size;
    end
end
mat_all_values = nan(max_size,num_plot);
num_samples = zeros(1,num_plot);
for m =1:num_plot
    temp_value = extractfield(values,char(fields{m}))';
    if(size(temp_value,2)>1)
        temp_value = temp_value';   % rotate to row vector
    end
    %     animal_value = arrayfun(@(x)nanmean(extractfield(values(x),fields{m})),1:size(values,2));
    mat_all_values(1:length(temp_value),m)= temp_value;
    num_samples(m) = length(temp_value);
end

% % stat stest on all combinations of fields
if num_plot>1
    all_cmp = nchoosek(1:num_plot,2);
else
    all_cmp =[1 1];
end
stats = struct();
max_value = max(mat_all_values(:));
this_plot_stats_hight = max_value;
plot_stats_step = 0.05*max_value;

if(~(strcmp(test_type,'anova')||(strcmp(test_type,'KW'))))
    % compare pairs - not corrected for multi comparison
    for c = 1:size(all_cmp,1)
        idx1 = all_cmp(c,1);
        idx2 = all_cmp(c,2);
        stats(c).field1 = fields{idx1};
        stats(c).field2 = fields(idx2);
        if(isempty(test_type))
            [P,h] = ranksum(mat_all_values(:,idx1),mat_all_values(:,idx2));
        else
            switch test_type
                case 'ttest'
                    [h,P]=ttest(mat_all_values(:,idx1),mat_all_values(:,idx2));
                case 'ranksum'
                    [P,h] = ranksum(mat_all_values(:,idx1),mat_all_values(:,idx2));
                case 'ttest2'
                    [h,P]=ttest2(mat_all_values(:,idx1),mat_all_values(:,idx2),'tail',tail);
                case 'signrank'
                    [P,h] = signrank(mat_all_values(:,idx1),mat_all_values(:,idx2),'tail',tail);
                case 'KS'
                    [h,P] = kstest2(mat_all_values(:,idx1),mat_all_values(:,idx2));
                    %                 IfHistogram = 1; % will plot histogram if chosen KS test
            end
        end
        stats(c).P = P;
        stats(c).h = h;
        
        % plot asterisks
        if(plot_stats)
            this_plot_stats_hight = this_plot_stats_hight+ plot_stats_step;
            plot([idx1 idx2],this_plot_stats_hight.*[1,1],'black');
            text(mean([idx1 idx2]),this_plot_stats_hight, get_stars_for_P(P), 'VerticalAlignment','bottom','HorizontalAlignment','center');
        end
        
    end
else
    % compare pairs - not corrected for multi comparison
    [~,~,mult_stats] = multcmp_stats_on_struct( values, 'test_type',test_type,'displayopt',displayopt );
   C = multcompare(mult_stats,'Display',displayopt);
    for c = 1:size(C,1)
       stats(c).P = C(c,end);
        idx1 = C(c,1);
        idx2 = C(c,2);
        % plot asterisks
        if(plot_stats)
            this_plot_stats_hight = this_plot_stats_hight+ plot_stats_step;
            plot([idx1 idx2],this_plot_stats_hight.*[1,1],'black');
            text(mean([idx1 idx2]),this_plot_stats_hight, get_stars_for_P( stats(c).P), 'VerticalAlignment','bottom','HorizontalAlignment','center');
        end
   end
end




if(isempty(test_type))
    [P,h] = ranksum(mat_all_values(:,1),mat_all_values(:,end));
else
    switch test_type
        case 'ttest'
            [h,P]=ttest(mat_all_values(:,1),mat_all_values(:,end));
        case 'ranksum'
            [P,h] = ranksum(mat_all_values(:,1),mat_all_values(:,end));
        case 'ttest2'
            [h,P]=ttest2(mat_all_values(:,1),mat_all_values(:,end));
        case 'signrank'
            [P,h] = signrank(mat_all_values(:,1),mat_all_values(:,end),'tail',tail);
    end
end


mean_values = nanmean(mat_all_values,1);
diff_value = 100.*(mean_values(:,end)-mean_values(:,1))./abs(mean_values(:,1));
se_values = nanstd( mat_all_values,0,1)/sqrt(size(mat_all_values,1));
sd_values = nanstd( mat_all_values,0,1);
sk_values = skewness(mat_all_values);
median_values = nanmedian(mat_all_values,1);
%% plot
x_positions = (1:num_plot).*x_interval+x_shift;
if ~IfPlotPrePostOnly   % plot all conditions
    if((~IfBoxWhisker)&&(~IfHistogram)) % mean +- se
        hold on
        if(connect_scatter)
            plot(repmat(x_positions,size(mat_all_values,1),1)',mat_all_values','color',[.7 .7 .7],'linewidth',0.5)
        end
        for i =1:num_plot
            temp_value = extractfield(values,char(fields{i}))';
            barwitherr(sd_values(i),mean_values(i),'XData',x_positions(i),'barwidth',barwidth,'FaceColor','none','EdgeColor',color_lut(i,:),'LineWidth',1);
            hold on
            if(~IfColorAnimals)
                if(~IfAddJitter)
                    scatter(x_positions(i).*ones(size(temp_value)),temp_value,50,'MarkerFaceColor','none','MarkerEdgeColor',color_lut(i,:));
                else
                    jitters = randi([-200 200],size(temp_value))./200.*x_interval.*0.1;
%                     jitters = linspace(-300,300,length(temp_value))'./200.*x_interval.*0.1;
                    scatter(x_positions(i).*ones(size(temp_value))+jitters,temp_value,15,'MarkerFaceColor',color_lut(i,:),'MarkerEdgeColor','none');
                end
            else
                for j = 1:numel(temp_value)
                    scatter(x_positions(i),temp_value(j),50,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',animal_colors(j,:),'MarkerFaceAlpha',0.5)
                end
            end
        end
        field_names = cellfun(@(x)strrep(x,'_',' '),fields,'UniformOutput',false);
        set(gca,'XTick',[x_positions(1)-x_interval, x_positions, x_positions(end)+x_interval],'XTickLabel',{'' char(field_names) ''});
    elseif IfBoxWhisker % plot interquantile range
        boxplot(mat_all_values,'Whisker',1,'colors',color_lut)
         field_names = cellfun(@(x)strrep(x,'_',' '),fields,'UniformOutput',false);
        set(gca,'XTick',[x_positions(1)-x_interval, x_positions, x_positions(end)+x_interval],'XTickLabel',{'' char(field_names) ''});
    elseif IfHistogram
        hold on
         for i =1:num_plot
             temp_value = extractfield(values,char(fields{i}))';
             if(strcmp(DisplayStyle,'stairs'))
                 hp(i) = histogram(temp_value,'Normalization',Normalization,'DisplayStyle',DisplayStyle,'facecolor','none', 'edgecolor',color_lut(i,:),'linewidth',2,'Binwidth',BinWidth);
             else
                 hp(i) = histogram(temp_value,'Normalization',Normalization,'DisplayStyle',DisplayStyle,'facecolor',color_lut(i,:), 'edgecolor','none','Binwidth',BinWidth,'FaceAlpha',1);
             end
         end  
         
        legend(hp,fields)
    end
    
    
    
else  % plot only the first and last set of values
    barwitherr(se_values(1),mean_values(1),'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',.5)
    barwitherr(se_values(end),mean_values(end),'FaceColor',[1 1 1],'EdgeColor',[0 0 0],'LineWidth',.5)
    
    num_plot = 2;
    plot([mat_all_values(:,1), mat_all_values(:,3)]','color',[.7 .7 .7],'linewidth',0.5)
    temp_value = values.pre';
    scatter(ones(1,length(temp_value)),temp_value,50,'MarkerFaceColor',color_lut(1,:),'MarkerEdgeColor',color_lut(i,:),...
        'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
    temp_value = values.post';
    scatter(2.*ones(1,length(temp_value)),temp_value,50,'MarkerFaceColor',color_lut(2,:),'MarkerEdgeColor',color_lut(i,:),...
        'MarkerFaceAlpha',.3,'MarkerEdgeAlpha',.3);
    set(gca,'XTick',0:3,'XTickLabel',{'' 'pre' 'post' ''});
    
end

if(~IfHistogram)
    ylabel(plot_name)
else
    plot_name = strrep(plot_name,'_',' ')
    title(plot_name)
    ylabel(Normalization)
end

% show stats if only two fields are compared
if(num_plot == 2)
    xlabel({[ test_type ' P=' num2str(P,'%10.1e') ' h= ' num2str(h)];...
        ['Difference: ' num2str(diff_value,'%10.3f')];...
        ['Mean: '  num2str(mean_values(1),'%10.3f') '; '  num2str(mean_values(end),'%10.3f')];...
        ['Median:' num2str(median_values(1),'%10.3f') '; '  num2str(median_values(end),'%10.3f')]
        ['SE: '  num2str(se_values(1),'%10.3f') '; '  num2str(se_values(end),'%10.3f') ];...
        ['SD: '  num2str(sd_values(1),'%10.3f') '; '  num2str(sd_values(end),'%10.3f') ];...
        ['Skewness: '  num2str(sk_values(1),'%10.3f') '; '  num2str(sk_values(end),'%10.3f')];...
        ['#Samples: ' num2str(num_samples)]})
else
    xlabel({
        ['Mean:  ' num2str(mean_values,'%10.3f')];...
        ['Median:' num2str(median_values,'%10.3f')];...
        ['SE:    '  num2str(se_values,'%10.3f')];...
        ['SD: '  num2str(sd_values,'%10.3f')];...
        ['Skew:  '  num2str(sk_values,'%10.3f')];...
        ['#Samples: ' num2str(num_samples)]})
    
end
set(gcf,'color','w')


y_pos = double(1.1*max(max(mat_all_values)));

end
