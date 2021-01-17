function [ r,hd,p_fit,s,y_fit,delta_fit,p,rr ] = plot_fit_linear( x,y,zz,color,varargin )
n = 1;
showError = 0;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'degree')
       n = varargin{v+1};
   elseif strcmpi(varargin{v},'showError')
        showError = varargin{v+1};
    end
end
    [p_fit,s] = polyfit(x,y,n);
    rr = 1 - s.normr^2 / norm(y-mean(y))^2;
    % y±delta contains at least 50% of the predictions of future observations at x
    [y_fit,delta_fit] = polyval(p_fit,zz,s);
    lm = fitlm(x,y,'linear');
    p = coefTest(lm);
    
    r = corrcoef(x,y); r = r(1,2);
    % plot 
    zz = sort(zz);
    [y_plot,delta_plot] = polyval(p_fit,zz,s);
    hd = plot(zz,y_plot,'Color',color);
    % plot error range
    if(showError)
        tempy = [y_plot+delta_plot];
        plot([zz(1) zz(end)],[tempy(1) tempy(end)],'Color',color,'linestyle','--');
        tempy = [y_plot-delta_plot];
        plot([zz(1) zz(end)],[tempy(1) tempy(end)],'Color',color,'linestyle','--');
    end

end

