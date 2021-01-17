function [ r,f,gof] = plot_fit_exp( x,y,varargin )
n = 1;
for v = 1:numel(varargin)
    if strcmpi(varargin{v},'degree')
       n = varargin{v+1};
    end
end
    [f,gof] = fit(x,y,['exp' num2str(n)]);
    plot(f,x,y,'oblack')
    r = gof.rsquare;
end

