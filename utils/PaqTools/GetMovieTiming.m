function varargout=GetMovieTiming(varargin)

filename=varargin{1};

try
    InitTime=varargin{2};
catch
    InitTime=0;
end

try
    NumHeaderRows=varargin{3};
catch
    NumHeaderRows=1;
end

MovieTime = dlmread(filename,'\t',NumHeaderRows,0);
MovieTime = MovieTime(:,3);
MovieTime = MovieTime + InitTime;

varargout{1}=MovieTime;