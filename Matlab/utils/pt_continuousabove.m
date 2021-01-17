function varargout=pt_continuousabove(data,baseline,abovethresh,mintime,maxtime,ignorepts)
% function varargout=continuousabove(data,baseline,abovethresh,mintime,maxtime,ignorepts);
%
% Finds periods in a linear trace that are above some baseline by some
% minimum amount (abovethresh) for between some minimum amount of time
% (mintime) and some maximum amount of time (maxtime).  Output is the
% indices of start and stop of those periods in data.

above=find(data>=baseline+abovethresh);

if isempty(above)||isempty(diff(diff(above)))
    aboveperiods=[];
elseif max(abs(diff(diff(above))))==0 && length(above)>=mintime && length(above)<=maxtime %if only 1 potential upstate found
    aboveperiods = [above(1) above(end)];
elseif ~isempty(above);%if many possible upstates
	ends=find(diff(above)~=1);%find breaks between potential upstates
    ends(end+1)=0;%for the purposes of creating lengths of each potential upstate
    ends(end+1)=length(above);%one of the ends comes at the last found point above baseline
    ends=sort(ends);
    if size(ends,2)==3%handle the rare case when there are two ends and the vector needs to be transposed AP 20090710
        ends=ends';
    end
    lengths=diff(ends);%length of each potential upstate
    ends(1)=[];%lose the 0 added before
    
    %% segment allows for ignoring small interruptions in above periods.
    % interruptions must be equal to or less than the value ignorepts
    % given as a function input
    e2=above(ends);
    l2=lengths-1;
    prelimstops = e2;
    
    try
        prelimstarts = e2-l2;
    catch
        prelimstarts = e2-l2';
    end
    if size(prelimstarts,1)~=size(e2,1) % 20190121 ZZ added for Matlab2017
        prelimstarts = e2-l2';
    end
    
    intereventintervals = prelimstarts(2:end)-prelimstops(1:end-1);
    shortenough = find(intereventintervals<=(ignorepts+1));
    for sidx = length(shortenough):-1:1;
        this = shortenough(sidx);
        lengths(this) = lengths(this) + lengths(this+1) + intereventintervals(this)-1;
        lengths(this+1) = [];
        
        ends(this) = [];
    end
    %%
    
    good=find(lengths>=mintime & lengths<=maxtime);%must be longer than 500ms but shorter than 15sec
    e3=reshape(above(ends(good)),[length(good) 1]);
    l3=reshape(lengths(good)-1,[length(good) 1]);
    aboveperiods(:,2)=e3;%upstate ends according to the averaged reading
    aboveperiods(:,1)=e3-l3;%upstate beginnings according to averaged reading
else
    aboveperiods=[];
end

varargout{1}=aboveperiods;
if nargout==2;
    if isempty(aboveperiods);
        belowperiods = [];
    else
        stops=aboveperiods(:,2);
        stops(2:end+1,1)=stops;
        stops(1)=0;
        stops=stops+1;

        starts=aboveperiods(:,1);
        starts(end+1,1)=length(data)+1;
        starts=starts-1;

        belowperiods=cat(2,stops,starts);

        if belowperiods(1,2)<=belowperiods(1,1);
            belowperiods(1,:)=[];
        end
        if belowperiods(end,2)<=belowperiods(end,1);
            belowperiods(end,:)=[];
        end
    end

    varargout{2}=belowperiods;
end