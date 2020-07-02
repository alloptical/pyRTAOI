function [ mean_dff_sta_traces,dff_sta,mean_sta_traces,raw_sta_traces,stim_frames,df_sta,mean_df_sta] = make_sta_from_traces( traces,frames_with_stim,pre_frames,post_frames,baseline_frames )
%
if(size(traces,2)<size(traces,1))
    traces = traces';
end

mean_dff_sta_traces = [];
dff_sta = [];
mean_sta_traces = [];
df_sta = [];
mean_df_sta = [];
raw_sta_traces = [];
stim_frames = [];


if (isempty(frames_with_stim))
    warning('no stim found, make_sta_from_traces return')
    return 
end

trace_length = size(traces,2);
stim_frames = frames_with_stim((frames_with_stim<trace_length-post_frames));
raw_sta_traces = make_sta_traces(traces, stim_frames, pre_frames, post_frames);

% get baselined sta traces
sta_baseline = mean (raw_sta_traces(:,:,baseline_frames),3);
% sta_baseline = min (raw_sta_traces(:,:,baseline_frames),[],3);

for m = 1:size(sta_baseline,1)
    for n = 1: size(sta_baseline,2)
        dff_sta(m,n,1:size(raw_sta_traces,3)) = (raw_sta_traces(m,n,:)-sta_baseline(m,n))./sta_baseline(m,n);
        df_sta(m,n,1:size(raw_sta_traces,3)) = (raw_sta_traces(m,n,:)-sta_baseline(m,n));

    end
end

mean_dff_sta_traces = squeeze(nanmean(dff_sta, 2));  % targets trace mean, averaged over stims
mean_sta_traces = squeeze(nanmean(raw_sta_traces, 2));
mean_df_sta = squeeze(nanmean(df_sta, 2));
dff_sta = squeeze(dff_sta);
df_sta = squeeze(df_sta);
raw_sta_traces = squeeze(raw_sta_traces); % added 20190803
% median_dff_sta_traces = squeeze(nanmedian(dff_sta, 2));


end

