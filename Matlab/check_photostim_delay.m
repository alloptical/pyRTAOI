function [delayed_phototrial_idx] = check_photostim_delay(caiman_data,trials,paq_photo_frames,opt)
    sens_stim_frames = caiman_data.sensory_stim_frames+caiman_data.t_init;
    photo_stim_frames =  caiman_data.online_photo_frames + caiman_data.t_init;
    photo_stim_frames(photo_stim_frames>sens_stim_frames(end)+opt.trial_length)=[];
    [photo_trial_idx] = get_trials_with_photostim( caiman_data.sensory_stim_frames, caiman_data.online_photo_frames );

    dummy_trial_idx = find(trials.trialVar==2);
    [dummy_photostim_trial_idx,idx]=intersect(photo_trial_idx,dummy_trial_idx);
    true_photo_trial_idx = setdiff(photo_trial_idx,dummy_photostim_trial_idx);
    dummy_photo_stim_frames = photo_stim_frames(idx);
    photo_stim_frames = setdiff(photo_stim_frames,dummy_photo_stim_frames);
    movie_photo_stim_frames = photo_stim_frames - caiman_data.t_init; % for checking with raw movie
    movie_dummy_photo_stim_frames = dummy_photo_stim_frames - caiman_data.t_init;
    try
        oppo_photo_frames = caiman_data.online_oppo_photo_frames + caiman_data.t_init;
        movie_oppo_photo_frames = oppo_photo_frames - caiman_data.t_init; % for checking with raw movie
        [oppo_photostim_trial_idx] = get_trials_with_photostim( caiman_data.sensory_stim_frames, caiman_data.online_oppo_photo_frames );

    catch
        oppo_photo_frames = [];movie_oppo_photo_frames = [];oppo_photostim_trial_idx = [];
    end
    movie_frames_skipped = caiman_data.frames_skipped;

    
    % check caiman recorded photo frame vs actual photo frame
    pre_frames = 30;
    check_photo_trace = zeros(1,movie_frames_skipped(end));
    check_photo_trace(paq_photo_frames) = 1;
    [~,~,~,paq_photo] = make_sta_from_traces(check_photo_trace,movie_photo_stim_frames,pre_frames,150,1:opt.sta_baseline_frames);
    check_photo_trace = zeros(1,movie_frames_skipped(end));
    check_photo_trace(paq_photo_frames) = 1;
    [~,~,~,paq_dummyphoto] = make_sta_from_traces(check_photo_trace,movie_dummy_photo_stim_frames,pre_frames,150,1:opt.sta_baseline_frames);
    check_photo_trace = zeros(1,movie_frames_skipped(end));
    check_photo_trace(paq_photo_frames) = 1;
    [~,~,~,paq_oppophoto] = make_sta_from_traces(check_photo_trace,movie_oppo_photo_frames,pre_frames,150,1:opt.sta_baseline_frames);
    
    figure('name','photostim delay check','position',[200 200 1200 500]);
    subplot(1,4,1);hold on;imagesc(paq_photo); colormap('gray')
    plot([30 30],ylim,'color','r');xlabel('Frames'); title('Photostim triggers'); axis square
    subplot(1,4,2);hold on;imagesc(paq_dummyphoto); colormap('gray')
    plot([30 30],ylim,'color','r');xlabel('Frames'); title('Dummy Photostim triggers'); axis square   
    subplot(1,4,3);hold on;imagesc(paq_oppophoto); colormap('gray')
    plot([30 30],ylim,'color','r');xlabel('Frames'); title('Oppo Photostim triggers'); axis square 
    
    % get trials with very long delay 
    max_photo_delay = 10; % frames
    binwidth = 3;
    [~,photo_delays] = find(paq_photo == 1); photo_delays = photo_delays-pre_frames;
    [~,dummyphoto_delays] = find(paq_dummyphoto == 1); dummyphoto_delays = dummyphoto_delays-pre_frames;
    [~,oppophoto_delays] = find(paq_oppophoto == 1); oppophoto_delays = oppophoto_delays-pre_frames;
     subplot(1,4,4); hold on
     histogram(photo_delays,'displaystyle','stairs','edgecolor','r','binwidth',binwidth)
     histogram(dummyphoto_delays,'displaystyle','stairs','edgecolor', tint([1,0,0],.5),'binwidth',binwidth)
     histogram(oppophoto_delays,'displaystyle','stairs','edgecolor',[0 0 0],'binwidth',binwidth); axis square
     xlabel('Delay frames'); ylabel('Num. trials')
     delayed_trials_idx = sort([true_photo_trial_idx(photo_delays>max_photo_delay) dummy_photostim_trial_idx(dummyphoto_delays>max_photo_delay) oppo_photostim_trial_idx(oppophoto_delays>max_photo_delay)]);
     delayed_phototrial_idx = sort([true_photo_trial_idx(photo_delays>max_photo_delay) oppo_photostim_trial_idx(oppophoto_delays>max_photo_delay)]);
     disp(['delayed trial idx:' num2str(delayed_trials_idx)])
     if ~isempty(delayed_phototrial_idx) % put these trials into cheated
         first_bug_trial = delayed_trials_idx(1);
         disp(['First delayed trial:' num2str(delayed_trials_idx(1))])
         suptitle(['First delayed trial:' num2str(delayed_trials_idx(1)) ', stim type:' num2str(trials.stim_type(first_bug_trial)) 'Var' num2str(trials.trialVar(first_bug_trial))...
             ' photo:'  num2str(trials.photostim(first_bug_trial)) ' fa:' num2str(trials.fa(first_bug_trial))])
         
     end

end

