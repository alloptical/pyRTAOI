%% 
shuf_accuracy = cell2mat(stim_struct.shuf_classif_accuracy');
num_shuf = size(shuf_accuracy,2);
num_frames = size(shuf_accuracy,1);
figure
hold on
for f = 1:num_frames
    scatter(ones(1,num_shuf).*f,shuf_accuracy(f,:),'MarkerEdgeColor','white','MarkerFaceColor',[.5 .5 .5])
end
title('stim')
%%
figure;
hold on
plot(stim_struct.framewise_hr,'color','black')
plot(stim_struct.framewise_fa,'color',[.5 .5 .5])
ylim([0 1])
title('stim decoder')

%%
figure;
hold on
plot(choice_struct.framewise_hr,'color','black')
plot(choice_struct.framewise_fa,'color',[.5 .5 .5])
ylim([0 1])
title('choice decoder')

