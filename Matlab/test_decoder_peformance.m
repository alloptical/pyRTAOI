function [test_decod_struct] = test_decoder_peformance(proj_struct,test_opt,decoder_name)
num_compares = round(numel(test_opt.fd_names)/2);

for i = 1:num_compares
    IF_REVERSE = 0;
    this_correct_fd = test_opt.fd_names{2*(i-1)+1};
    this_incorrect_fd = test_opt.fd_names{2*(i-1)+2};
    this_cmp_fds = {this_correct_fd,this_incorrect_fd};
    this_stim_type = strsplit(this_correct_fd,'_');
    this_stim_type = this_stim_type{2};
    % reverse for stim type2
    if strcmp(this_stim_type,'2')||strcmp(this_stim_type,'3')
        IF_REVERSE = 1;
    end

        [  test_decod_struct{i} ] =  get_binary_decoder_disc_time( proj_struct, struct(),...
            this_cmp_fds,test_opt,'IF_FRAMEWISE',0,'threshold',0,'IF_REVERSE',IF_REVERSE);
        test_decod_struct{i}.correct_fd = this_correct_fd;
        test_decod_struct{i}.incorrect_fd = this_incorrect_fd;
        
        figure('name',[decoder_name ' decoder performance tests'],'position',[100 100 500 500])
        plot_binary_decoder(test_decod_struct{i},test_opt)
        suptitle([{[decoder_name ' decoder : ']}; strrep(this_cmp_fds{1},'_','-'),' vs ' strrep(this_cmp_fds{2},'_','-')])
        

end
end

