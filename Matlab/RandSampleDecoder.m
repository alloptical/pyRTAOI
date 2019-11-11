% train decoder on randomly sammpled cells
% resample based on value in resamp_fd
resamp_fd = 'ipsi_corr';
resamp_cell_idx = [];
%% fa
[traces_in,fa_trial_idx,num_trials] = get_input_seq(cell_struct,resmp_cell_idx,...
    fa_opt.fd_names,fa_opt.bin_size,'IF_MEDFILT',fa_opt.IF_MEDFILT);%%
fa_struct = struct();
fa_struct.mean = mean(traces_in);
fa_struct.std = std(traces_in);
[fa_struct.lambda,fa_struct.psi,fa_struct.T,fa_struct.stats,fa_struct.F] = factoran(traces_in,fa_opt.m,'Xtype','data','Maxit',1000);
invsqrtPsi = diag(1 ./  sqrt(fa_struct.psi)); % get transition matrix (multiply to zscored data)
fa_struct.transmat = invsqrtPsi/(fa_struct.lambda'*invsqrtPsi);
%% decoder
