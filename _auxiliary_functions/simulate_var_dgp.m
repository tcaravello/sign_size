function  [report_sign_1, report_sign_2, report_sign_3, report_size_1, report_size_2, report_size_3, report_lin, extra_stuff]  = simulate_var_dgp(type,order_sim)

global T_to_use coefs_y_1_aux coefs_y_2_aux T_burn C_mat n_lags_lp hor_lp_max N_sim N_report opt hors_report

opt_1 = opt;
opt_2 = opt;
opt_3 = opt;

opt_1.transformation = 1;
opt_2.transformation = 2;
opt_3.transformation = 3;

% collector
irf_collector_sign_1 = zeros(N_sim,N_report,2, length(T_to_use));
se_collector_sign_1 = zeros(N_sim,N_report,2, length(T_to_use));
t_collector_sign_1 = zeros(N_sim,N_report,2, length(T_to_use));
irf_collector_sign_2 = zeros(N_sim,N_report,2, length(T_to_use));
se_collector_sign_2 = zeros(N_sim,N_report,2, length(T_to_use));
t_collector_sign_2 = zeros(N_sim,N_report,2, length(T_to_use));
irf_collector_sign_3 = zeros(N_sim,N_report,2, length(T_to_use));
se_collector_sign_3 = zeros(N_sim,N_report,2, length(T_to_use));
t_collector_sign_3 = zeros(N_sim,N_report,2, length(T_to_use));


irf_collector_size_1 = zeros(N_sim,N_report,2, length(T_to_use));
se_collector_size_1 = zeros(N_sim,N_report,2, length(T_to_use));
t_collector_size_1 = zeros(N_sim,N_report,2, length(T_to_use));
irf_collector_size_2 = zeros(N_sim,N_report,2, length(T_to_use));
se_collector_size_2 = zeros(N_sim,N_report,2, length(T_to_use));
t_collector_size_2 = zeros(N_sim,N_report,2, length(T_to_use));
irf_collector_size_3 = zeros(N_sim,N_report,2, length(T_to_use));
se_collector_size_3 = zeros(N_sim,N_report,2, length(T_to_use));
t_collector_size_3 = zeros(N_sim,N_report,2, length(T_to_use));


irf_collector_lin = zeros(N_sim,N_report,2, length(T_to_use));
se_collector_lin = zeros(N_sim,N_report,2, length(T_to_use));
t_collector_lin = zeros(N_sim,N_report,2, length(T_to_use));


% report

report_size_1 = zeros(3,N_report, length(T_to_use),2);
report_size_2 = zeros(3,N_report, length(T_to_use),2);
report_size_3 = zeros(3,N_report, length(T_to_use),2);

report_sign_1 = zeros(3,N_report, length(T_to_use),2);
report_sign_2 = zeros(3,N_report, length(T_to_use),2);
report_sign_3 = zeros(3,N_report, length(T_to_use),2);
report_lin = zeros(3,N_report, length(T_to_use),2);


for i_ts = 1:length(T_to_use)
T_s = T_to_use(i_ts); %sample size
quants = [tinv(0.95,T_s - 3); tinv(0.975,T_s - 3); tinv(0.995,T_s - 3)];
coefs_y_1 = [coefs_y_1_aux; zeros(T_s,1)];
coefs_y_2 = [coefs_y_2_aux; zeros(T_s,1)];
shocks_2 = randn(T_burn+T_s,N_sim);
mat_coefs_y_1=  toeplitz(coefs_y_1, zeros(T_s+T_burn,1));
mat_coefs_y_2=  toeplitz(coefs_y_2, zeros(T_s+T_burn,1));

[shocks_1,f_shock_1] = generate_shocks(type,order_sim,T_s);

y_data_1_aux = mat_coefs_y_1 * C_mat(1,1) * f_shock_1;
y_data_2_aux = mat_coefs_y_2 * (0.5 * [zeros(1,N_sim);y_data_1_aux(2:end,:)] + C_mat(2,1) * f_shock_1 + C_mat(2,2) * shocks_2);
% generate series T_s to regress.
for i_sim = 1:N_sim
    if floor(i_sim/100) == i_sim/100
        fprintf('Starting iter.. %4.0f ',i_sim)
        tic
    end
    y_data = [y_data_1_aux(T_burn+1:end,i_sim),y_data_2_aux(T_burn+1:end,i_sim)];

% Size

% order 1
[irf_coefs_size_1,irf_se_size_1,irf_t_stats_size_1,~,~] =...
    run_non_linear_local_projections(y_data,shocks_1(:,i_sim),y_data,n_lags_lp,hor_lp_max, 'size', opt_1);
% order 2
[irf_coefs_size_2,irf_se_size_2,irf_t_stats_size_2,~,~] =...
    run_non_linear_local_projections(y_data,shocks_1(:,i_sim),y_data,n_lags_lp,hor_lp_max, 'size', opt_2);
%order 3
[irf_coefs_size_3,irf_se_size_3,irf_t_stats_size_3,~,~] =...
    run_non_linear_local_projections(y_data,shocks_1(:,i_sim),y_data,n_lags_lp,hor_lp_max, 'size', opt_3);

% Sign

% order 1
[irf_coefs_sign_1,irf_se_sign_1,irf_t_stats_sign_1,~,~] =...
    run_non_linear_local_projections(y_data,shocks_1(:,i_sim),y_data,n_lags_lp,hor_lp_max, 'sign', opt_1);
% order 2
[irf_coefs_sign_2,irf_se_sign_2,irf_t_stats_sign_2,~,~] =...
    run_non_linear_local_projections(y_data,shocks_1(:,i_sim),y_data,n_lags_lp,hor_lp_max, 'sign', opt_2);
%order 3
[irf_coefs_sign_3,irf_se_sign_3,irf_t_stats_sign_3,~,~] =...
    run_non_linear_local_projections(y_data,shocks_1(:,i_sim),y_data,n_lags_lp,hor_lp_max, 'sign', opt_3);

%Linear
[irf_coefs_lin,irf_se_lin,irf_t_stats_lin,~,~] =...
    run_non_linear_local_projections(y_data,shocks_1(:,i_sim),y_data,n_lags_lp,hor_lp_max, 'linear');

% Size

irf_collector_size_1(i_sim, :,:, i_ts) = irf_coefs_size_1([2 4],hors_report)';
se_collector_size_1(i_sim, :,:, i_ts) = irf_se_size_1([2 4],hors_report)';
t_collector_size_1(i_sim, :,:, i_ts) = irf_t_stats_size_1([2 4],hors_report)';

irf_collector_size_2(i_sim, :,:, i_ts) = irf_coefs_size_2([2 4],hors_report)';
se_collector_size_2(i_sim, :,:, i_ts) = irf_se_size_2([2 4],hors_report)';
t_collector_size_2(i_sim, :,:, i_ts) = irf_t_stats_size_2([2 4],hors_report)';

irf_collector_size_3(i_sim, :,:, i_ts) = irf_coefs_size_3([2 4],hors_report)';
se_collector_size_3(i_sim, :,:, i_ts) = irf_se_size_3([2 4],hors_report)';
t_collector_size_3(i_sim, :,:, i_ts) = irf_t_stats_size_3([2 4],hors_report)';

% Sign

irf_collector_sign_1(i_sim, :,:, i_ts) = irf_coefs_sign_1([2 4],hors_report)';
se_collector_sign_1(i_sim, :,:, i_ts) = irf_se_sign_1([2 4],hors_report)';
t_collector_sign_1(i_sim, :,:, i_ts) = irf_t_stats_sign_1([2 4],hors_report)';

irf_collector_sign_2(i_sim, :,:, i_ts) = irf_coefs_sign_2([2 4],hors_report)';
se_collector_sign_2(i_sim, :,:, i_ts) = irf_se_sign_2([2 4],hors_report)';
t_collector_sign_2(i_sim, :,:, i_ts) = irf_t_stats_sign_2([2 4],hors_report)';

irf_collector_sign_3(i_sim, :,:, i_ts) = irf_coefs_sign_3([2 4],hors_report)';
se_collector_sign_3(i_sim, :,:, i_ts) = irf_se_sign_3([2 4],hors_report)';
t_collector_sign_3(i_sim, :,:, i_ts) = irf_t_stats_sign_3([2 4],hors_report)';

% Linear

irf_collector_lin(i_sim, :,:, i_ts) = irf_coefs_lin([1 2],hors_report)';
se_collector_lin(i_sim, :,:, i_ts) = irf_se_lin([1 2],hors_report)';
t_collector_lin(i_sim, :,:, i_ts) = irf_t_stats_lin([1 2],hors_report)';

    if floor((i_sim+1)/100) == (i_sim+1)/100
        toc
        clc
    end

end


% Do tables with this.
for i_var=1:2
report_sign_1(:,:, i_ts,i_var) = [mean(abs(t_collector_sign_1(:,:,i_var, i_ts))>=quants(1)); mean(abs(t_collector_sign_1(:,:,i_var, i_ts))>=quants(2));  mean(abs(t_collector_sign_1(:,:,i_var, i_ts))>=quants(3))];
report_sign_2(:,:, i_ts,i_var) = [mean(abs(t_collector_sign_2(:,:,i_var, i_ts))>=quants(1)); mean(abs(t_collector_sign_2(:,:,i_var, i_ts))>=quants(2));  mean(abs(t_collector_sign_2(:,:,i_var, i_ts))>=quants(3))];
report_sign_3(:,:, i_ts,i_var) = [mean(abs(t_collector_sign_3(:,:,i_var, i_ts))>=quants(1)); mean(abs(t_collector_sign_3(:,:,i_var, i_ts))>=quants(2));  mean(abs(t_collector_sign_3(:,:,i_var, i_ts))>=quants(3))];

report_size_1(:,:, i_ts,i_var) = [mean(abs(t_collector_size_1(:,:,i_var, i_ts))>=quants(1)); mean(abs(t_collector_size_1(:,:,i_var, i_ts))>=quants(2));  mean(abs(t_collector_size_1(:,:,i_var, i_ts))>=quants(3))];
report_size_2(:,:, i_ts,i_var) = [mean(abs(t_collector_size_2(:,:,i_var, i_ts))>=quants(1)); mean(abs(t_collector_size_2(:,:,i_var, i_ts))>=quants(2));  mean(abs(t_collector_size_2(:,:,i_var, i_ts))>=quants(3))];
report_size_3(:,:, i_ts,i_var) = [mean(abs(t_collector_size_3(:,:,i_var, i_ts))>=quants(1)); mean(abs(t_collector_size_3(:,:,i_var, i_ts))>=quants(2));  mean(abs(t_collector_size_3(:,:,i_var, i_ts))>=quants(3))];

report_lin(:,:, i_ts,i_var) = [mean(abs(t_collector_lin(:,:,i_var, i_ts))>=quants(1)); mean(abs(t_collector_lin(:,:,i_var, i_ts))>=quants(2));  mean(abs(t_collector_lin(:,:,i_var, i_ts))>=quants(3))];

end

end

extra_stuff = struct();

extra_stuff.irf_collector_sign_1 = irf_collector_sign_1;
extra_stuff.se_collector_sign_1 = se_collector_sign_1;
extra_stuff.t_collector_sign_1 = t_collector_sign_1;
extra_stuff.irf_collector_sign_2 = irf_collector_sign_2;
extra_stuff.se_collector_sign_2 = se_collector_sign_2;
extra_stuff.t_collector_sign_2 = t_collector_sign_2;
extra_stuff.irf_collector_sign_3 = irf_collector_sign_3;
extra_stuff.se_collector_sign_3 = se_collector_sign_3;
extra_stuff.t_collector_sign_3 = t_collector_sign_3;

extra_stuff.irf_collector_size_1 = irf_collector_size_1;
extra_stuff.se_collector_size_1 = se_collector_size_1;
extra_stuff.t_collector_size_1 = t_collector_size_1;
extra_stuff.irf_collector_size_2 = irf_collector_size_2;
extra_stuff.se_collector_size_2 = se_collector_size_2;
extra_stuff.t_collector_size_2 = t_collector_size_2;
extra_stuff.irf_collector_size_3 = irf_collector_size_3;
extra_stuff.se_collector_size_3 = se_collector_size_3;
extra_stuff.t_collector_size_3 = t_collector_size_3;


extra_stuff.irf_collector_lin = irf_collector_lin ;
extra_stuff.se_collector_lin = se_collector_lin;
extra_stuff.t_collector_lin = t_collector_lin;

end