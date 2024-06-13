%% Simulate from simple DGP

% Tomas E. Caravello and Pedro Martinez Bruera

% this vesion June 13, 2024
%% HOUSEKEEPING
clear all
close all
clc

%change the path to the correponding one in your computer.

path = '/Users/tomyc/Dropbox (MIT)/asymmetric LP/code/GitHub/sign_size';
vintage = '';
task = '/simulation';

addpath([path vintage task])
addpath([path vintage '/_auxiliary_functions'])
cd([path vintage task])

% initialize random number generator to be able to replicate results exactly
rng default

% Set text interpreter for figures to latex
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

save_sim = 0; % 1 if you want to save simulation results, 0 otherwise.

%% Calibration
global prop_big coef_linear coef_non_linear cut_off sigma_oil T_burn N_sim

global hor_coefs var_x_mom prop_zeros hor_lp_max N_report coefs_y_1_aux coefs_y_2_aux
global  opt T_to_use n_lags_lp hors_report

T_burn = 300; %truncation horizon for IRFs
T_to_use = 400; %months
hors_report = [0;6;12;24;48]+1; %horizons to report

N_T = size(T_to_use,1);
N_sim = 5000; %number of simulations 
sigma_oil = 1; %standard deviation of the oil shock

% Simualtion settings

type = 'size';
distribution_simul = 'normal'; % 'normal', or 'exponential';
order_sim = 1; %order of the true DGP

% coefficients on the true DGP
coef_linear = 1;
coef_non_linear = 1;

prop_big = 0.4;
prop_zeros = 0.0;
cut_off  = -norminv(prop_big/2);

% LP estimation settings
opt = struct();
opt.transformation = 1; %order to use
opt.bound = sigma_oil*cut_off;
opt.order = 1:2;
opt.har = 1;
n_lags_lp = 1; %number of lags in the local projections
hor_lp_max = max(hors_report);
N_report = length(hors_report);
%%

hor_coefs = (0:1:T_burn-1)';
n_x_mom  = 1000000;

% generate a simulated sample to compute moments

obs_shock = rand(1000000,1);
obs_shock_2 = obs_shock>prop_zeros;

if strcmp(distribution_simul,'normal')
    x_moments_aux = sigma_oil * randn(n_x_mom,1);
    x_moments = obs_shock_2.*x_moments_aux;
elseif strcmp(distribution_simul,'exponential')
    x_moments_aux = sigma_oil * exprnd(ones(n_x_mom,1))-1;
    x_moments = obs_shock_2.*x_moments_aux;
end


if strcmp(type,'sign')
if order_sim == 1
    f_shock = coef_linear * x_moments + coef_non_linear * abs(x_moments);
else
    f_shock = coef_linear * x_moments + coef_non_linear*(abs(x_moments).^(order_sim));
end
elseif strcmp(type,'size')
    if order_sim == 1
        bound = cut_off*sigma_oil;
        fun_size = @(x) ((x>bound).*(x-bound)+(x<-bound).*(x+bound));
        f_shock = coef_linear * x_moments + coef_non_linear * fun_size(x_moments);
    else
        f_shock = coef_linear * x_moments + coef_non_linear * x_moments.*(abs(x_moments).^(order_sim-1));
    end
else
    f_shock = x_moments;
end

var_x_mom = var(f_shock);
%% DGP
global C_mat
rho = 0.95;
a = 1/sqrt(var_x_mom);
rho_e = 0.3 * a;
C_mat = [a 0;
        rho_e/a sqrt(1-(rho_e/a)^2)];

A_mat = [rho, 0;
        0.5, rho ];

coefs_y_1_aux = A_mat(1,1).^(hor_coefs);
coefs_y_2_aux = A_mat(2,2).^(hor_coefs);


% get IRFs on the first shock
irf_population = zeros(2,T_burn);
irf_population(:,1) = C_mat * [1;0];
for tt = 2:T_burn
irf_population(:,tt) = A_mat * irf_population(:,tt-1);
end

%%

types_use = {'linear','sign','size'};
orders_dgp_use = [1,2,3];

n_types = length(types_use);
n_orders = length(orders_dgp_use);

collector_sign_1 = zeros(3,N_report, length(T_to_use),2,n_types,n_orders);
collector_sign_2 = zeros(3,N_report, length(T_to_use),2,n_types,n_orders);
collector_sign_3 = zeros(3,N_report, length(T_to_use),2,n_types,n_orders);

collector_size_1 = zeros(3,N_report, length(T_to_use),2,n_types,n_orders);
collector_size_2 = zeros(3,N_report, length(T_to_use),2,n_types,n_orders);
collector_size_3 = zeros(3,N_report, length(T_to_use),2,n_types,n_orders);

collector_lin =  zeros(3,N_report, length(T_to_use),2,n_types,n_orders);

extra_stuff_collect = cell(3,3);
for i_type = 1:n_types
    if  strcmp(types_use{i_type},'linear')
        [collector_sign_1(:,:,:,:,i_type,1), collector_sign_2(:,:,:,:,i_type,1), collector_sign_3(:,:,:,:,i_type,1), ...
            collector_size_1(:,:,:,:,i_type,1), collector_size_2(:,:,:,:,i_type,1), collector_size_3(:,:,:,:,i_type,1),...
            collector_lin(:,:,:,:,i_type,1), extra_stuff_collect{i_type,1}]  = simulate_var_dgp(types_use{i_type},orders_dgp_use(1));
    else
        for i_order = 1:n_orders
                [collector_sign_1(:,:,:,:,i_type,i_order), collector_sign_2(:,:,:,:,i_type,i_order), collector_sign_3(:,:,:,:,i_type,i_order), ...
            collector_size_1(:,:,:,:,i_type,i_order), collector_size_2(:,:,:,:,i_type,i_order), collector_size_3(:,:,:,:,i_type,i_order),...
            collector_lin(:,:,:,:,i_type,i_order),extra_stuff_collect{i_type,i_order}]  = simulate_var_dgp(types_use{i_type},orders_dgp_use(i_order));
        end
    end
end

%% Save stuff
if save_sim == 1
save _results\table_simulations_goncalves.mat collector_lin collector_sign_1 collector_sign_2 collector_sign_3 collector_size_1 collector_size_2 collector_size_3 extra_stuff_collect
end
