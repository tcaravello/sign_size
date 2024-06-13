%% Create table from the simulation.

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

%% Load
load _results\table_simulations_goncalves.mat  collector_lin collector_sign_1 collector_sign_2 collector_sign_3 collector_size_1 collector_size_2 collector_size_3 extra_stuff_collect

collector_total = zeros([size(collector_sign_1),7]);
collector_total(:,:,:,:,:,:,1) = collector_lin;
collector_total(:,:,:,:,:,:,2) = collector_sign_1;
collector_total(:,:,:,:,:,:,3) = collector_sign_2;
collector_total(:,:,:,:,:,:,4) = collector_sign_3;
collector_total(:,:,:,:,:,:,5) = collector_size_1;
collector_total(:,:,:,:,:,:,6) = collector_size_2;
collector_total(:,:,:,:,:,:,7) = collector_size_3;


% first dim: 10,5,1 significance
% second dim: horizons
% third dim: sample sizes
% fourth dim: variables
% fifth: type of non-linearity in the DGP
% sixth: order of non-linearity in the DGP
% seventh: regressor

% for extra_stuff_collect: rows of the structure indicate the type of
% non-linearity in the DGP, the columns indicate the order of the
% polynomial used.
%% Plot Histogram of SE for the sign 3, size 3 case.

s_use = extra_stuff_collect{2,3}; %DGP: sign 3;
series_plot = s_use.t_collector_size_3(:,1,2);
%%
mean(series_plot)
quants = 1.64;
algo1 = mean(series_plot>=quants);
algo2 = mean(series_plot<=quants);

%% Create Latex table
sig_level = 1; % 1 for 10%, 2 for 5%, 3 for 10%
hors_report = [0;6;12;24;48];
sample_size = 1;
var_use = 2; %order of variable to get table for.
types_use = {'Linear','Sign','Size'};
reg_use = {'Linear', 'Sign - 1', 'Sign - 2', 'Sign - 3', 'Size - 1', 'Size - 2', 'Size - 3'};
orders_dgp_use = [1,2,3];
clc
fprintf('\\begin{table}[h!] \n')
fprintf('\\centering \n')
fprintf('\\begin{tabular}{|cc|ccccccccc|} \\hline \n')
fprintf('Regressor & Horizon &  \\multicolumn{9}{|c|}{DGP} \\\\ \n')
fprintf('- Order &  & Linear & & Sign - 1 & Sign - 2 & Sign - 3 & & Size - 1 & Size - 2 & Size - 3 \\\\ \\hline \n')
for i_reg = 1:size(collector_total,7)
        fprintf(reg_use{i_reg})
        for i_h = 1:length(hors_report)
        hor_aux = hors_report(i_h);
        fprintf('& $%2.0f$ & $%1.3f$ & & $%1.3f$ & $%1.3f$ & $%1.3f$ & & $%1.3f$ & $%1.3f$ & $%1.3f$ \\\\ \n',...
            hor_aux, collector_total(sig_level,i_h,sample_size,var_use,1,1,i_reg),...
            collector_total(sig_level,i_h,sample_size,var_use,2,1,i_reg),...
            collector_total(sig_level,i_h,sample_size,var_use,2,2,i_reg),...
            collector_total(sig_level,i_h,sample_size,var_use,2,3,i_reg),...
            collector_total(sig_level,i_h,sample_size,var_use,3,1,i_reg),...
            collector_total(sig_level,i_h,sample_size,var_use,3,2,i_reg),...
            collector_total(sig_level,i_h,sample_size,var_use,3,3,i_reg));
        end
        fprintf(' \\multicolumn{11}{|c|}{} \\\\ \n')
        if i_reg == 1 | i_reg == 4 | i_reg == 7
        fprintf('\\hline \n')
        end
end
fprintf('\\end{tabular} \n')
fprintf('\\caption{} \n')
fprintf('\\end{table}')