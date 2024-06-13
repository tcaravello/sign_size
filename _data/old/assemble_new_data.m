clear all, close all, clc
%%
filename = 'macro_vars';
table_data_fred = readtable(filename);
table_proxy = readtable('oilSupplyNewsShocks_2022M06', Sheet = 'Monthly');
%%
% define variables explained in Kaenzig
proxy = table_proxy{:,2}; %shock
svar_proxy = table_proxy{:,3};
DATES = table_proxy{:,1};
CPI = table_data_fred{:,2}; %cpi
CPI_CORE = table_data_fred{:,3}; %core cpi
FF = table_data_fred{:,4}; %fed funds rate
M3RATE = table_data_fred{:,5};
Y10RATE = table_data_fred{:,6};
IP = table_data_fred{:,7};
UNEMP = table_data_fred{:,8};
RPCE = table_data_fred{:,9};
WAGES = table_data_fred{:,10};
POIL = table_data_fred{:,11};


%%
save('oil_data_new')






