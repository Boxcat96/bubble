% MAIN FILE 
% Bubble Growth model 
% PaGq - RJ - TH
% Feb 10, 2017
% Revision August 2019
% Use model with stock prices

addpath('\04_Bubble\')

% close all
clear

% if you reload data, choose 1.
pick_reload = 1;

% 0. Get data and Figure 7
tref = 1 ;
DATA;
datosnew ;

% 1. Get Parameters
get_parameters_bubble

% 2. Setup parameters to estimate
setup_bubble

% 3. Wrap parameters
bubbledo

% 4. Compute posterior medians
% load chains from Metropolis Hasting
load lineardata
xx  = quantile(rr',0.5)' ;
% print table: 2 5percentile median 95percentile
disp("Table 2: Estimated Parameters -- Posterior")
table2 = (quantile(rr([5 6 7 8 2 4 1 3],1:100000)',[0.05 0.5 0.95]))';
disp(table2);
xlswrite("graph/table2.xlsx",table2)

% 5. Generate figures
% Figures 8 & 9
growth_max_plot(xx,par,zzd,tref,0) ;
% close all
% Figure 10
load data_start
growth_max(xx,par,zzd,tref,0) ;

% 6. extract data
hfig = open('bubbly_regime.fig');
fig_data = hfig.Children;
result_bubbly_regime= fig_data.Children.YData';
xlswrite("bubbly_regime.xlsx",result_bubbly_regime);
close gcf;

hfig = open('productivity.fig');
fig_data = hfig.Children;
result_productivity= fig_data.Children.YData';
close gcf;

hfig = open('preference.fig');
fig_data = hfig.Children;
result_preference= fig_data.Children.YData';
close gcf;

close all;