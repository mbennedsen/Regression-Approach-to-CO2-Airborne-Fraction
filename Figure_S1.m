%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Figure S1 of the Supporting Information. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (c) Mikkel Bennedsen (2024)
%
% This code can be used, distributed, and changed freely. Please cite Bennedsen,
% Hillebrand, and Koopman (2024): "A Regression-Based Approach to the CO2 Airborne Fraction: Enhancing Statistical Precision and Tackling Zero Emissions".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('Data');
%% Init
filenam = 'AF_data.xlsx';

start_year = 1959;
end_year = 2022;

%% Load data
dat = xlsread(filenam,1);

%% Construct data

%%% GHG %%%
N1 = sum(dat(:,1)<start_year)+1;
N2 = sum(dat(:,1)<end_year)+1;

t       = dat(N1:N2,1);
FF_GCP  = dat(N1:N2,4);
y_ATM   = dat(N1:N2,5);
LUC_GCP = dat(N1:N2,6);
LUC_HC  = dat(N1:N2,7);
LUC_NEW = dat(N1:N2,8);

ENSO = dat(N1:N2,10);
VAI = dat(N1:N2,9);

n = length(t);

x_E = FF_GCP + LUC_GCP;
AF = y_ATM./x_E;

%% Estimate b
b_hat = mean(diff(x_E));
b_std = sqrt(var(diff(x_E))/(length(diff(x_E))-1));

disp(['bhat = ',num2str(b_hat),' (SE = ',num2str(b_std),')']);

res = diff(x_E)-b_hat;
[~,pval] = jbtest(res);
disp(['p-value of JB test = ',num2str(pval)]);


%% plot
fig2 = figure;
subplot(1,2,1);
plot(t(2:end),diff(x_E),'b-','LineWidth',1.5), hold on
title('a) $\Delta E_t$','FontSize',8,'Interpreter','latex');
ylabel('GtC','FontSize',8)
set(gca,'FontSize',8)
grid on

subplot(1,2,2);
autocorr(diff(x_E)); 
title('b) Autocorrelation of $\Delta E_t$','FontSize',8,'Interpreter','latex');
set(gca,'FontSize',8)
grid on


