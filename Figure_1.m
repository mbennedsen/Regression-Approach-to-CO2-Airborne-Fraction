%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Figure 1 of the main paper. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (c) Mikkel Bennedsen (2023)
%
% This code can be used, distributed, and changed freely. Please cite Bennedsen,
% Hillebrand, and Koopman (2023): "A New Approach to the CO2 Airborne Fraction: Enhancing Statistical Precision and Tackling Zero Emissions".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
%% Init
filenam = 'AF_data.xlsx';

start_year = 1959;
end_year = 2021;

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
LUC_HN  = dat(N1:N2,7);
LUC_NEW = dat(N1:N2,8);

ENSO = dat(N1:N2,10);
VAI = dat(N1:N2,9);

n = length(t);

x_E = FF_GCP + LUC_GCP;
AF = y_ATM./x_E;

%% Estimate a and get residuals
a_mean = mean(AF);
res_mean = AF-a_mean;
    
%%% OLS estimator
yy = y_ATM;
XX = x_E;
a_ols = (XX'*XX)\XX'*yy;
res_ols = yy-a_ols*XX;

%% plot
fig1 = figure;
subplot(2,2,1)
plot(t,y_ATM,'b-','LineWidth',1.5), hold on
ylabel('GtC','FontSize',8)
title('a) Atmospheric growth','FontSize',8,'Interpreter','latex');
axis tight;
set(gca,'FontSize',8)
grid on


subplot(2,2,2)
plot(t,x_E,'b-','LineWidth',1.5), hold on
title('b) Emissions','FontSize',8,'Interpreter','latex');
ylabel('GtC','FontSize',8)
axis tight;
set(gca,'FontSize',8)
grid on

subplot(2,2,3)
plot(t,AF,'b-','LineWidth',1.5), hold on
title('c) Atmospheric growth / Emissions','FontSize',8,'Interpreter','latex');
axis tight;
set(gca,'FontSize',8)
grid on

subplot(2,2,4)
plot(t,res_ols,'b-','LineWidth',1.5), hold on
title('d) Residuals','FontSize',8,'Interpreter','latex');
axis tight;
set(gca,'FontSize',8)
grid on

