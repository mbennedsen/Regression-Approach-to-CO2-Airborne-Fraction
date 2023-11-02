%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Figure 2 of the main paper. 
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


%% plot
fig2 = figure;
subplot(1,2,1);
plot(t(2:end),diff(x_E),'b-','LineWidth',1.5), hold on
title('a) $\Delta E_t$','FontSize',8,'Interpreter','latex');
ylabel('GtC','FontSize',8)
%axis tight;
%lgd = legend('Atm. growth','Emissions','Interpreter','latex','Location','NorthWest');
%lgd.FontSize = 8;
%legend('boxoff');
set(gca,'FontSize',8)
grid on

subplot(1,2,2);
autocorr(diff(x_E));%,'b-','LineWidth',1.5), hold on
title('b) Autocorrelation of $\Delta E_t$','FontSize',8,'Interpreter','latex');
%axis tight;
%lgd = legend('Atm. growth','Emissions','Interpreter','latex','Location','NorthWest');
%lgd.FontSize = 8;
%legend('boxoff');
set(gca,'FontSize',8)
grid on

