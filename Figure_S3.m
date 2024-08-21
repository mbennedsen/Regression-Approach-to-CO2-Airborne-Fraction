%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Figure S3 of the Supporting Information. 
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
addpath(genpath('Functions/'));
addpath('Data');
%% Init
filenam = 'AF_data.xlsx';

rng(666);
SSP = 2;

start_year = 2023;
end_year   = 2100;

tau_window = [50,20,10,5,1]; % Window-size

%% Load GCB data
dat = xlsread(filenam,1);

%%% Assign data values
y_FF = dat(:,4);
y_LUC = dat(:,6);
y_ATM = dat(:,5);

%%% Other data
conc_1750 =  278; % ppm
%conc_1850 =  284.7; % ppm
conc_1959 = 315.39; %ppm
C0  = conc_1959*2.127;
C00 = conc_1750*2.127;

tmp = y_ATM; tmp(isnan(tmp)) = 0;
y_C = C0 + cumsum(tmp);
y_C(isnan(y_ATM)) = nan;
E_GCB = y_FF+y_LUC;

t_GCB = dat(:,1);

sig_af = std(y_ATM./E_GCB);
sig_e  = std(diff(E_GCB));
sig_a  = std(diff(y_ATM));

%% Load RCP Data
if SSP == 1
    load('SSP119_output_v2');
    disp('Analysing RCP119 data...')

    plStr = 'SSP119';
elseif SSP==2
    load('SSP126_output_v2');
    disp('Analysing RCP126 data...')

    plStr = 'SSP126';
elseif SSP==3
    load('SSP245_output_v2');
    disp('Analysing RCP245 data...')

    plStr = 'SSP245';
elseif SSP==4
    load('SSP370_output_v2');
    disp('Analysing RCP370 data...')

    plStr = 'SSP370';
elseif SSP==5
    load('SSP434_output_v2');
    disp('Analysing RCP434 data...')

    plStr = 'SSP434';
elseif SSP==6
    load('SSP585_output_v2');
    disp('Analysing RCP585 data...')

    plStr = 'SSP585';
else
    error('RCP scenario not implemented...');
end


%% PLot AF
G_ATM = diff(C_magicc);
E = E_magicc;
t  = t_m;
te = t_e;

G_ATM(t(2:end)<start_year) = [];
E(te<start_year) = [];
t(t<start_year) = [];
te(te<start_year) = [];

G_ATM(t>end_year) = [];
E(te>end_year) = [];
t(t>end_year) = [];
te(te>end_year) = [];

AF = G_ATM./E;

E_noise = E + (sig_e*randn(length(t),1));
G_noise = G_ATM + (sig_a*randn(length(t),1));
AF_noise2 = AF + sig_af*randn(length(t),1);
AF_noise = G_noise./E_noise;
AF_noise3 = cumsum(G_noise)./cumsum(E_noise);

%% Make CAFs
fig1 = figure;
CAF_noise_tot = cumsum(G_noise)./cumsum(E_noise);

subplot(3,2,1);
plot(t,CAF_noise_tot,'b-','LineWidth',1.5), hold on
text(0.02, 0.98, 'a)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
set(gca,'FontSize',9)

ylabel('Unitless','FontSize',8,'Interpreter','latex');
xlabel('Year','FontSize',8,'Interpreter','latex');
title(['Cumulative airborne fraction ($w = ',num2str(length(t)),'$)'],'FontSize',8,'Interpreter','latex');
axis([2022,end_year,-1,2])


for j = 1:5
    tau0 = tau_window(j);

    CAF_noise = nan(length(t),1);
    for i = 1:length(t)
        if i < tau0
            CAF_noise(i) = sum(G_noise(1:i))./sum(E_noise(1:i));
        else
            CAF_noise(i) = sum(G_noise((i-tau0+1):i))./sum(E_noise((i-tau0+1):i));
        end
    end

    subplot(3,2,j+1);
    plot(t,CAF_noise,'b-','LineWidth',1.5), hold on
    

    axis([2022,end_year,-1,2])
    if j ==1
        text(0.02, 0.98, 'b)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    elseif j == 2
        text(0.02, 0.98, 'c)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    elseif j == 3
        text(0.02, 0.98, 'd)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    elseif j == 4
        text(0.02, 0.98, 'e)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    elseif j == 5
        text(0.02, 0.98, 'f)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
    end
    set(gca,'FontSize',9)

    ylabel('Unitless','FontSize',8,'Interpreter','latex');
    xlabel('Year','FontSize',8,'Interpreter','latex');
    title(['Cumulative airborne fraction ($w = ',num2str(tau0),'$)'],'FontSize',8,'Interpreter','latex');  
end


