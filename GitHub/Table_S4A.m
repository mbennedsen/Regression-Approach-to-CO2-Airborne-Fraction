%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Table S4 (left panel) of the Supporting Information. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (c) Mikkel Bennedsen (2024)
%
% This code can be used, distributed, and changed freely. Please cite Bennedsen,
% Hillebrand, and Koopman (2024): "A Regression-Based Estimator of the CO2 Airborne Fraction: Enhancing Statistical Precision and Tackling Zero Emissions".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath('Data');
%% Init
filenam = 'AF_data.xlsx';

title_str = {'Data: GCP','Data: H&C','Data: vMa'};

start_year = 1959;
end_year = 2022;

indx = 1:64; % 1:64 will select entire sample (1959-2022)

detrend_ENSO = 1; % if = 1, then detrend ENSO before running analysis.

delta_grid = [0.2,0.5,1,2,5];
%% Load data
dat = xlsread(filenam,1);


%% Construct data
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
%% Detrend ENSO data
if detrend_ENSO == 1
    X = [ones(length(ENSO),1),t-t(1)];
    y = ENSO;
    bhat = (X'*X)\X'*y;

    ENSO = ENSO - X*bhat;
end


%% Index data to requested sub-sample
y_ATM = y_ATM(indx);
t = t(indx);
ENSO = ENSO(indx);
VAI = VAI(indx);

n = length(t);
%% Analyze data
for j = 1:3
    if j == 1 % Data: GCP (raw)
        disp('Analyzing: CGP AF data...');
        LUC = LUC_GCP;
    elseif j == 2 % Data: H&N (raw)
        disp(' ');
        disp(' ');
        disp('Analyzing: H&C AF data...');
        LUC = LUC_HN;
    elseif j == 3 % Data: New (raw)
        disp(' ');
        disp(' ');
        disp('Analyzing: vMa AF data...');
        LUC = LUC_NEW;
    end   
    x_E = FF_GCP + LUC;
    x_E = x_E(indx);
    AF = y_ATM./x_E;

    %% Deming estimator -- no intercept
    x = x_E;
    y = y_ATM;

    Axx = mean(x.^2);
    Ayy = mean(y.^2);
    Axy = mean(x.*y);
    beta = nan(length(delta_grid),1);
    for iD = 1:length(delta_grid)
        delta = delta_grid(iD);
    
        beta(iD) = (Ayy - delta*Axx + sqrt( (Ayy-delta*Axx)^2 + 4*delta*Axy^2))/2/Axy;
    end

    
    output = beta'; 
    
    %% Create table
    tab_res = [output];
    
    
    %% Print
    disp(' ')
    disp(title_str{j});
    disp(['    delta = ',num2str(delta_grid)])
    disp(round(tab_res,4));
    


end

