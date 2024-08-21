%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Panel A of Table 1 (left panel) of the main paper. 
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

indx = 1:64; % 1:64 will select entire sample (1959-2022)

detrend_ENSO = 1; % if = 1, then detrend ENSO before running analysis.

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
disp('Analyzing: GCP AF data...');
LUC = LUC_GCP;
x_E = FF_GCP + LUC;
x_E = x_E(indx);
AF = y_ATM./x_E;

%% Common estimator
a_hat0 = mean(AF);
s20 = var(AF);

% Use HAC estimator
EstCov = hac(ones(n,1),AF,'display','off','intercept',false);
se_HAC0 = sqrt(EstCov(1,1));

% 95% confidence interval (assuming Gaussianity)
CI00 = [a_hat0-1.96*se_HAC0;a_hat0+1.96*se_HAC0];

R2 = 1-sum((AF-a_hat0).^2)/sum( (AF-mean(AF)).^2);

output1 = [a_hat0;se_HAC0;se_HAC0/se_HAC0;CI00;sqrt(s20);R2];




%% New estimator
yy = y_ATM;
XX = x_E;
btmp = (XX'*XX)\XX'*yy;
s21 = sum( (yy-XX*btmp).^2 )/(n-length(btmp));

a_hat1 = btmp(1);

% w HAC
EstCov = hac(XX,yy,'display','off','intercept',false);
se_HAC1 = sqrt(EstCov(1,1));

% 95% confidence interval (valid by asymptotic arguments)
CI11 = [a_hat1-1.96*se_HAC1;a_hat1+1.96*se_HAC1];

R2 = 1-sum((yy-XX*btmp).^2)/sum( (yy-mean(yy)).^2);

output2 = [a_hat1;se_HAC1;se_HAC1/se_HAC0;CI11;sqrt(s21);R2];




%% Common estimator - with VAI and ENSO
yy = AF;
XX = [ones(n,1),ENSO,VAI];
btmp = (XX'*XX)\XX'*yy;
s20_EV = sum( (yy-XX*btmp).^2 )/(n-length(btmp));

a_hat0_EV = btmp(1);

% w HAC
EstCov = hac(XX,yy,'display','off','intercept',false);
se_HAC0_EV = sqrt(EstCov(1,1));

% 95% confidence interval (assuming Gaussianity)
CI00_EV = [a_hat0_EV-1.96*se_HAC0_EV;a_hat0_EV+1.96*se_HAC0_EV];

R2 = 1-sum((yy-XX*btmp).^2)/sum( (yy-mean(yy)).^2);

output3 = [a_hat0_EV;se_HAC0_EV;se_HAC0_EV/se_HAC0;CI00_EV;sqrt(s20_EV);R2];




%% New estimator - with VAI and ENSO
yy = y_ATM;
XX = [x_E,ENSO,VAI];
btmp = (XX'*XX)\XX'*yy;
s21_EV = sum( (yy-XX*btmp).^2 )/(n-length(btmp));

a_hat1_EV = btmp(1);

% w HAC
EstCov = hac(XX,yy,'display','off','intercept',false);
se_HAC1_EV = sqrt(EstCov(1,1));

% 95% confidence interval (valid by asymptotic arguments)
CI11_EV = [a_hat1_EV-1.96*se_HAC1_EV;a_hat1_EV+1.96*se_HAC1_EV];

R2 = 1-sum((yy-XX*btmp).^2)/sum( (yy-mean(yy)).^2);


output4 = [a_hat1_EV;se_HAC1_EV;se_HAC1_EV/se_HAC0;CI11_EV;sqrt(s21_EV);R2];


%% Create table
tab_res = [output1,output2,output3,output4];


%% Print
disp(' ')
disp('Table 1A (full sample, 1959-2022): ')
disp('    Eq. (1)   Eq. (2)   Eq. (3)   Eq. (4)')
disp(round(tab_res,4));




