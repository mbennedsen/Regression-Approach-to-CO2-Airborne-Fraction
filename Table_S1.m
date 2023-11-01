%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Table S1 of the Supporting Information. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% (c) Mikkel Bennedsen (2023)
%
% This code can be used, distributed, and changed freely. Please cite Bennedsen,
% Hillebrand, and Koopman (2023): "A New Approach to the CO2 Airborne Fraction: Tackling Zero Emissions and Enhancing Statistical Precision".
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear; close all;
%% Init
filenam = 'AF_data.xlsx';

start_year = 1959;
end_year = 2021;

indx = 1:63; % 1:63 will select entire sample (1959-2021)

maxLags = 5;

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

%% Get DGP values
x_E = FF_GCP + LUC_GCP;
AF = y_ATM./x_E;

%% index
y_ATM = y_ATM(indx);
t = t(indx);
ENSO = ENSO(indx);
VAI = VAI(indx);
n = length(t);
x_E = x_E(indx);
AF = AF(indx);


%%

ADF_tab = nan(6,maxLags+1);
EG_tab = nan(2,maxLags+1);
%%  ADF

for i = 0:maxLags
    [h,pValue,stat,cValue] = adftest(y_ATM,'Lags',i,'Model','AR');
    ADF_tab(1,i+1) = pValue;

    [h,pValue,stat,cValue] = adftest(y_ATM,'Lags',i,'Model','ARD');
    ADF_tab(2,i+1) = pValue;

    [h,pValue,stat,cValue] = adftest(y_ATM,'Lags',i,'Model','TS');
    ADF_tab(3,i+1) = pValue;

    [h,pValue,stat,cValue] = adftest(x_E,'Lags',i,'Model','AR');
    ADF_tab(4,i+1) = pValue;

    [h,pValue,stat,cValue] = adftest(x_E,'Lags',i,'Model','ARD');
    ADF_tab(5,i+1) = pValue;

    [h,pValue,stat,cValue] = adftest(x_E,'Lags',i,'Model','TS');
    ADF_tab(6,i+1) = pValue;

    [h,pValue,stat,cValue] = egcitest([y_ATM,x_E],'Lags',i,'RReg','adf');
    EG_tab(1,i+1) = pValue;

    [h,pValue,stat,cValue] = egcitest([y_ATM,x_E],'Lags',i,'RReg','pp');
    EG_tab(2,i+1) = pValue; 
end


%%
res_tab = [ADF_tab;EG_tab];
disp(' ');
disp('Table S1:')
disp(res_tab)