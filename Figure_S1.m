%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Figure S1 of the Supporting Information. 
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

indx = 1:63; % 1:63 will select entire sample (1959-2021)

detrend_ENSO = 1; % if = 1, then detrend ENSO.

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
FF_GCP = FF_GCP(indx);
LUC_GCP = LUC_GCP(indx);
LUC_HN = LUC_HN(indx);
LUC_NEW = LUC_NEW(indx);

n = length(t);

%% plot
fig1 = figure(1);
subplot(2,2,1);
plot(t,y_ATM,'b-','LineWidth',1.5), hold on
plot(1992*ones(100,1),linspace(0,6.2,100),'k--','LineWidth',1), hold on
title('a) Atmospheric CO2 changes ($G_t$)','FontSize',8,'Interpreter','latex');
ylabel('GtC/yr','FontSize',8,'Interpreter','latex');
grid on
set(gca,'FontSize',8)
axis tight;

subplot(2,2,2);
plot(t,FF_GCP,'b-','LineWidth',1.5), hold on
plot(t,LUC_GCP,'r-','LineWidth',1.5), hold on
plot(t,LUC_HN,'g-','LineWidth',1.5), hold on
plot(t,LUC_NEW,'c-','LineWidth',1.5), hold on
plot(1992*ones(100,1),linspace(0,10.5,100),'k--','LineWidth',1), hold on
title('b) CO2 emissions ($E_t$)','FontSize',8,'Interpreter','latex');
lgd = legend('FF (GCP)','LULCC (GCP)','LULCC (H\&N)','LULCC (vMa)','Interpreter','latex','Location','NorthWest');
lgd.FontSize = 6;
legend('boxoff');
ylabel('GtC/yr','FontSize',8,'Interpreter','latex');
grid on
axis tight;
set(gca,'FontSize',8)

subplot(2,2,3);
plot(t,ENSO,'b-','LineWidth',1.5), hold on
plot(1992*ones(100,1),linspace(-1.1,2.2,100),'k--','LineWidth',1), hold on
title('c) ENSO','FontSize',8,'Interpreter','latex');
grid on
ylabel('Dimensionless','FontSize',8,'Interpreter','latex');
axis tight;
set(gca,'FontSize',8)

subplot(2,2,4);
plot(t,VAI,'b-','LineWidth',1.5), hold on
plot(1992*ones(100,1),linspace(0,0.175,100),'k--','LineWidth',1), hold on
title('d) VAI','FontSize',8,'Interpreter','latex');
grid on
ylabel('Dimensionless','FontSize',8,'Interpreter','latex');
axis tight;
set(gca,'FontSize',8)
