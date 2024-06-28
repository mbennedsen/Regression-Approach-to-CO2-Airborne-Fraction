%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Figure 1 of the main paper. 
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


%% plot
fig1 = figure;
subplot(3,2,1)
plot(t,y_ATM,'k-','LineWidth',1.5), hold on
ylabel('GtC','Interpreter','latex','FontSize',8)
title('a) Atmospheric growth','FontSize',8,'Interpreter','latex');

lgd = legend('Data','Interpreter','latex','Location','NorthWest');
lgd.FontSize = 6;
legend('boxoff');

axis tight;
set(gca,'FontSize',8)
grid on


subplot(3,2,2)
plot(t,x_E,'k-','LineWidth',1.5), hold on
title('b) Emissions','FontSize',8,'Interpreter','latex');

lgd = legend('Data','Interpreter','latex','Location','NorthWest');
lgd.FontSize = 6;
legend('boxoff');

ylabel('GtC','Interpreter','latex','FontSize',8)
axis tight;
set(gca,'FontSize',8)
grid on


%% Ratio-based 
a1_hat = mean(AF); %
res1 = a1_hat-AF;

EstCov = hac(ones(n,1),AF,'display','off','intercept',false);
se = sqrt(EstCov(1,1));
CI1 = a1_hat*ones(n,1) - 1.96*se;
CI2 = a1_hat*ones(n,1) + 1.96*se;
t2 = [t;flipud(t)];

yy = AF;
XX = [ones(n,1),ENSO,VAI];
btmp = (XX'*XX)\XX'*yy;
a3_hat = btmp(1);
res3 = XX*btmp-AF;

EstCov = hac(XX,AF,'display','off','intercept',false);
SIG = EstCov;%sqrt(EstCov(1,1));
CI1_2 = nan(n,1);
CI2_2 = nan(n,1);
for i = 1:n
    CI1_2(i) = XX(i,:)*btmp - 1.96*sqrt(XX(i,:)*SIG*XX(i,:)');
    CI2_2(i) = XX(i,:)*btmp + 1.96*sqrt(XX(i,:)*SIG*XX(i,:)');
end
%t2 = [t;flipud(t)];


subplot(3,2,3)

plot(t,AF,'k-','LineWidth',1.5), hold on
plot(t,a1_hat*ones(length(t),1),'b-','LineWidth',1.5), hold on
plot(t,XX*btmp,'r-.','LineWidth',1.5), hold on
patch(t2,[CI1',flip(CI2')], 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none'), hold on
patch(t2,[CI1_2',flip(CI2_2')], 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none'), hold on
plot(t,AF,'k-','LineWidth',1.5), hold on
plot(t,a1_hat*ones(length(t),1),'b-','LineWidth',1.5), hold on
plot(t,XX*btmp,'r-.','LineWidth',1.5), hold on

ylabel('Atmospheric growth / Emissions','Interpreter','latex','FontSize',8)
title('c) Ratio-based estimators','FontSize',8,'Interpreter','latex');
axis([t(1),t(end),0.1,1.1]);

lgd = legend('Data',['Ratio-based (intercept = ',num2str(a1_hat,2),')'],['With ENSO+VAI (intercept = ',num2str(a3_hat,2),')'],'Interpreter','latex','Location','NorthWest');
lgd.FontSize = 6;
legend('boxoff');

set(gca,'FontSize',8)
grid on

subplot(3,2,4)
plot(t,res1,'b-','LineWidth',1.5), hold on
plot(t,res3,'r-.','LineWidth',1.5), hold on
title('d) Ratio-based residuals','FontSize',8,'Interpreter','latex');
axis([t(1),t(end),-0.35,0.35]);
set(gca,'FontSize',8)
grid on


%% Regression-based
yy = y_ATM;
XX2 = x_E;

[tmp,indx] = sort(x_E);
XX2 = XX2(indx);
yy = yy(indx);

btmp2 = (XX2'*XX2)\XX2'*yy;
a2_hat = btmp2(1);
res2 = x_E*btmp2 - y_ATM;

EstCov = hac(XX2,yy,'display','off','intercept',false);
se = sqrt(EstCov(1,1))*XX2;
CI1 = XX2*btmp2 - 1.96*se;
CI2 = XX2*btmp2 + 1.96*se;
t2 = [XX2(:,1);flipud(XX2(:,1))];

yy = y_ATM;
XX4 = [x_E,ENSO,VAI];

[tmp,indx] = sort(x_E);
XX4 = XX4(indx,:);
yy = yy(indx);

btmp4 = (XX4'*XX4)\XX4'*yy;
a4_hat = btmp4(1);
res4 = [x_E,ENSO,VAI]*btmp4 - y_ATM;

EstCov = hac(XX4,yy,'display','off','intercept',false);
SIG = EstCov;%sqrt(EstCov(1,1));
CI1_2 = nan(n,1);
CI2_2 = nan(n,1);
for i = 1:n
    CI1_2(i) = XX4(i,:)*btmp4 - 1.96*sqrt(XX4(i,:)*SIG*XX4(i,:)');
    CI2_2(i) = XX4(i,:)*btmp4 + 1.96*sqrt(XX4(i,:)*SIG*XX4(i,:)');
end


subplot(3,2,5)
scatter(x_E,y_ATM,'ko'), hold on
plot(XX2(:,1),XX2*btmp2,'b-','LineWidth',1.5), hold on
plot(XX4(:,1),XX4*btmp4,'r-.','LineWidth',1.5), hold on
patch(t2,[CI1',flip(CI2')], 'b', 'FaceAlpha', 0.5, 'EdgeColor', 'none'), hold on
patch(t2,[CI1_2',flip(CI2_2')], 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none'), hold on
scatter(x_E,y_ATM,'ko'), hold on
plot(XX2(:,1),XX2*btmp2,'b-','LineWidth',1.5), hold on
plot(XX4(:,1),XX4*btmp4,'r-.','LineWidth',1.5), hold on

xlabel('Emissions','Interpreter','latex','FontSize',8)
ylabel('Atmospheric growth','Interpreter','latex','FontSize',8)
title('e) Regression-based estimators','FontSize',8,'Interpreter','latex');
axis([min(x_E)*0.95,max(x_E)*1.05,0,9]);
%title('e) Atmospheric growth / Emissions','FontSize',8,'Interpreter','latex');
lgd = legend('Data',['Regression-based (slope = ',num2str(a2_hat,2),')'],['With ENSO+VAI (slope = ',num2str(a4_hat,2),')'],'Interpreter','latex','Location','NorthWest');
lgd.FontSize = 6;
legend('boxoff');

%axis tight;
set(gca,'FontSize',8)
grid on

subplot(3,2,6)
plot(t,res2,'b-','LineWidth',1.5), hold on
plot(t,res4,'r-.','LineWidth',1.5), hold on
title('f) Regression-based residuals','FontSize',8,'Interpreter','latex');
%axis tight;
axis([t(1),t(end),-2.5,2.5]);
set(gca,'FontSize',8)
grid on

%% JB test on residuals from regression-based estimator
[~,pval] = jbtest(res2);
disp(['p-value of JB test = ',num2str(pval)]);

