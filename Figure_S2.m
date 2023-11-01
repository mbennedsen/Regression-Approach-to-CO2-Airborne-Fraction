%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Figure S2 of the Supporting Information. 
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
M = 1e6;

start_year = 1959; % For retrieving DGP values
end_year = 2021;

Tend   = 2100-1959+1; 
Tstart = 63;

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

%% Get DGP values
x_E = FF_GCP + LUC_GCP;
AF = y_ATM./x_E;

%%% Conventional estimator
a1 = mean(AF);
s21 = var(AF);
    
%%% OLS estimator
yy = y_ATM;
XX = x_E;
a2 = (XX'*XX)\XX'*yy;
s22 = sum( (yy-XX*a2).^2 )/(n-length(a2));

%%% DGP for emissions
yy = x_E;
XX = [ones(n,1),(1:n)'];
b_E = (XX'*XX)\XX'*yy;
s2E = sum( (yy-XX*b_E).^2 )/(n-length(b_E));

%% Run MC study
a1_sim = nan(M,Tend-Tstart+1);
a2_sim = nan(M,Tend-Tstart+1);
for i = 1:M
    %%% Simulate whole trajectory
    eps_t = randn(Tend,1); % Common random numbers
    E_sim = [ones(Tend,1),(1:Tend)']*b_E + sqrt(s2E)*randn(Tend,1);

    AF_sim  = a1 + sqrt(s21)*eps_t;
    G_sim   = a1*E_sim + sqrt(s22)*eps_t;


    for t = Tstart:Tend
        %%% Conventional estimator
        yy_sim = AF_sim(1:t);
        a1_sim(i,t-Tstart+1) = mean(yy_sim);

        %%% OLS estimator
        yy_sim = G_sim(1:t);
        XX_sim = E_sim(1:t);
        a2_sim(i,t-Tstart+1) = (XX_sim'*XX_sim)\XX_sim'*yy_sim;
    end

end


%% Calculate RMSE
a1_rmse = sqrt( mean( (a1_sim - a1).^2 ) );
a2_rmse = sqrt( mean( (a2_sim - a1).^2 ) );

%% plot
fig1 = figure(1);
subplot(1,2,1);
plot(Tstart:Tend,a1_rmse,'b-','LineWidth',1.5), hold on
plot(Tstart:Tend,a2_rmse,'r-.','LineWidth',1.5), hold on
lgd = legend('RMSE($\alpha_1$)','RMSE($\alpha_2$)','Interpreter','latex','Location','NorthEast');
lgd.FontSize = 8;
legend('boxoff');
ylabel('Root mean squared error (RMSE)','FontSize',8);
xlabel('T','FontSize',8);
axis tight;
set(gca,'FontSize',8)

subplot(1,2,2);
plot(Tstart:Tend,a2_rmse./a1_rmse,'b-','LineWidth',1.5), hold on
plot(Tstart:Tend,ones(Tend-Tstart+1,1),'k--'), hold on
lgd = legend('RMSE($\alpha_2$)/RMSE($\alpha_1$)','Interpreter','latex','Location','NorthEast');
lgd.FontSize = 8;
legend('boxoff');
ylabel('Relative RMSE','FontSize',8);
xlabel('T','FontSize',8);
axis tight;
set(gca,'FontSize',8)

