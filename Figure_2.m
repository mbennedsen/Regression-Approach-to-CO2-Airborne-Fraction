%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This script will produce Figure 2 of the main paper. 
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

%% Load GCB data
dat = xlsread(filenam,1);

%%% Assign data values
y_FF = dat(:,4);
y_LUC = dat(:,6);
y_ATM = dat(:,5);

%%% Other data
conc_1750 =  278; % ppm
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
AF_noise = G_noise./E_noise;

%% State space modelling
y = [G_noise];
x = [E_noise];
t2 = t;

yesDiffuse = 0;
yes_refine = 0;
params0 = [0;0];

if sum(x<0) >0
    tmp = find(x<0);
    break_time = tmp(1);
end
Opts = optimset('Display','iter','TolFun',1e-12,'MaxFunEvals',1e6,'MaxIter',1e6);

    
%% Setup SS model
if yesDiffuse == 1
    if sum(x<0) >0
        Mdl = dssm(@(params)helpfct_TVP_regression_reflect_v01(params,length(x),x,break_time,yesDiffuse));
    else
        Mdl = dssm(@(params)helpfct_TVP_regression_v01(params,length(x),x,yesDiffuse));
    end
else
    if sum(x<0) >0
        Mdl = ssm(@(params) helpfct_TVP_regression_reflect_v01(params,length(x),x,break_time,yesDiffuse));
    else
        Mdl = ssm(@(params) helpfct_TVP_regression_v01(params,length(x),x,yesDiffuse));
    end
end


%% Estimate SS model
if yes_refine == 1
    %% Find good initial values with refine function
    refine_output = refine(Mdl,y,params0);
    logL = cell2mat({refine_output.LogLikelihood})';
    [~,maxLogLIndx] = max(logL);
    refinedParams0 = refine_output(maxLogLIndx).Parameters;
    
    [EstMdl,estParams0,EstParamCov,logL,output_est] = estimate(Mdl,y,refinedParams0,'Options',Opts); 
else
    [EstMdl,estParams0,EstParamCov,logL,output_est] = estimate(Mdl,y,params0,'Options',Opts);
end

%% Look at estimated parameters
sigu_hat = exp(estParams0(1));
sigu_std = sqrt(EstParamCov(1,1))*exp(estParams0(1)); % Delta rule

siga_hat = exp(estParams0(2));
siga_std = sqrt(EstParamCov(2,2))*exp(estParams0(2)); % Delta rule


ML_est = [sigu_hat;siga_hat];
ML_std = [sigu_std;siga_std];
ML_tstat = ML_est./ML_std;
%% Print parameters to screen
disp(' ');
disp('Parameters estimated by ML:    sig_u    sig_a');
disp(['Estimates:                :   ',num2str(ML_est',3)]);
disp(['(Std. Errs)               :   ',num2str(ML_std',3)]);
disp(['t-stats                   :   ',num2str(ML_tstat',3)]);



%% Get smoothed states
stVal = 1;
[x_smooth,logL_sm,Output_sm]   = smooth(EstMdl,y);

x_sm_cov = cell2mat({Output_sm.SmoothedStatesCov})'; % Covar matrix of filtered states.
x_sm_cov = [nan(size(x_smooth,2)*(stVal-1),size(x_smooth,2));x_sm_cov];

% Construct covar matrix of smoothed states
std_x_smooth = nan(size(x_smooth));
cov_smooth= nan(size(x_smooth,2),size(x_smooth,2),size(x_smooth,1));
for i = 1:size(x_smooth,1)
    cov_smooth(:,:,i) = x_sm_cov( ((i-1)*size(x_smooth,2)+1):i*size(x_smooth,2),: );
    std_x_smooth(i,:) = real( sqrt(diag(cov_smooth(:,:,i)))' );
end


%%
fig2 = figure;
subplot(3,1,1)
plot(t_GCB,y_ATM,'k-','LineWidth',1.5), hold on
plot(t,G_noise,'b-','LineWidth',1.5), hold on
plot(t,G_ATM,'m-','LineWidth',1.5), hold on
plot(t_GCB,y_ATM,'k-','LineWidth',1.5), hold on
plot(t,G_noise,'b-','LineWidth',1.5), hold on
plot(t,G_ATM,'m-','LineWidth',1.5), hold on
title('Atmospheric changes','FontSize',8,'Interpreter','latex');
ylabel('GtC/yr','FontSize',8,'Interpreter','latex');
xlabel('Year','FontSize',8,'Interpreter','latex');
if SSP == 6 || SSP == 4
    lgd = legend('Historical data',[plStr,' data (perturbed)'],[plStr,' data (original)'],'Interpreter','latex','Location','NorthWest');
else
    lgd = legend('Historical data',[plStr,' data (perturbed)'],[plStr,' data (original)'],'Interpreter','latex','Location','SouthWest');
end
text(0.02, 0.98, 'a)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
lgd.FontSize = 9;
legend('boxoff');
set(gca,'FontSize',9)
axis([1959,end_year,-5.5,11])
set(lgd, 'Units', 'normalized');
set(lgd, 'Position', [0.23 0.13+.6 0.05 0.05]); 


subplot(3,1,2)
plot(t_GCB,E_GCB,'k-','LineWidth',1.5), hold on
plot(t,E_noise,'b-','LineWidth',1.5), hold on
plot(t,E,'m-','LineWidth',1.5), hold on
plot(t_GCB,E_GCB,'k-','LineWidth',1.5), hold on
plot(t,E_noise,'b-','LineWidth',1.5), hold on
plot(t,E,'m-','LineWidth',1.5), hold on
title('Emissions','FontSize',8,'Interpreter','latex');
ylabel('GtC/yr','FontSize',8,'Interpreter','latex');
xlabel('Year','FontSize',8,'Interpreter','latex');
if SSP == 6 || SSP == 4
    lgd = legend('Historical data',[plStr,' data (perturbed)'],[plStr,' data (original)'],'Interpreter','latex','Location','NorthWest');
else
    lgd = legend('Historical data',[plStr,' data (perturbed)'],[plStr,' data (original)'],'Interpreter','latex','Location','SouthWest');
end
text(0.02, 0.98, 'b)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');    
lgd.FontSize = 9;
legend('boxoff');
set(gca,'FontSize',9)
axis([1959,end_year,-5.5,12])
set(lgd, 'Units', 'normalized');
set(lgd, 'Position', [0.23 0.13+.3 0.05 0.05]); 


subplot(3,1,3)

if sum(x<0)>0
    plot(t_GCB,y_ATM./E_GCB,'k-','LineWidth',1.5), hold on
    plot(t,AF_noise,'b-','LineWidth',1.5), hold on
    plot(t2(stVal:break_time-1),x_smooth(stVal:break_time-1,1),'r-.','LineWidth',2), hold on
    plot(t2(break_time:end),x_smooth(break_time:end,1),'r-.','LineWidth',2), hold on
    t3 = [t2(stVal:break_time-1);flipud(t2(stVal:break_time-1))];
    CI1 = x_smooth(stVal:break_time-1,1) - 1.96*std_x_smooth(stVal:break_time-1,1);
    CI2 = x_smooth(stVal:break_time-1,1) + 1.96*std_x_smooth(stVal:break_time-1,1);
    patch(t3,[CI1',flip(CI2')], 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none'), hold on
    
    t3 = [t2(break_time:end);flipud(t2(break_time:end))];
    CI1 = x_smooth(break_time:end,1) - 1.96*std_x_smooth(break_time:end,1);
    CI2 = x_smooth(break_time:end,1) + 1.96*std_x_smooth(break_time:end,1);
    patch(t3,[CI1',flip(CI2')], 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none'), hold on
    
    
    plot(t_GCB,y_ATM./E_GCB,'k-','LineWidth',1.5), hold on
    plot(t,AF_noise,'b-','LineWidth',1.5), hold on
    plot(t2(stVal:break_time-1),x_smooth(stVal:break_time-1,1),'r-.','LineWidth',2), hold on
    plot(t2(break_time:end),x_smooth(break_time:end,1),'r-.','LineWidth',2), hold on
else
    plot(t_GCB,y_ATM./E_GCB,'k-','LineWidth',1.5), hold on
    plot(t,AF_noise,'b-','LineWidth',1.5), hold on
    plot(t2(stVal:end),x_smooth(stVal:end),'r-.','LineWidth',2), hold on
    t3 = [t2;flipud(t2)];
    CI1 = x_smooth(stVal:end) - 1.96*std_x_smooth(stVal:end,1);
    CI2 = x_smooth(stVal:end) + 1.96*std_x_smooth(stVal:end,1);
    patch(t3,[CI1',flip(CI2')], 'r', 'FaceAlpha', 0.5, 'EdgeColor', 'none'), hold on
    
    plot(t_GCB,y_ATM./E_GCB,'k-','LineWidth',1.5), hold on
    plot(t,AF_noise,'b-','LineWidth',1.5), hold on
    plot(t2(stVal:end),x_smooth(stVal:end),'r-.','LineWidth',2), hold on

end

axis([1959,end_year,-1,2])
ylabel('Unitless','FontSize',8,'Interpreter','latex');
xlabel('Year','FontSize',8,'Interpreter','latex');
title('Atmospheric changes / Emissions','FontSize',8,'Interpreter','latex');
if SSP == 6
    lgd = legend('Historical ratio, $G_t/E_t$','Ratio-based estimate, $\widehat \alpha_{1,t} = G_t/E_t$','Regression-based estimate, $\widehat \alpha_{2,t}$','Interpreter','latex','Location','NorthWest');
else
    lgd = legend('Historical ratio, $G_t/E_t$','Ratio-based estimate, $\widehat \alpha_{1,t} = G_t/E_t$','Regression-based estimate, $\widehat \alpha_{2,t}$','Interpreter','latex','Location','SouthWest');
end
text(0.02, 0.98, 'c)', 'Units', 'normalized', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top');
lgd.FontSize = 9;
legend('boxoff');
set(gca,'FontSize',9)

set(lgd, 'Units', 'normalized');
set(lgd, 'Position', [0.26 0.13 0.05 0.05]);




