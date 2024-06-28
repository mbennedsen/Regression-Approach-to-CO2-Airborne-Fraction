function [A,B,C,D,Mean0,Cov0,StateType] = helpfct_TVP_regression_v01(params,T,E,yesDiffuse)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Matlab ss representation:
%%%
%%% x_{t+1} = A*x_t + B*eta_t    (eta_t std. normal)
%%% y_t     = C*x_t + D*eps_t    (eps_t std. normal)
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load in parameters
sigu  = exp(params(1));
siga  = exp(params(2));

%% Set mean and covariance for initial states. 
StateType = [2]; % 0: Stationary; 1: unity; 2: non-stationary.

if yesDiffuse == 1
    Mean0 = 0.44;%[0.44]; % Mean of initial states
    Cov0  = Inf;%blkdiag(Inf*eye(2),0); % Cov. matrix of initial states
else
    %Mean0 = [0;0;1]; % Mean of initial states
    %Cov0  = blkdiag(1e6*eye(2),0); % Cov. matrix of initial states

    Mean0 = 0.44;%[0;0;1]; % Mean of initial states
    Cov0  = 1e3;%blkdiag(VAR_dT,Var_e,0); % Cov. matrix of initial states
end
%% Construct transition matrices

A = 1;

%% Matrix C in obs. eq. 
C=[];
for i = 1:T
    C = [C;{E(i)}]; % Time-varying regressor
end

%% Construct error structure in state eq.
% Note: Matlab wants the (lower diagonal) sqrt-matrix of the var-covar matrix
R = 1;

B = R*siga;

%% Construct error structure in obs eq.
% Note: Matlab wants the (lower diagonal) sqrt-matrix of the var-covar matrix
D = sigu;

%% The matrices not defined to be time-varying must be "forced" to be time-varying (Matlab idiosyncracity) 
A = repmat({A},T,1);
B = repmat({B},T,1);
     %repmat({B2},T2,1)];
%C = repmat({C},T,1);
D = repmat({D},T,1);

