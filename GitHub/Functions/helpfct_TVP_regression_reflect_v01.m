function [A,B,C,D,Mean0,Cov0,StateType] = helpfct_TVP_regression_reflect_v01(params,T,E,break_time,yesDiffuse)
%
% _reflect: Go to 1-alpha in break year
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
StateType = [2;1]; % 0: Stationary; 1: unity; 2: non-stationary.

if yesDiffuse == 1
    Mean0 = [0.44;1];%[0.44]; % Mean of initial states
    Cov0  = blkdiag(Inf,0);%blkdiag(Inf*eye(2),0); % Cov. matrix of initial states
else
    %Mean0 = [0;0;1]; % Mean of initial states
    %Cov0  = blkdiag(1e6*eye(2),0); % Cov. matrix of initial states

    Mean0 = [0.44;1];%[0;0;1]; % Mean of initial states
    Cov0  = blkdiag(1e3,0);%blkdiag(VAR_dT,Var_e,0); % Cov. matrix of initial states
end
%% Construct transition matrices

Atmp1 = [1,0;
         0,1];

Atmp2 = [-1,1;
         0,1];
%% Matrix C in obs. eq. 
A=[];
C=[];
for i = 1:T
    C = [C;{[E(i),0]}]; % Time-varying regressor

    if i == break_time
        A = [A;{Atmp2}];
    else
        A = [A;{Atmp1}];
    end
end

%% Construct error structure in state eq.
% Note: Matlab wants the (lower diagonal) sqrt-matrix of the var-covar matrix
R = [1;0];

B = R*siga;

%% Construct error structure in obs eq.
% Note: Matlab wants the (lower diagonal) sqrt-matrix of the var-covar matrix
D = sigu;

%% The matrices not defined to be time-varying must be "forced" to be time-varying (Matlab idiosyncracity) 
%A = repmat({A},T,1);
B = repmat({B},T,1);
     %repmat({B2},T2,1)];
%C = repmat({C},T,1);
D = repmat({D},T,1);

