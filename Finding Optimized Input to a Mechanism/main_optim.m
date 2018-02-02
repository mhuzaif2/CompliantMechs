
% An optimization formulation that calculates the required force input to a
% mechanism. The formulation uses the 'fmincon' function in MATLAB. 

% Note: No nonlinear constraints are given but can be added as a separate
% script.

% Umer Huzaifa
% 1st Feb 2018

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Mechanism Definitions

[nS, nE, param] = mech_def;
f0 = 100;  % An initial applied load
pinput = f0; % optimization variable

   
% Give linear equality and inequality constraints
Aineq = []; Bineq = []; % we have no linear inequality constraints
Aeq = []; Beq = []; % we have no equality constraints
LB = 0; UB = Inf; 

options = optimset('display','iter','diffmaxchange',1.1*1e-5, ...
    'diffminchange',1e-5,'MaxFunEvals',20000);

% Find the input force achieving the objective
[presult,optfval] = ...
    fmincon(@obj,pinput,Aineq,Bineq,Aeq,Beq,LB,UB,[],options,param);

% Compute the curvature and plot the deflected mechanism
postprocess(presult,param);

