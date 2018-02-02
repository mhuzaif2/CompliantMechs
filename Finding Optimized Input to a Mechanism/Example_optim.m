
% An optimization formulation that calculates the required force input to a
% 5 bar linkage mechanism giving a required displacement.

% Umer Huzaifa
% 1st Feb 2018


% Mechanism Definitions

f0 = 4500;  % An initial applied load
pinput = f0; % optimization variable

   
Aineq = []; Bineq = []; % we have no linear inequality constraints
Aeq = []; Beq = []; % we have no equality constraints
LB = 0; UB = Inf; 

options = optimset('display','iter','diffmaxchange',1.1*1e-5, ...
    'diffminchange',1e-5,'MaxFunEvals',20000);

% Find the input force achieving the objective
[presult,optfval] = ...
    fmincon(@obj_ex3,pinput,Aineq,Bineq,Aeq,Beq,LB,UB,[],options,param);

% Compute the curvature and plot the deflected mechanism
postprocess_ex3(presult,param);

