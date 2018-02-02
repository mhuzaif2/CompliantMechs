# CompliantMechs
I have made this repository to share code snippets that I am making while taking SE 598 Compliant Systems Design at University of Illinois at Urbana-Champaign. 


1 - Finding Optimized Input to a Mechanism
------------------------------------------

A general script that can find loading (force) for a compliant mechanism using an FEA solver.
Your mechanism should be defined inside 'mech_def.m' using following way:

Node Set:
 nodeSet = [x y; ...]

Element Set: (a graph representing connections of nodes)
 elemSet = [i j; ...]

Boundary Conditions:
 bc=[a b c; ...];
 a = node number
 b=degree of freedom of the corresponding node. x=1,y=2,theta=3
 c=value of the displacement

Material Properties:
% Out-of-plane thickness of all elements
 B = ...

% In-plane thickness of all elements
 H = ...

% Young's Modulus for all elements
 Y = ...

Solver parameter:
 inc = number of increments for nonlinear loading computation



Pre-Requisite:

- Optimization Toolbox
- MATLAB R2006 and above


Notes:
DemoNonlinearCode.m shared in the SE 598 Compliant Systems Design by Prof. Girish Krishnan at University of Illinois at Urbana-Champaign.
