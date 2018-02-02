# CompliantMechs
I have made this repository to share code snippets that I am making while taking SE 598 Compliant Systems Design at University of Illinois at Urbana-Champaign.


1 - Finding Optimized Input to a Mechanism
------------------------------------------

A general script that can find loading (force) for a compliant mechanism using an FEA solver. Currently it is finding input vertical loading force
that gives a desired vertical displacement by the internal element of a 5-bar linkage mechanism.

If you wish to change it to your own problem, your mechanism should be defined inside 'mech_def.m' using following way:

**Node Set:**
 nodeSet = [x y; ...]

**Element Set:** (a graph representing connections of nodes)

 elemSet = [i j; ...]

**Boundary Conditions:**

 bc=[a b c; ...];
 a = node number
 b=degree of freedom of the corresponding node. x=1,y=2,theta=3
 c=value of the displacement

**Material Properties:**

% Out-of-plane thickness of all elements
 B = ...

% In-plane thickness of all elements
 H = ...

% Young's Modulus for all elements
 Y = ...

**Solver parameter:**

 inc = number of increments for nonlinear loading computation

Furthermore, when specifying desired displacements in obj.m, understand the structure of 'd' vector before writing the cost function 'cost'. <Will add details to that later>

**Running:**

Just play main_optim.m and wait for the computation to finish.


**Pre-Requisite:**

- Optimization Toolbox
- MATLAB R2006 and above


**Notes:**

DemoNonlinearCode.m shared in the SE 598 Compliant Systems Design by Prof. Girish Krishnan at University of Illinois at Urbana-Champaign.
