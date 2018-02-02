function [nodeSet, elemSet, param] = mech_def
    
    
    % Load the mechanism parameters
    % load('mech_params.mat');    
    % OR enter the parameters here
    
    % Enter the Element Connectivity
    elemSet=[1 2;
        2 3;
        3 4;    
        4 1;
        4 5];

    %Enter the nodal co-ordinates
    nodeSet=[0 0;
        10 20;
        0 40;    
        -10 20;
        0 20];

    %Enter the out-of-plane thickness for all the elements
    B=5*ones(length(elemSet(:,1)),1);

    %Enter the in-plane thickness for all the elements
    H=5*ones(length(elemSet(:,1)),1);

    
    %Number of load increments
    inc = 20;

    %Enter the Young's Modulus for all elements
    Y=200e2*ones(length(elemSet(:,1)),1);
    %%Enter applied boundary conditions in the following format
    %f=[a b c];
    %a = node number
    %b=degree of freedom of the corresponding node. x=1,y=2,theta=3
    %c=value of the displacement 

    bc=[1 1 0;
        1 3 0;
        2 1 0;
        3 2 0;
        3 3 0];
    %     5 1 0];    
    %     3 3 0];
    %     3 3 0];
    %     3 2 0];
    
    % Saving the parameters in a variable accessible in other snippets
    param.B = B; param.H = H; param.Y = Y; param.bc = bc; param.inc = inc;
    
end
