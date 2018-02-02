% Here we define the objective for our optimization
% We wish to move the endpoint of one of the links to a desired horizontal
% distance.

function [ cost ] = obj_ex3( pinput, param )
    


    % Mechanism Definition
    [nS, eS, param] = mech_def;
    nElem = size(eS,1);
    nNode = size(nS,1);
    B = param.B; H = param.H; Y = param.Y; inc = param.inc;  
    bc = param.bc;
    
    %Enter applied forces in the following format
    %f=[a b c];
    %a = node number
    %b=degree of freedom of the corresponding node. x=1,y=2,theta=3
    %c=value of the force
    f = [1 2 pinput(1)];
    
    %Run Nonlinear FEA
    [u, Ri,alpha] = DemoNonlinearCode(nS,eS,Y,H,B,f,bc,inc,'off');
    
    % A vector initialization for storing positions of deflected nodes
    d = []; 
    for j=1:inc
        for i=1:nElem,
            m=eS(i,1);
            n=eS(i,2);
            d = [d; ...
                 nS([m n],1)+[u(m,j) u(n,j)]' nS([m n],2)+[u(m+nNode,j) u(n+nNode,j)]'
                  ];    
            axis equal            
        end
    end    
    % Separating out the vertical position of end point of the fifth link
    d = d(200,2);         
    cost = (d - 25)^2;
end

