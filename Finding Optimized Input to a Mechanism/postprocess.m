function postprocess(presult,param)
    
    [nS, eS, param] = mech_def;
    B = param.B; H = param.H; Y = param.Y; inc = param.inc;  
    bc = param.bc;
    nNode = size(nS,1);
    nElem = size(eS,1);
    
    % Plotting the Mechanisms
    
    figure 
    clf
    subplot(121)
    % Plotting the undeflected mechanism
    for i=1:nElem
         m=eS(i,1);
         n=eS(i,2);
         plot(nS([m n],1),nS([m n],2),'Linewidth',2); hold on
    end
    title('The Undeflected Beam')
    
    % The optimal loading condition:
    sprintf('The optimal loading is: \n ')
    f=[1 2 presult]
    [u, Ri,alpha]=DemoNonlinearCode(nS,eS,Y,H,B,f,bc,inc,'off');
    
    % Plotting the deflected mechanism
   
    subplot(122)
    
    d = [];
   for j=1:20
        for i=1:nElem,
            m=eS(i,1);
            n=eS(i,2);
            plot(nS([m n],1),nS([m n],2),'Linewidth',2); hold on
            plot(nS([m n],1)+[u(m,j) u(n,j)]',nS([m n],2)+[u(m+nNode,j) u(n+nNode,j)]','r','Linewidth',2); hold on   
            d = [d; ...
                 nS([m n],1)+[u(m,j) u(n,j)]' nS([m n],2)+[u(m+nNode,j) u(n+nNode,j)]'
                  ];    

            axis equal
        end
   end
  title('The Deflected Beam')
  sprintf('The position of deflected internal link is: \n')
d (end,2)    
end