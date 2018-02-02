function [u, Ri,alpha]=DemoNonlinearCode(nodeSet,elemSet,Y,H,B,F,bc,numSubstep,cmdout)
%f=zeros(3*nNode,1)
%bc=dof that is constrained
switch cmdout
    case 'on'
        header = sprintf('\n\t Substep\t\tIter#\t\t  ConvParm\t\t\t  MaxDisp');
        formatstr = '\t  %3d\t\t    %3d\t\t  %16.4f\t\t  %6.3f';
        disp(header);
end
nDOF = 3;
nNode=length(nodeSet(:,1));
totalDOF=nDOF*nNode;
nElem = length(elemSet(:,1));

BC=nan(nNode*3,1);
for i=1:length(bc(:,1)),
    BC(bc(i,1)+nNode*(bc(i,2)-1))=bc(i,3);
end

f=zeros(nNode*3,1);
for i=1:length(F(:,1))
    f(F(i,1)+nNode*(F(i,2)-1))=F(i,3);
end

Kglobal = getStiffMatrix(nodeSet,elemSet,Y,H,B);
K = Kglobal;
bcwt=mean(diag(K));
freeBC=isnan(BC);
cnstBC=~freeBC;
f = externalLoadReduction(f,K,BC,freeBC,cnstBC,bcwt);
K = StiffnessReduction(K,cnstBC,BC,bcwt);
uLinear=K\f;
Freaction=Kglobal*uLinear;
F=Freaction;
u=zeros(totalDOF,numSubstep); %Displacement
Ri=zeros(totalDOF,numSubstep); %Internal Force,Ri
lambda=zeros(numSubstep,1);
K=zeros(totalDOF,totalDOF,numSubstep); %Tangential Stiffness
R=zeros(totalDOF,numSubstep); %Residual Force,R
alpha=zeros(nElem,1); %Rigid body rotation for each element

for substep=1:numSubstep
    lambda(substep)=substep/numSubstep;
    iter=0;
    if substep==1 && iter==0
        K(:,:,substep)=full(Kglobal);
        f=F*lambda(substep);
        bcwt=mean(diag(K(:,:,substep)));
        fred = externalLoadReduction(f,K(:,:,substep),BC,freeBC,cnstBC,bcwt);
        Kred = StiffnessReduction(K(:,:,substep),cnstBC,BC,bcwt);
        u(:,substep)=Kred\fred;
    else
        f=F*lambda(substep);
        u(:,substep)=u(:,substep-1);
    end
    convParam=1;
    while convParam>1e-8 && iter<10000
        iter=iter+1;
        [K(:,:,substep) ,Ri(:,substep),alpha] =...
            getGTStiffness(elemSet,nodeSet,Y,H,B,u(:,substep),alpha);
        R(:,substep)=f-Ri(:,substep);
        convParam=norm(R(freeBC,substep))^2/(1+norm(f(freeBC))^2);
        bcwt=mean(diag(K(:,:,substep)));
        Rred = externalLoadReduction(R(:,substep),K(:,:,substep),BC,freeBC,cnstBC,bcwt);
        Kred = StiffnessReduction(K(:,:,substep),cnstBC,BC,bcwt);
        u(:,substep)=u(:,substep)+Kred\Rred;
        
        switch cmdout
            case 'on'
                MaxDisp=max(abs(u(:,substep)));
                currOutput = sprintf(formatstr,substep,iter,convParam,MaxDisp);
                disp(currOutput);
        end
    end
    switch cmdout
        case 'on'
            disp(sprintf(' '));
    end
end
end

function K = StiffnessReduction(K,cnstBC,BC,bcwt)
K(cnstBC,:)=0;
K(:,cnstBC)=0;
K(cnstBC,cnstBC)=bcwt*speye(size(BC(cnstBC),1));%speye(length(BC(cnstBC)));
end

%--------------------------------------------------------------------------
function f = externalLoadReduction(f,K,BC,freeBC,cnstBC,bcwt)
f(freeBC) = f(freeBC)-(K(freeBC,cnstBC)*BC(cnstBC));
f(cnstBC) = bcwt*BC(cnstBC);
end

function GM = getStiffMatrix(nodeSet,elemSet,Y,H,B)

%elemSet = elemSet(:);nodeSet = nodeSet(:);
nElem=size(elemSet,1);nNode=size(nodeSet,1);
nDOF=3;
GM = sparse(nNode*nDOF,nNode*nDOF);
elemCnnt = elemSet;
for i=1:nElem
    sctr = elemCnnt(i,:);
    %                 sctrVct=[];
    sctrVct = [sctr, sctr+nNode, sctr+2*nNode];
    %                 for indexDOF=1:nDOF
    %                     sctrVct = cat(2,sctrVct, sctr+(indexDOF-1)*nNode);
    %                 end
    
    EMatrix=getEMatrixStiff(elemSet(i,:),nodeSet,Y,H,B);
    
    EMatrixTrans=getEMatrixTrans(elemSet(i,:),nodeSet);
    GM(sctrVct,sctrVct) = GM(sctrVct,sctrVct) +...
        EMatrixTrans'*EMatrix*EMatrixTrans;
end
end

function Ke = getEMatrixStiff(elemSet,nodeSet,Y,H,B)

%           obj=obj(:);nodeSet=nodeSet(:);matSet=matSet(:);sectSet=sectSet(:);
nElem=length(elemSet(:,1)) ;
Ke=zeros(6,6,nElem);
for i=1:nElem
    
    in=elemSet(i,1);
    jn=elemSet(i,2);
    E=Y(i);
    A=B(i)*H(i);
    I=B(i)*H(i)^3/12;
    L=norm(nodeSet(jn,:)-nodeSet(in,:));
    L2=L^2;
    IdL=6*I/L;
    IdL2=12*I/L2;
    I_4=4*I;
    I_2=2*I;
    Ke(:,:,i) = E./L.*[...
        A   -A    0       0       0      0   ;...
        -A   A    0       0       0      0   ;...
        0    0    IdL2   -IdL2    IdL    IdL ;...
        0    0   -IdL2    IdL2   -IdL   -IdL ;...
        0    0    IdL    -IdL     I_4    I_2 ;...
        0    0    IdL    -IdL     I_2    I_4 ];
end
end

function Te = getEMatrixTrans(elemSet,nodeSet)

nElem=length(elemSet(:,1)) ;
Te=zeros(6,6,nElem);
for i=1:nElem
    %                 in=find(nodeNumbers==obj(i,1).connectivity(1), 1);
    %                 jn=find(nodeNumbers==obj(i,1).connectivity(2), 1);
    in=elemSet(i,1);
    jn=elemSet(i,2);
    theta=atan2(nodeSet(jn,2)-nodeSet(in,2),...
        nodeSet(jn,1)-nodeSet(in,1));
    Ecos=cos(theta);Esin=sin(theta);
    Te(:,:,i) = [...
        Ecos  0     Esin  0     0   0;...
        0     Ecos  0     Esin  0   0;...
        -Esin  0     Ecos  0     0   0;...
        0    -Esin  0     Ecos  0   0;...
        0     0     0     0     1   0;...
        0     0     0     0     0   1];
end
end

function [K, f, alpha] = getGTStiffness(elemSet,nodeSet,Y,H,B,u,alpha)
nDOF = 3;
nNode = length(nodeSet(:,1));
nElem = length(elemSet(:,1));
elemCnnt=elemSet;
totDOF=nNode*nDOF;
K = sparse(totDOF,totDOF);
f=zeros(totDOF,1);
nodeLoc=nodeSet;
for i=1:nElem
    sctr = elemCnnt(i,:); %[nodeI, nodeJ] for elemIndex
    sctrVct = [sctr,sctr+nNode,sctr+nNode*2];
    Pg=u(sctrVct);
    [Kg Fg alpha(i)] =...
        getEMatrixTStiff(elemSet(i,:),nodeSet,Y(i),H(i),B(i),Pg,alpha(i));
    K(sctrVct,sctrVct) = K(sctrVct,sctrVct) + Kg;
    f(sctrVct,1)=f(sctrVct,1)+Fg;
end
end

function [Kg Fg alpha] = getEMatrixTStiff(elemSet,nodeSet,Y,H,B,Pg,alpha0)
% Initializing Matrixes and vector


%Kg = sparse(6);
u1=Pg(1);w1=Pg(3);theta1=Pg(5);u2=Pg(2);w2=Pg(4);theta2=Pg(6);
% Other Parameters
E=Y;
A=B*H;
I=B*H^3/12;
nodeI=elemSet(1,1);
nodeJ=elemSet(1,2);
x1=nodeSet(nodeI,1);
x2=nodeSet(nodeJ,1);
y1=nodeSet(nodeI,2);
y2=nodeSet(nodeJ,2);

L0=norm([x2-x1 y2-y1]);
%             Ln=norm([x1 y1]+[u1 w1]-[x2 y2]-[u2 w2]);
Ln=norm([x1+u1-x2-u2 y1+w1-y2-w2]);
co=(x2-x1)/L0;%obj.trig(1); %cos beta0
so=(y2-y1)/L0;%obj.trig(2); %sin beta0
c=(x2+u2-x1-u1)/Ln; %cos beta
s=(y2+w2-y1-w1)/Ln; %sin beta
sin_a=co*s-so*c;
cos_a=co*c+so*s;
alpha=atan2(sin_a,cos_a);

if abs(alpha0)-pi/2>0
    switch sign(alpha0)
        case 1 % previous alpha>0
            switch sign(alpha)
                case -1
                    alpha=alpha+2*pi;
            end
        case -1
            switch sign(alpha)
                case 1
                    alpha=alpha-2*pi;
            end
    end
end
% local nodal displacement vector Pl
[u_hat theta1_hat theta2_hat] = PlCalculation(Ln,L0,theta1,theta2,alpha);
Pl=[u_hat,theta1_hat,theta2_hat]';
% local internal force vector, Fl
[N M1 M2] = FlCalculation(Pl,E,A,L0,I);
Fl=[N;M1;M2];
% global internal force vector, Fg
[Fg B] = FgCalculation(Fl,s,c,Ln);
% global tagent Stiffness Matrix, TStiffMatrix

% r = [-c -s 0 c s 0]';
% q = [s -c 0 -s c 0]';
r = [-c  c -s s 0 0]';
q = [ s -s -c c 0 0]';
Kl=KlCalculation(Pl,E,A,L0,I);
Kg=B'*Kl*B+(q*q'*N+(r*q'+q*r')*(M1+M2)/Ln)/Ln;

end


function [u_hat theta1_hat theta2_hat] = PlCalculation(Ln,L0,theta1,theta2,alpha)
u_hat=Ln-L0;
theta1_hat=theta1-alpha;
theta2_hat=theta2-alpha;
end

function [N M1 M2] = FlCalculation(Pl,E,A,L0,I)
%Pl=[u_hat,theta1_hat,theta2_hat]'
%without modified strain measure
% N=E*A/L0;
% M1=E*I/L0*E*I/L0*(4*Pl(2)+2*Pl(3));
% M2=E*I/L0*E*I/L0*(2*Pl(2)+4*Pl(3));
%with modified strain measure

% f1=Pl(1)/L0+(Pl(2)^2-Pl(2)*Pl(3)/2+Pl(3)^2)/15;
% f2=2/15*Pl(2)-Pl(3)/30;
% f3=2/15*Pl(3)-Pl(2)/30;
v215=0.133333333333333;%2/15;
f1=Pl(1)/L0+(Pl(2)*Pl(2)-Pl(2)*Pl(3)/2+Pl(3)*Pl(3))/15;
f2=v215*Pl(2)-Pl(3)/30;
f3=v215*Pl(3)-Pl(2)/30;
EA=E*A;
N=EA*f1;
EAL0f1=N*L0;
EIL0=E*I/L0;
M1=EAL0f1*f2+EIL0*(4*Pl(2)+2*Pl(3));
M2=EAL0f1*f3+EIL0*(2*Pl(2)+4*Pl(3));
end

%-------------------------------------------------------------------------
function [Fg B] = FgCalculation(Fl,s,c,Ln)
sLn=s/Ln;
cLn=c/Ln;
B=[ -c   c   -s   s   0 0;
    -sLn sLn cLn -cLn 1 0;
    -sLn sLn cLn -cLn 0 1];
Fg=B'*Fl;
end

%-------------------------------------------------------------------------
function Kl=KlCalculation(Pl,E,A,L0,I)
%without modified strain measure
% r = [-c -s 0 c s 0]';
% q = [s -c 0 s c 0]';
% b1 = r;
% b2 = [0 0 1 0 0 0]'-q/Ln;
% b3 = [0 0 0 0 0 1]'-q/Ln;
% Kl = E/L0*[A 0 0;
%     0 4*I 2*I;
%     0 2*I 4*I];
%with modified strain measure
v215=0.133333333333333;%2/15
f1=Pl(1)/L0+(Pl(2)*Pl(2)-Pl(2)*Pl(3)/2+Pl(3)*Pl(3))/15;
f2=v215*Pl(2)-Pl(3)/30;
f3=v215*Pl(3)-Pl(2)/30;
EA=E*A;
IdL0=I/L0;
AL0=A*L0;
f1215=v215*f1;%2/15*f1;
Kl=zeros(3,3);
Kl(1,1)=EA/L0;
Kl(1,2)=EA*f2;
Kl(1,3)=EA*f3;
Kl(2,1)=Kl(1,2);
% Kl(2,2)=E*A*L0*f2^2+2/15*E*A*L0*f1+4*E*I/L0;
Kl(2,2)=E*(AL0*(f2*f2+f1215)+4*IdL0);
% Kl(2,3)=E*A*L0*f2*f3   -   E*A*L0/30*f1   +    2*E*I/L0;
Kl(2,3)=E*(AL0*(f2*f3-f1/30)+2*IdL0);
Kl(3,1)=Kl(1,3);
Kl(3,2)=Kl(2,3);
% Kl(3,3)=E*A*L0*f3^2+2/15*E*A*L0*f1+4*E*I/L0;
Kl(3,3)=E*(AL0*(f3*f3+f1215)+4*IdL0);
end




