function [lambdaP,gap]=LMI_Aut15_th1(Ai,B,C,g,kStar,h,Mstar,M1)
% This MATLAB program checks the feasibility of LMIs from Theorem 1 of the paper 
% A. Selivanov, E. Fridman, and A. Fradkov, "Passification-based adaptive control: Uncertain input and output delays," Automatica, vol. 54, pp. 107–113, 2015.

% The program uses YALMIP parser (http://users.isy.liu.se/johanl/yalmip/)
% and SeDuMi solver (http://sedumi.ie.lehigh.edu/)

% Input: 
% Ai        - the vetices of the polytope (4) in the from Ai=[A1 A2 A3 ...]; 
% B,C       - the parameters of the system (3); 
% g         - the vector from Assumption 2; 
% kStar     - the nominal controller gain from (1); 
% Mstar, M1 - the controller gain bounds from (13). 

% Output: 
% lambdaP   - the minimum eigenvalue of P; 
% gap       - the maximum element of abs(PB-C'g). 

n=size(Ai,1); 
N=size(Ai,2)/n; % the number of vertices from (4) 

%% Decision variables 
P=sdpvar(n); 
R=sdpvar(n); 
S=sdpvar(n); 
G1=sdpvar(n,n,'f'); 
G2=sdpvar(n,n,'f'); 
G3=sdpvar(n,n,'f'); 

%% The LMIs for H
LMIs=[]; 
for i=1:N     
    A=Ai(:,n*(i-1)+1:n*i); % Choose a vertex
    for a=[-Mstar Mstar]
        for b=[-Mstar Mstar]
            for c=[-M1 M1]
                H=blkvar; 
                H(1,1)=P*(A-B*kStar*g'*C)+(A-B*kStar*g'*C)'*P+S; 
                H(1,2)=c*P*B*g'*C; 
                H(1,4)=kStar*h*P*B*g'*C; 
                H(1,5)=kStar*h*P*B*g'*C-a*h*P*B*g'*C; 
                H(1,6)=h*A'*R; 
                H(2,2)=-R; 
                H(2,3)=R; 
                H(2,4)=a*h*C'*g*B'*P-h*G2; 
                H(2,5)=-h*G1; 
                H(2,6)=h*b*C'*g*B'*R-h*kStar*C'*g*B'*R; 
                H(3,3)=-(S+R); 
                H(3,4)=h*G2; 
                H(3,5)=h*G1; 
                H(4,4)=-h^2*R; 
                H(4,5)=a*h^2*P*B*g'*C-h^2*G3'; 
                H(5,5)=-h^2*R; 
                H(6,6)=-R; 
                H=sdpvar(H); 
                LMIs=[LMIs, H<=0];  %#ok<AGROW>
            end
        end
    end
end

%% Park's conditions
Park1=[R G1; G1' R]; 
Park2=[R G2; G2' R]; 
Park3=[R G3; G3' R]; 

LMIs=[LMIs, P>=0, R>=0, S>=0, Park1>=0, Park2>=0, Park3>=0, P*B==C'*g]; 

%% Solution of LMIs
options=sdpsettings('solver','sedumi','verbose',0);
sol=optimize(LMIs,[],options); 

lambdaP=[]; 
if sol.problem==0
    primal=check(LMIs); 
    if min(primal(1:end-1))>=0 && min(primal(1:end-4))>0
        lambdaP=min(eig(double(P))); 
        gap=max(abs(double(P)*B-C'*g)); 
    end
else
    yalmiperror(sol.problem); 
end
