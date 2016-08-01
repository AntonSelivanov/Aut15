% This MATLAB program checks the feasibility of LMIs from Theorems 1 of the paper 
% A. Selivanov, E. Fridman, and A. Fradkov, "Passification-based adaptive control: Uncertain input and output delays," Automatica, vol. 54, pp. 107–113, 2015.

%% Parameters from Example 1
a1=[.1 1.5]; a2=[27 52]; 
Ai=[a1(1) 1 0   a1(1) 1 0   a1(2) 1 0   a1(2) 1 0; 
    a2(1) 1.3 0 a2(2) 1.3 0 a2(1) 1.3 0 a2(2) 1.3 0; 
    0 1 0       0 1 0       0 1 0       0 1 0]; 
B=[19/15; 19; 0]; 
C=[0 1 0; 0 0 1]; 
g=[1;1]; 
kStar=4.61; 

% Theorem 1
h=1e-3; Mstar=5; M1=.4; 
[lambdaP,gap]=LMI_Aut15_th1(Ai,B,C,g,kStar,h,Mstar,M1); 

if ~isempty(lambdaP)
    disp(['The LMIs are feasible with max|PB-C''g|=' num2str(gap)]); 
    disp(['h1 bound=' num2str(M1*lambdaP/(Mstar^2*norm(g'*C)^2))]); 
else
    disp('The LMIs are not feasible'); 
end