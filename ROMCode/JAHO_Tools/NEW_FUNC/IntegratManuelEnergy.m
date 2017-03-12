function [z,w]=IntegratManuelEnergy(PHI_ENERGY,ponder_factors,nGP)

%clc
%clear all
%format long g
%close all

% Reduced-order quadrature of polynomials 
% J.A. Hernández, 3-Dec-2015
%------------------------------------
% Construction of matrix J and vector b  for energy functions 
% M = total number of integration points 
% pPOL = Polynomial order
% --------------------------------------

%setPoints=[];
%weights=[];

%nameWS = 'EnergyModes.mat' ; 
factorLEQ = 1;  
 
TOL = 1e-10 ; 

%load(nameWS) ; 

Xf =  PHI_ENERGY;
 %Lambda = sqrt(Sing_Val_DIS(1:150) ) ; 
 %Lambda(1:6) = Lambda(7)/Lambda(1)*Lambda(1:6);
W = ponder_factors;
%W =1e3*W ; 

%m = 148 ; % If m is large, then the only criterion guiding the optimization is TOL 

%Xf = bsxfun(@times,Xf',Lambda); 

[J,b,W,INTexac] = Calc_J_b_Manuel(Xf,W,factorLEQ); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Location of optimal points 
% Box 6.1 
% -------
% Initial data: J, b , TOL, m 
%----
% Set of candidate points 
M = length(W);
y=1:M ; 
% Set of choosen points 
z = [] ; 
% Vector of weights 
alpha = [] ;
% Residual vector 
r = b ; 
% Number of iterations
k = 1; 
% Number of nonzero weights 
mPOS = 0 ; 
errorGLO = [] ; 

% Jnorm = sqrt(sum(J.*J,1));

while  norm(r)/norm(b) >TOL && mPOS <=nGP  
    % 1. Compute new point
    ObjFun = J(:,y)'*r ;
    % ObjFun = ObjFun./Jnorm(y)/norm(r);
    [~, iLOC] = max(ObjFun)  ; iLOC = iLOC(1);    
    i = y(iLOC) ; 
    % 2. Move i from set y to set z 
    z = [z;i] ; 
    y(iLOC) = [] ; 
    % 3. Determine alpha by solving a NNLS 
    [alpha,~,~,EXITFLAG,outputL] =lsqnonneg(J(:,z),b);
    
    %alpha=J(:,z)\b;
    
    % 4. Update the residual 
    r = b-J(:,z)*alpha ; 
    errorGLO(k) = norm(r)/norm(b)*100 ; 
    % 5. Update mPOS and k 
    mPOS = length(find(alpha>0)) ; 
    disp(['k = ',num2str(k),', mPOS=',num2str(mPOS),' ,','; error (%) = ',...
            num2str(errorGLO(k))]) ; ;   
    
    k = k + 1; 
%m = 1e6 ; % If m is large, then the only criterion guiding the optimization is TOL 
end

INDzero = find(alpha==0) ;
alpha(INDzero) = [] ; 
z(INDzero)=[] ; 
% Reduced set of weights 
w = alpha.*sqrt(W(z)) ; 

% figure(1)
% hold on 
% xlabel('x')
% ylabel('w')
% lBAR = W(1);
% bar(x(z),w,0.1)
% AXX=[axis]; 
% AXX([1,2]) = [-1,1] ; 
% 
% axis(AXX)

%disp('WEIGHTS= ')
%w
%%%Xf = Xf' ; 
%INT_APROX = (Xf(z,:)'*w) ; 

%ERROR = norm(INTexac-INT_APROX)./norm(INTexac) ; 


%%
 

 


 


 