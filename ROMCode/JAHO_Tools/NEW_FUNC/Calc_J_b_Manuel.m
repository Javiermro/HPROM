function [J,b,weigE,INTexac] = Calc_J_b_Manuel(SNAPfNW,weigE,factorLEQ)
% Construction of matrix J and vector b  for polynomial functions
% M = total number of integration points 
% pPOL = Polynomial order
% --------------------------------------
 
% --------------
H  =sqrt(weigE) ; 
% Matrix of snapshots of the integrand  (with no weights)
 %SNAPfNW
% "Exact" integrals (numerical evaluation)
INTexac = SNAPfNW'*weigE ;
% Total  
V = sum(weigE) ;
% Matrix of snapshots of snapshots with zero integral 
Xf = zeros(size(SNAPfNW)) ; 
for i=1:size(SNAPfNW,2)
    Xf(:,i) =   SNAPfNW(:,i) - INTexac(i)/V ; 
    Xf(:,i) =  bsxfun(@times, Xf(:,i), sqrt(weigE)) ;    
end
 
% SVD of Xf 
[Lambda,SValues,~] = svd(Xf,0) ;
% -------------------------------
SValues = diag(SValues) ; 
% Rank 
tol = max(size(Xf)*eps(max(SValues))); % defition of rank
RankXf = sum(SValues> tol);  % Rank Xf
Lambda = Lambda(:,1:RankXf);

        % Matrix J appearing in the minimization problem 
Jhat = Lambda' ; 
n=1 ; 
Jw = factorLEQ*sqrt(weigE'*n)/sqrt(V) ; 
J = [Jhat; Jw] ; 
% Vector b
b = zeros(size(J,1),1) ; 
b(end) = factorLEQ*sqrt(n*V);