function [z,w,errorGLO,b,J]= OptimizedQuadratureMultiREST(SNAPfNW,W,TOL,m,DATA,z,J,b,DATA_CUBATURE)

%dbstop('4')
if nargin == 0
    load('tmp1.mat')
end
% --------------------------------------
COMPTIME = tic ;
if isempty(J)
    DATA.TOL = TOL ;
    [J,b,DATA] = Jmin_bmin_calc(SNAPfNW,W,DATA,DATA_CUBATURE);
end

% Initial data: J, b , TOL, m
%----
% Set of candidate points
M = size(J,2) ;
y=1:M ;
y(z) = []  ;
% Vector of weights
alpha = [] ;
% Residual vector
r = b ;
% Number of iterations
k = 1;
% Number of nonzero weights
mPOS = 0 ;
errorGLO = [] ;

 if DATA.STANDARD_APPROACH ==1
     DATA.NNEGLSQ_ONLY_WHEN_NECESSARY = 0 ;
 end

acum_iter = 0 ; 


if isempty(m)
    m=1e9 ; 
end

while  norm(r)/norm(b) >TOL && mPOS <=m
    % 1. Compute new point
    ObjFun = J(:,y)'*r ;
    %  dbstop('38')
    if DATA.DivideByNorm == 1
        ObjFun = ObjFun./DATA.Jnorm(y)'/norm(r) ;
    end
    
    [maxLOC, indSORT] = max(ObjFun)  ;
    if DATA.NNEGLSQ_ONLY_WHEN_NECESSARY ~=0
        i = y(indSORT(1)) ;
        z = [z;i] ;
        y(indSORT(1)) = [] ;
        outputL.iterations = 1 ;
        idxNEGalpha = [] ;
        alpha =  J(:,z)\b ;
        idxNEGalpha = find(alpha<=0) ;    %    nmaxNEG = 20 ;
        if ~isempty(idxNEGalpha) & DATA.NNEGLSQ_ONLY_WHEN_NECESSARY==1
         %  dbstop('52')
            [alpha,~,~,EXITFLAG,outputL] =lsqnonneg(J(:,z),b);
        elseif ~isempty(idxNEGalpha) & DATA.NNEGLSQ_ONLY_WHEN_NECESSARY==-1
            warning(['Negative weights = ',num2str(length(idxNEGalpha))])
        end
    elseif DATA.NNEGLSQ_ONLY_WHEN_NECESSARY ==0
        i = y(indSORT(1)) ;
        z = [z;i] ;
        y(indSORT(1)) = [] ;
        [alpha,~,~,EXITFLAG,outputL] =lsqnonneg(J(:,z),b);
    end
    
    % 4. Update the residual
    r = b-J(:,z)*alpha ;
    
    % 5. Update mPOS and k
    mPOS = length(find(alpha>0)) ;
  
    if DATA.NOSHOWITERS == 0
        disp(['k = ',num2str(k),', mPOS=',num2str(mPOS),' ,','; error n(res)/n(b) (%) = ',...
            num2str(norm(r)/norm(b)*100),'; error max  ',num2str(max(abs(r(1:end-1)))),'; niterLOC =',num2str(outputL.iterations)]) ; ;
    end
    acum_iter = acum_iter + outputL.iterations ; 
    
    %  if  ~isempty(idxNEGalpha)
    INDzero = find(alpha==0) ;
    alpha(INDzero) = [] ;
    % dbstop('76')
    y = [y z(INDzero)'];
    z(INDzero)=[] ;
    
     
    comp = length(z) ;
  %  dbstop('93')
      errorGLO(comp) = norm(r)/norm(b)*100 ;
    
    % end
    k = k + 1;
end

if DATA.NOSHOWITERS == 1
    disp(['k = ',num2str(k-1),', mPOS=',num2str(mPOS),' ,','; error n(res)/n(b) (%) = ',...
        num2str(errorGLO(k-1)),'; error max  ',num2str(max(abs(r(1:end-1)))),'; niterLOC =',num2str(outputL.iterations)]) ; ;
end
INDzero = find(alpha==0) ;
alpha(INDzero) = [] ;
z(INDzero)=[] ;
% Reduced set of weights
switch DATA.vectWEIGHTS_inJ
    case 'SQUARE'
        w = alpha.*sqrt(W(z)) ;
    case {'NONE','NO'}
        w = alpha ; %.*sqrt(W(z)) ;
    case 'ALL'
        w = alpha.*(W(z)) ;
    otherwise
        error('Option not implemented')
end
COMPTIME = toc(COMPTIME) ; %  = tic ;
disp(['COMP. TIME =',num2str(COMPTIME), 's'])
disp(['Total number of iterations =',num2str(acum_iter)])
