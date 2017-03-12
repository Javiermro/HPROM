function [J,b,DATA] = Jmin_bmin_calc(SNAPfNW,weigE,DATA,DATA_CUBATURE)
% Construction of matrix J and vector b
% --------------------------------------
%dbstop('5')
if nargin == 0
    load('tmp1.mat')
end
nsnap = size(SNAPfNW,2);
INTexac = zeros(1,nsnap) ;
M = size(SNAPfNW,1) ;

INTexac = SNAPfNW'*weigE ;

% Total
V = sum(weigE) ;

% Matrix of snapshots of snapshots with zero integral
%Xf = zeros(size(SNAPfNW)) ;

if DATA.STANDARD_APPROACH == 1
    DATA.vectWEIGHTS_inJ = 'NO' ;
end

%dbstop('25')
switch   DATA.vectWEIGHTS_inJ
    case 'NO'
        Xf = SNAPfNW ;
        Jw = []; ; %
    case 'SQUARE'
        Xf = bsxfun(@minus,SNAPfNW',INTexac/V)' ;
        Xf = bsxfun(@times,Xf,sqrt(weigE)) ;
       
        % Vector b
        wMUlT = sqrt(weigE) ;
        if DATA.unitaryBvector==1
            b_end = 1 ; 
           Jw = sqrt(weigE')/V;
        else
             Jw = sqrt(weigE')/sqrt(V) ;
             b_end = sqrt(V);

        end
    case 'ALL' % ; 'SQUARE' ;  % ALL, NONE
        %         for i=1:size(SNAPfNW,2)
        %             average = INTexac(:,i)/V  ;
        %           %  average = repmat(average,M,1) ;
        %             Xf(:,i) =   SNAPfNW(:,i) -average(i) ;
        %             Xf(:,i) =   Xf(:,i).*(weigE) ;
        %         end
        Xf = bsxfun(@minus,SNAPfNW',INTexac/V)' ;
        Xf = bsxfun(@times,Xf,weigE) ;
        
        Jw =weigE' ;
        b_end = V;
        wMUlT = ones(size(Jw')) ;
    case 'NONE' % ;
        % dbstop('48')
        %           average = INTexac/V  ;
        %         for i=1:size(SNAPfNW,2)
        %
        %            % average = repmat(average,M,1) ;
        %             Xf(:,i) =   SNAPfNW(:,i) -average(i) ;
        %             %  Xf(:,i) =   Xf(:,i); %.*sqrt(weigE) ;
        %         end
        Xf = bsxfun(@minus,SNAPfNW',INTexac/V)' ;
        %  Xf = bsxfun(@times,Xf,weigE)' ;
        
        Jw = ones(size(weigE')); ; %  sqrt(weigE'*n)/sqrt(V) ;
        wMUlT = weigE ;
        % Vector b
        %b = zeros(size(J,1),1) ;
        b_end= V;
        
end

DATA = DefaultField(DATA,'TOL_REL_J',1e-12) ;

DATA = DefaultField(DATA,'MAKESVD',1);

 DATA = DefaultField(DATA,'nmodesFint_from_snapMfNL',[]) ;
 
if DATA.STANDARD_APPROACH ==1 | DATA.MAKESVD==0
    Lambda = Xf ;
    nmodesFINTex = size(Lambda,2) ;
else
    disp('SVD snapshot matrix')
  %  dbstop('75')
    [Lambda,S,~] = svd(Xf,0) ;
    % -------------------------------
    S = diag(S) ;
    %dbstop('78')
    ENERGYCAL = 0 ;
   
   
        if  ENERGYCAL == 1
            %  Energy associated to each norm
            Scuad_rev = S(length( S ):-1:1).^2 ;
            Scuad_norm = sum(Scuad_rev) ;
            CSS =  cumsum(Scuad_rev)/Scuad_norm ;
            enngy = sqrt(CSS(end:-1:1))*100;
            % Trunctation level %
            RankXf = length(find(enngy>DATA.TOL_REL_J)) ;
        else
            RankXf = length(find((S/S(1)*100)>DATA.TOL_REL_J));
            
        end
        if  ~isempty(DATA.nmodesFint_from_snapMfNL)
            
            RankXf = min(DATA.nmodesFint_from_snapMfNL,RankXf) ; 
        end
%     
%     if ~isempty(DATA.npoints)
%         RankXf = min(DATA.npoints,RankXf) ; 
%     end
    
    
    Lambda = Lambda(:,1:RankXf);
    nmodesFINTex = RankXf ;
    
    %dbstop('114')
    if DATA_CUBATURE.MULTIPLY_BASISF_by_SINGVAL ==1 
        S =S(1:RankXf) ; 
       Lambda =  bsxfun(@times,Lambda',S)' ; 
    end
    
end

%dbstop('107')
INTexac = zeros(1,nmodesFINTex) ;
INTexac  =  Lambda'*weigE ;



% Matrix J appearing in the minimization problem


Jmin = Lambda' ;

%dbstop('103')
if ~isempty(Jw)
    J = [Jmin; Jw] ;
    b = zeros(size(J,1),1) ;
    %b(1:end-1) = Jmin*wMUlT ;
    b(end) = b_end ;
else
    J = Jmin ;
    b = INTexac(:) ;
end

% 
if  DATA.DivideByNorm ==1 
    DATA.Jnorm = sqrt(sum(J.*J,1)) ; 
else
    DATA.Jnorm = [] ; 
end