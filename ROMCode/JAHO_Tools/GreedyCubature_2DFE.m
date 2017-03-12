function [setPoints,weights,norm_resid,TIME_POINTS,setPointsALL,Jmin,bMIN] = GreedyCubature_2DFE(bMIN,Jmin,ngausG,ngaussSET,nmodesU,...
    BOUND_WEIGHTS_m_mg,DATA_CUBATURE,...
    UNSTRESSED_POINTS,weigE,setPoints)


%%%function [setPoints,weights,norm_resid,TIME_POINTS,setPointsALL,Jmin,bMIN] = GreedyCubature_2DFE(bMIN,Jmin,ngausG,nmodesU,...
%%%    BOUND_WEIGHTS_m_mg,NAME_SAVE_BasisS,NAME_SAVE_BASIS_DISP,...
%%%    nmodesSelas,TYPE_SVD_RESIDUAL,nmodesUelas,TYPE_SVD,DATA_CUBATURE,...
%%%    UNSTRESSED_POINTS,weigE,setPoints,TIME_POINTS)

format long g

%dbstop('8')
if nargin==0
    load('tmp70.mat')
    %%%ngausG = 20 ;
    % nmodesU =3 ;
    % BOUND_WEIGHTS_m_mg = 0 ;
end


%%%
%ngausG = ngausG +1 ;
%dbstop('19')


ngaus = size(Jmin, 2);

nmodes = size(Jmin,1)/nmodesU ;
if  DATA_CUBATURE.VOLUME_CONSTRAINT ==1
    nmodesF = nmodes - nmodesU ;
elseif DATA_CUBATURE.VOLUME_CONSTRAINT ==0
    nmodesF = nmodes ;
end


%%%WEIGHTED_MODES = DATA_CUBATURE.WEIGHTED_MODES ;
%%%INCLUDE_SVDISPL = DATA_CUBATURE.INCLUDE_SVDISPL ;
%%%INTENSITY_ELASTIC_FACTOR_FORCE = DATA_CUBATURE.INTENSITY_ELASTIC_FACTOR_FORCE ;
%%%INTENSITY_ELASTIC_FACTOR_DISP = DATA_CUBATURE.INTENSITY_ELASTIC_FACTOR_DISP ;
%%%VOLUME_CONSTRAINT = DATA_CUBATURE.VOLUME_CONSTRAINT  ;

%%%TYPE_SVD_FORCES = TYPE_SVD_RESIDUAL;
%%%[Jmin intfALL bMIN]= MatrixMinimization_2DFE(WEIGHTED_MODES,INCLUDE_SVDISPL,nmodes,NAME_SAVE_BASIS_DISP,...
%%%    NAME_SAVE_BasisS,Jmin,nmodesSelas,TYPE_SVD_FORCES,nmodesUelas,TYPE_SVD,nmodesU,...
%%%    INTENSITY_ELASTIC_FACTOR_FORCE,INTENSITY_ELASTIC_FACTOR_DISP,bMIN,nmodesF,VOLUME_CONSTRAINT) ;


%
%
% if ~isempty(setPoints)
%
%     Gmin = Jmin(:,setPoints) ;
%     if BOUND_WEIGHTS_m_mg == 0
%         weights = lsqnonneg(Gmin,bMIN);
%     else
%         AA = -eye(length(setPoints)) ;
%         bb = -BOUND_WEIGHTS_m_mg*ones(length(setPoints),1) ;
%         weights = lsqlin(Gmin,bMIN,AA,bb);
%     end
%
%     norm_resid = Gmin*weights-bMIN;
%     PLOT_RESID = 0 ;
%     %  dbstop('119')
%     %     if DATA_CUBATURE.REMOVE_ZERO_WEIGHTS ==1
%     %
%     %         nnnn= find(weights == 0);
%     %         setPoints(nnnn) = [] ;
%     %         weights(nnnn)=[] ;
%     %     end
%     %
% else
%
%
 
resid = bMIN ;
if size(Jmin,1) == 2
    resid =ones(size(bMIN)) ;
end
 
normb = norm(bMIN);

%%%
norm_resid = norm(resid);

resid_VECT =[] ;
PLOT_RESID = 1 ;

resid_VECTG = zeros(size(Jmin,1),ngausG) ;
%dbstop('62')surf
tol_bound =BOUND_WEIGHTS_m_mg  ;

AllPoints = 1:ngaus ;
CandPoints = ones(ngaus,1);



igaus=1 ;
nmodesUloc = nmodesU ;
% % sumJMIN
% if DATA_CUBATURE.SINGLE_ROW_GminH_FACTOR  > 0
%     nmodesUloc = 1 ;
% end
% if DATA_CUBATURE.CONSTRAINED_VOLUME  > 0
%     nmodesUloc = 0 ;
% end
TOL_ZEROS = 1e-14 ;
% if VOLUME_CONSTRAINT == 1
%     iini = size(Jmin,1)-nmodesUloc^2 ;
% else
    iini =  size(Jmin,1) ;
% end

% JminF =  Jmin(1:iini,:) ;
% sumJMIN = sum(abs(JminF),1);
% indZEROS = find(sumJMIN <TOL_ZEROS);
% CandPoints(indZEROS) =0 ;
setPoints = [] ;
% UNSTRESSED POINTS
% CandPoints(UNSTRESSED_POINTS) = 0 ;

rankBasisSg = 0 ;
TIME_POINTS = tic ;
% dbstop('127')
normresidINI = norm_resid ;

%TOL=10e-6; errorINT=1;

igausEFFECT=0;

error_INTEGRATION_nz = zeros(ngausG,1) ;   % Integration error with nonzero weights
igausEFF = zeros(ngausG,1) ;
while igaus <=ngausG && igausEFFECT<ngaussSET %% && errorINT>TOL %&& norm_resid>TOL
%while igaus <=500 %& norm_resid>TOL
    maxTARG = -1e20 ;
    %%% Looking for new candidates
    %disp(['igaus = ',num2str(igaus)]) ;
    
    % Jmin of candidate points
    indCandPoints = find(CandPoints==1);
    targFUNglo = Jmin(:,indCandPoints)'*resid ;
    [max_value index_max_value] = max(targFUNglo)  ;
    igausCAND = indCandPoints(index_max_value(1));
    %%%% End
    CandPoints(igausCAND) = 0 ;
    setPoints = [setPoints igausCAND] ;
    Gmin = Jmin(:,setPoints) ;
    
    if BOUND_WEIGHTS_m_mg == 0
        %OPTION.tolX = tolX;
        %weights = lsqnonneg(Gmin,bMIN,OPTION);
        weights = lsqnonneg(Gmin,bMIN);
    else
        AA = -eye(length(setPoints)) ;
        bb = -tol_bound*ones(length(setPoints),1) ;
        weights = lsqlin(Gmin,bMIN,AA,bb);
    end
    resid = bMIN-Gmin*weights ;
    resid_VECTG(:,igaus) = resid ;
    
    norm_resid = norm(resid);
    resid_VECT(end+1) = norm_resid;    
    errorINT = norm_resid/normresidINI*100 ;        
    
    igausEFFECT =   length(find(weights>0)) ;    
    igausEFF(igaus) = igausEFFECT ;
    error_INTEGRATION_nz(igaus) = errorINT ;
    disp(['igaus = ',num2str(igaus),'(nz=',num2str(igausEFFECT),')','; error  = ',num2str(errorINT),' %']) ;
        
    igaus = igaus+1 ;


end


TIME_POINTS = toc(TIME_POINTS) ;

%
% end


if  PLOT_RESID == 1
     PlotErrorEffectiveIntegration(igausEFF,error_INTEGRATION_nz)       
end



%dbstop('141')
setPointsALL = setPoints;
if DATA_CUBATURE.REMOVE_ZERO_WEIGHTS ==1
    
    nnnn= find(weights == 0);
    setPoints(nnnn) = [] ;
    weights(nnnn)=[] ;
end



%plotJmin = 1;
%dbstop('224')

%weights = weights.*sqrt(weigE(setPoints)) ;

%dbstop('173')



% sortPOINTS = 0 ;
% if sortPOINTS == 1
%     [weights IIXX] =    sort(weights,'descend') ;
%     setPoints = setPoints(IIXX);
%     figure(7890)
%     bar(weights)
%     xlabel('Points')
%     ylabel('Weights')
%
%
% end





%%%
