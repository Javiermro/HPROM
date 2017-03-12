function [setPoints,weights,TIME_Jmin,TIME_POINTS,wSTRESSES] =...
    EmpiricalCubatureMAnuel(npointsINPUT,BBfRI,BasisS,DATA_CUBATURE,...
    NAME_BASIS,weigM,ntens,SingVal) 

%load('DATAMANUEL/EmpiricalCubatureDATA.mat')

npoints = npointsINPUT ; 
DATA_CUBATURE.INTEGRATE_STRESSES   =1 ; 

% -----------------------------------
% SNAPSHOT MATRIX (INTERNAL FORCES COMPONENT BY COMPONENT)
%- -----------------------------------
%dbstop('18')
%DATAGEN.RELOAD_BASIS_FORCES = 0 ;
setPointsini = [] ;
%BasisS = [] ;
% dbstop('29')
if  DATA_CUBATURE.STRAIN_ENERGY_METHOD == 0
    [ SNAPfNW W BasisS setPointsini] =     SnapMatrixFINTcomp(BasisS,BBfRI,...
        weigM,ntens,SingVal) ;
elseif  DATA_CUBATURE.STRAIN_ENERGY_METHOD == 1
    % Xavier's method (strain energy)
    W = weigM;
    SNAPfNW = BasisS;
    %[ SNAPfNW W ] =     SnapMatrixENER(DATA_CUBATURE,NAME_BASIS,...
    %    NAME_SAVE_BasisS,NAME_SAVE_BASIS_DISP,DATA_CUBATURE.INCLUDE_SVDISPL,DATAGEN,...
    %    DATA1st) ;
end

%DATAGEN.NAME_SAVE_BasisS = NAME_SAVE_BasisS ;
% end
%dbstop('45')
%name_WS_POINTS = ['DATA/',NAME_BASIS{1},'_POINTSquad_',num2str(size(SNAPfNW,2)),'.mat'];


%if exist(name_WS_POINTS) ==2 && DATA_CUBATURE.RELOAD_POINTS == 1
%    load(name_WS_POINTS,'setPoints','weights','J','b') ;
%else
DATA.zINI = setPointsini ;

DATA.vectWEIGHTS_inJ = 'SQUARE' ;
DATA_CUBATURE = DefaultField(DATA_CUBATURE,'NNEGLSQ_ONLY_WHEN_NECESSARY',1);
DATA.NNEGLSQ_ONLY_WHEN_NECESSARY = DATA_CUBATURE.NNEGLSQ_ONLY_WHEN_NECESSARY;
ndom = 1; nperm = [] ;
TOL = 1e-12 ;
DATA_CUBATURE = DefaultField(DATA_CUBATURE,'TOL_REL',1e-10);
DATA.TOL_REL_J = DATA_CUBATURE.TOL_REL ;
DATA.npoints  = npoints ;
DATA_CUBATURE = DefaultField(DATA_CUBATURE,'DivideByNorm',1) ;
DATA.DivideByNorm = DATA_CUBATURE.DivideByNorm  ;
DATA_CUBATURE = DefaultField(DATA_CUBATURE,'AN2009_APPROACH',0) ;
DATA.STANDARD_APPROACH =  DATA_CUBATURE.AN2009_APPROACH  ;
DATA_CUBATURE = DefaultField(DATA_CUBATURE,'nmodesFint_from_snapMfNL',[]) ;
DATA.nmodesFint_from_snapMfNL = DATA_CUBATURE.nmodesFint_from_snapMfNL ;
DATA_CUBATURE = DefaultField(DATA_CUBATURE,'unitaryBvector',0) ;
% DATA.MULTIPLY_BASISF_by_SINGVAL =     DATA_CUBATURE.MULTIPLY_BASISF_by_SINGVAL ;

DATA_CUBATURE = DefaultField(DATA_CUBATURE,'MULTIPLY_BASISF_by_SINGVAL',0) ;
[setPoints,weights,errorGLO,J,b] = OptimizedQuadraturePARTITIONED(SNAPfNW,ndom,nperm,W,DATA,TOL,DATA_CUBATURE) ;

%save(name_WS_POINTS,'setPoints','weights') ;
%save(NAME_SAVE_BasisS,'J','-append') ;
%end
TIME_Jmin=[] ; TIME_POINTS =[] ;

%DATAGEN.name_WS_POINTS = name_WS_POINTS ;

%dbstop('85')
if DATA_CUBATURE.PLOT_WEIGHTS  ==1
    [weightsSORT IIXX] =    sort(weights,'descend') ;
    setPointsSORT = setPoints(IIXX);
    %  weightsSORT = weights ;
    % setPointsSORT = setPoints ;
    figure(1078)
    hold on
    xlabel('Point')
    ylabel('Weight')
    % dbstop('94')
    bar(weights)
    % dbstop('96')
    setElem =   funFPLOTJMIN_2DFE(nmodesU,NAME_BASIS,weightsSORT,setPointsSORT)  ;
    clipboard('copy',num2str(setElem));
    disp(setElem)
    disp(['weights=',num2str(weightsSORT')])
    disp('Press enter to continue... (after pasting the set of number of elements in Gid window)')
    pause ;
    
    
end

%INTEGRATION OF STRESSES
%dbstop('100')
DATA_CUBATURE = DefaultField(DATA_CUBATURE,'INTEGRATE_STRESSES',0) ;
if DATA_CUBATURE.INTEGRATE_STRESSES  ==1  %& DATA_CUBATURE.STRAIN_ENERGY_METHOD == 0
    % Exact integral of BasisS
    nstrain = size(BasisS,1)/length(W) ;
    BasisS = ArrangeByComponentSVD(BasisS,nstrain,[])   ;
    BasisSw = bsxfun(@times,BasisS,sqrt(W)) ;
    b = BasisSw'*sqrt(W) ;
    [alpha,~,~,EXITFLAG,outputL] =lsqnonneg(BasisSw(setPoints,:)',b);
    resid = norm(b-BasisSw(setPoints,:)'*alpha);
    disp(['ERROR INTEGRATION STRESS %',num2str(resid/norm(b)*100)])
    wSTRESSES = alpha.*sqrt(W(setPoints)) ;
    %save(name_WS_POINTS,'wSTRESSES','-append')
    %DATAGEN.name_WS_POINTS = name_WS_POINTS ;
    %else
    %   error('Option not implemented')
end


