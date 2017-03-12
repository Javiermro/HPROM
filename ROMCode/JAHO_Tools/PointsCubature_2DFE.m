function [setPoints,weights,TIME_Jmin,TIME_POINTS] =PointsCubature_2DFE(nmodesU,npoints,setPoints,...
    BOUND_WEIGHTS_m_mg,DATA_CUBATURE,ngaus,weigE,Gmin,bMIN)

dbstop('6')
if nargin==0
    load('tmp2.mat')
end

TIME_POINTS = [] ;
TIME_Jmin = tic ;
% Matrix Jmin
%ngaus = size(BasisS,1)/4 ;

% Hid operator
% ---------------------------------------
%name_WS_FE = ['DATA/',NAME_BASIS{1},'_FEoper.mat'];
%load(name_WS_FE,'OPERfe') ;
%weigE = OPERfe.weigEM ;

% PONER PESOS EN LOS PUNTOS DE GAUSS!!!!
% weigE = ....


vol= sum(weigE) ;
Hid = zeros(nmodesU*ngaus,nmodesU);
disp('Constructing Hid')
VECTORIZED =0 ;
Hid = repmat(eye(nmodesU),ngaus,1) ;

weigEMaux =   repmat(weigE',nmodesU,1) ;
SQweigEM  = reshape(sqrt(weigEMaux),nmodesU*ngaus,1);
Hid = bsxfun(@times,Hid,SQweigEM) ;
 
if  DATA_CUBATURE.VOLUME_CONSTRAINT ==1
    Hid  = (1/vol*Hid)  ;
end
disp('Done')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jmin 
%if DATAGEN.RELOAD_BASIS_FORCES == 0
    
    %BasisF =     BasisF_fromStress_2DFE(BasisS,BBfRI,DATA_CUBATURE,NAME_BASIS) ;
%    save(NAME_SAVE_BasisS,'BasisF','-append')
%else
    
%    load(NAME_SAVE_BasisS,'BasisF')
%end
        
    
disp(' Calculating Jmin ... ')

if  DATA_CUBATURE.VOLUME_CONSTRAINT ==1
    BasisFex = [BasisF Hid] ;
    nmodesFINTex = size(BasisFex,2) ;
    
elseif DATA_CUBATURE.VOLUME_CONSTRAINT ==0
    BasisFex = [BasisF ] ;
    nmodesFINTex = size(BasisFex,2) ;
    
else
    error('Option not implemented')
end

Jmin = zeros(nmodesU*nmodesFINTex,ngaus);

for imode=1:nmodesU
    imodeG = imode:nmodesU:nmodesU*ngaus ;
    JminLOC = BasisFex(imodeG,:)' ;
    
    igausG = imode:nmodesU:nmodesU*nmodesFINTex ;
    Jmin(igausG,:) = JminLOC ;
    
end

disp('Done')
%
% end
TIME_Jmin = toc(TIME_Jmin)  ;

% Calculating rank
%dbstop('72')
if  DATA_CUBATURE.VOLUME_CONSTRAINT ==1
    ngausGmin =  nmodesU*(nmodesFINTex-nmodesU)+1  ;
else
    ngausGmin =  nmodesU*nmodesFINTex ;
end

if isempty(npoints)
    ngausG =   ngausGmin ;
else
    ngausG = npoints ; % max(ngausGmin,npoints) ;
end
%TOLERANCE = 1e-10 ;

%aTOL = [] ;
% bMIN
if  DATA_CUBATURE.VOLUME_CONSTRAINT ==1
    bMIN = zeros(size(Jmin,1),1);
    iini = size(Jmin,1)-nmodesU^2+1 ;
    EE = eye(nmodesU);
    EE =reshape(EE,nmodesU^2,1);
    bMIN(iini:end)  = EE ;
elseif DATA_CUBATURE.VOLUME_CONSTRAINT ==0
     bMIN = zeros(size(Jmin,1),1);
     for imode = 1:nmodesFINTex
         iini = (imode-1)*nmodesU+1 ; 
         ifin = nmodesU*imode ;
         bMIN(iini:ifin) = Hid'*BasisFex(:,imode); 
     end 
else
    error('Option not implemented')
end

%UNSTRESSED_POINTS =  find(sum(abs(BasisS),2)<1e-16);
UNSTRESSED_POINTS =  [];

% dbstop('167')
% switch POINTS_CUBATURE
%     case 'GREEDY'

% % name_WS_POINTS = ['DATA/',NAME_BASIS{1},'_POINTSquad_',num2str(ngausG),'.mat'];
% % name_WS_POINTSmin = ['DATA/',NAME_BASIS{1},'_POINTSquad_',num2str(ngausGmin),'.mat'];
% % %dbstop('108')
% % if  isempty(setPoints)
% %     if exist(name_WS_POINTS) ==2 & DATA_CUBATURE.RELOAD_POINTS == 1
% %         load(name_WS_POINTS,'setPointsALL','TIME_POINTS')
% %         setPoints = setPointsALL;
% %     elseif exist(name_WS_POINTSmin) ==2 & DATA_CUBATURE.RELOAD_POINTS == 1 & ngausG <= ngausGmin
% %         load(name_WS_POINTSmin,'setPointsALL','TIME_POINTS')
% %         setPoints = setPointsALL(1:ngausG);
% %     
% %     end
% % end
    
[setPoints,weights,norm_resid,TIME_POINTS,setPointsALL,Jmin,bMIN] = ...
    GreedyCubature_2DFE(bMIN,Jmin,ngausG,nmodesU,BOUND_WEIGHTS_m_mg,DATA_CUBATURE,UNSTRESSED_POINTS,...
    weigE,setPoints,TIME_POINTS);

   % dbstop('132')
    save(name_WS_POINTS,'setPointsALL','TIME_POINTS')

%% Checking how good is the approximation of the total volume of the body
%dbstop('195')

disp('------------------------------------')
disp('ERROR IN APPROXIMATING TOTAL VOLUME')
disp('------------------------------------')
errorV = abs((sum(weights)-sum(weigE))/sum(weigE))*100 ;
disp(['ERROR =',num2str(errorV),' %']);
disp('---------------------------------------')
pause(1.5)

%%%%
 

%%%%

%CalculateApproxVolume_2DFE(NAME_BASIS,setPoints,weights,DATA_CUBATURE,weigE)

 
% % if DATA_CUBATURE.PLOT_WEIGHTS  ==1
% %     [weightsSORT IIXX] =    sort(weights,'descend') ;
% %     setPointsSORT = setPoints(IIXX);
% %     %  dbstop('176')
% %     setElem =   funFPLOTJMIN_2DFE(nmodesU,NAME_BASIS,weightsSORT,setPointsSORT)  ;
% %     clipboard('copy',num2str(setElem));
% %     disp(setElem)
% %     disp(['weights=',num2str(weightsSORT')])
% %     disp('Press enter to continue... (after pasting the set of number of elements in Gid window)')
% %     pause ;
% %     
% %     
% % end



save('points.mat','setPoints')
% end


