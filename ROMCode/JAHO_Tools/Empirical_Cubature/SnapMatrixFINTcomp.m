function    [ SNAPfNW weigM BasisS zSTRESSES]=SnapMatrixFINTcomp(BasisS,BBfRI,...
    weigM,nstrain,SingValStr)

%dbstop('5')
if nargin == 0
    load('tmp3.mat')   %
end
zSTRESSES = [] ; 
%
nmodesU = size(BBfRI,2) ;
nmodesS = size(BasisS,2) ;

% Singular values for elastic/inelastic
%[SSd,SS] =RetrieveSingV(NAME_SAVE_BASIS_DISP,NAME_SAVE_BasisS,INCLUDE_SVDISPL) ;

SSd=SingValStr.SV_EPS;
%SingVal = SS(1:nmodesS) ;
SingVal=SingValStr.SV_STRESS(1:nmodesS);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

%name_WS_FE = ['DATA/',NAME_BASIS{1},'_FEoper.mat'];
%load(name_WS_FE,'OPERfe') ;
%weigM = OPERfe.weigEM ;
%nstrain = OPERfe.nstrain ;
ngaus = size(BasisS,1)/nstrain;
volF = sum(weigM) ;

%%% BasisS multiplied by singular values
% for imodes = 1:nmodesS
%     BasisS(:,imodes) = SingVal(imodes)* BasisS(:,imodes)  ;
% end

BasisF = zeros(ngaus*nmodesU,nmodesS) ;
BBfRIt = BBfRI' ;
%BstT_BasisS  = BBfRIt*BasisS ;
disp('--------------------------------------------')
disp('COMPUTING BASIS FOR INTERNAL FORCES....(arranged as vectors)')
disp('--------------------------------------------')

%DATAGEN = DefaultField(DATAGEN,'INCLUDE_Bs_in_FINT',0) ;
%matRef=zeros(15,10);
%PgMAT = zeros(size(BBfRIt,1),size(BasisS,2));
%SNAPfNW = zeros(nmodesU*nmodesS,ngaus);
for igaus=1:ngaus
    iini = (igaus-1)*nmodesU+1;
    ifin =  igaus*nmodesU ;
    ig = (igaus-1)*nstrain+1 ;
    fg = (igaus*nstrain);
    
    %BasisF(iini:ifin,:) =(1/weigM(igaus))*BBfRIt(:,ig:fg)*BasisS(ig:fg,:) ;
    BasisF(iini:ifin,:) =BBfRIt(:,ig:fg)*BasisS(ig:fg,:) ;
    
    %PgMAT=(1/weigM(igaus))*BBfRIt(:,ig:fg)*BasisS(ig:fg,:) ;
    %PgMAT=BBfRIt(:,ig:fg)*BasisS(ig:fg,:);
    
    %SNAPfNW(:,igaus) = PgMAT(:);
    %matRef = matRef + weigM(igaus)^2*BBfRIt(:,ig:fg)*BasisS(ig:fg,:);
    aaaa =mod(igaus,100) ;
    if aaaa==0
        disp(['igaus=',num2str(igaus),' of ',num2str(ngaus)])
    end
end

%% Now we re-arrange the above matrix in a componentwise fashion
% --------------------------------------------------------------
%dbstop('77')
SNAPfNW = ArrangeByComponentSVD(BasisF,nmodesU,SSd)   ;
%SNAPfNW = SNAPfNW';

disp('DONE')

