function ROM_II = HROM_OFFLINE_PROC(e_DatSet,ModoPHI_REG,ModoPHI_DIS,BASIS_DATA,...
    BASIS_EXPANSION,NBASIS_ENER_PREV_REG,NBASIS_ENER_PREV_DIS,ALTERNATE_BASIS_STRESS_B,NGAUSS_GAPPY_loc,...
    file,GREEDY_ALGOR,SNAPSHOTS_ORIGIN,CHATUR_NONREPEATED,INDICES_GIVEN_LOC,SET_GPS,PARAMETER_GREEDY,...
    LOAD_RESULT_GREEDY,COMPARE_WITH_SNAPSHOTS,nModesEPS_bar,nModesEPS_bar2,PARTIAL_RECONSTR,WHERE_WEIGHTS,...
    LEVEL_TRUN_INTERN,EnergySETS,Set_BANDAS,e_VG)

% Reconocimiento del set multiescala
% **********************************
%conshyp_GEN = zeros(e_VG.nSet,1) ;
%for iConshyp = 1:e_VG.nSet
%    conshyp_GEN(iConshyp) = e_DatSet(iConshyp).e_DatMat.conshyp ;
%end

% **********************************
% Definicion del set multiescala 
% Nota: en este caso inicialmente solo se considerara un posible tipo de microcelda..
% **********************************
%ME_Set=find(conshyp_GEN==53);
e_VG_ME = e_DatSet.e_DatMat.e_VG;
nSet = e_DatSet.e_DatMat.e_VG.nSet;
ntens = e_VG_ME.ntens ;
nGTOT=0;
for iSet = 1:nSet
    nGTOT = nGTOT + e_DatSet.e_DatMat.e_DatSet(iSet).e_DatElem.npg*...
        e_DatSet.e_DatMat.e_DatSet(iSet).nElem;
end

BmatRIstT  = zeros(e_VG_ME.nModesEPS_TOT,ntens*nGTOT);
BmatRIst_W = cell(nGTOT,1); % (% BmatRIst, not affected by WEIGHTS)

IntPhiTGen_REG = zeros(e_VG_ME.nModesEPS_bar,ntens-1);
IntPhiTGen_DIS = zeros(e_VG_ME.nModesEPS_bar2,ntens-1);
%IntPhiTGen  = zeros(e_VG_ME.nModesEPS_TOT,ntens-1);

% General B Matrix (for stabilization purposes)
BmatGEN = zeros(ntens*nGTOT,e_VG_ME.ndoft);
BmatGEN_W = zeros(ntens*nGTOT,e_VG_ME.ndoft);

KG_WEAK = zeros(e_VG_ME.ndoft);

dof_FLAG = false(e_VG_ME.ndoft,nSet);
SET_ENER_COMP = false(nGTOT,nSet);

ponder_factors = zeros(nGTOT,1) ;
ponder_factors_SET = zeros(nGTOT,nSet);

% DEFINICION DE SET DE LAS BANDAS
Set_BANDAS = 1;

COMP_STR_REG = []; ElemREG = [];
COMP_STR_DIS = []; ElemDIS = [];

K_CONVEX_REG = zeros(e_VG_ME.nModesEPS_bar) ;
K_CONVEX_DIS = zeros(e_VG_ME.nModesEPS_bar2) ;

for iSet = 1:nSet
    e_Dat_iSet = e_DatSet.e_DatMat.e_DatSet(iSet);
    nElem = e_Dat_iSet.nElem;
    eltype = e_Dat_iSet.e_DatElem.eltype;
    npg = e_Dat_iSet.e_DatElem.npg;

    % B Matrix by set
    m_BT_SET = e_DatSet.e_DatMat.e_DatSet(iSet).m_BT;
    
    % Element id
    m_NumElem =  e_DatSet.e_DatMat.e_DatSet(iSet).m_NumElem ;
    
    %for integration of the weak problem (recover displacements from strains)
    m_DetJT =   e_DatSet.e_DatMat.e_DatSet(iSet).m_DetJT;
    wg = e_DatSet.e_DatMat.e_DatSet(iSet).e_DatElem.wg;
    m_pesoPG = bsxfun(@times,m_DetJT,wg);
    
    % Loop over the elements in each Set
    for iElem=1:nElem
        dof = f_DofElem(e_Dat_iSet.conec(iElem,:),e_VG_ME.ndn);
        ponder_factors((m_NumElem(iElem)-1)*npg+1:m_NumElem(iElem)*npg) = e_Dat_iSet.m_VolElemGP(:,iElem);
        
        % Por set para determinacion de integradores reducidos de tensiones
        ponder_factors_SET((m_NumElem(iElem)-1)*npg+1:m_NumElem(iElem)*npg,iSet) = e_Dat_iSet.m_VolElemGP(:,iElem);

        % Set sobre el tipo de elementos en la microescala
        switch eltype
            case 2
            %    for ipg = 1:npg
            %        ilocp = (m_NumElem(iElem)-1)*npg+ipg ;
            %        %BmatRIT =   ModoPHI(dof,:)'*e_Dat_iSet.m_BT(:,:,ipg,iElem)' ;
            %        
            %        % Definint the ilocp entry of BmatRIstT
            %        iini = (ilocp-1)*ntens +1;
            %        ifin = ilocp*ntens ;
            %        
            %        BmatRIT = PHI_GEN(iini:ifin,:)' ;                    
            %        %  WHERE_WEIGHTS = 'WinB' ;
            %        % WinB: Weighting factors are included in B-matrix (default)
            %        % WinS: Weighting factors are included in SNAPSHOTS
            %        BmatRIstT(:,iini:ifin) = BmatRIT*e_Dat_iSet.m_VolElemGP(ipg,iElem);
            %        BmatRIst_W{(m_NumElem(iElem)-1)*npg+ipg} = BmatRIT';
            %        % CAVEAT: There is duplication of information. Take this
            %        % into consideration in future, optimized versions
            %    end                
            case {4,31,32}
                
                e_PGsID = (m_NumElem(iElem)-1)*ntens*npg+1:m_NumElem(iElem)*ntens*npg;
                e_PGsID2 = (m_NumElem(iElem)-1)*npg+1:m_NumElem(iElem)*npg ;
                
                %PHI_GEN = [ModoPHI_REG(e_PGsID,:) ModoPHI_DIS(e_PGsID,:) ModoPHI_ELAS2(e_PGsID,:)] ;
                PHI_GEN = [ModoPHI_REG(e_PGsID,:) ModoPHI_DIS(e_PGsID,:)] ;
                
                PHI_REG = ModoPHI_REG(e_PGsID,:);
                PHI_DIS = ModoPHI_DIS(e_PGsID,:);

                % Identificador de componentes de grados de libertad para cada uno de los sets
                dof_FLAG(dof,iSet)= true;
                
                % Energy components in each Set at the microscale (for Energy Decomposition)
                SET_ENER_COMP(e_PGsID2,iSet) = true;
                
                e_BT = m_BT_SET(:,:,:,iElem);
                
                if iSet==Set_BANDAS
                    ElemDIS = [ElemDIS;m_NumElem(iElem)] ;
                else
                    ElemREG = [ElemREG;m_NumElem(iElem)] ;
                end
                                
                for ipg = 1:npg
                    ilocp = (m_NumElem(iElem)-1)*npg+ipg ;
                    %BmatRIT =   ModoPHI(dof,:)'*e_Dat_iSet.m_BT(:,:,ipg,iElem)';

                    % Definint the ilocp entry of BmatRIstT
                    iini = (ilocp-1)*ntens +1;
                    ifin = ilocp*ntens ;
                    
                    %BmatRIT = PHI_GEN(iini:ifin,:)' ;
                    iPG_ID = (ipg-1)*ntens+1:ipg*ntens ;
                    
                    % For stabilization purposes
                    BmatGEN(e_PGsID(iPG_ID),dof) = e_BT(:,:,ipg) ;
                    
                    BmatGEN_W(e_PGsID(iPG_ID),dof) = e_BT(:,:,ipg)*m_pesoPG(ipg,iElem);
                    
                    BmatRIT = PHI_GEN(iPG_ID,:)' ;
                    
                    %  WHERE_WEIGHTS = 'WinB' ;
                    % WinB: Weighting factors are included in B-matrix (default)
                    % WinS: Weighting factors are included in SNAPSHOTS
                    BmatRIstTW = BmatRIT*e_Dat_iSet.m_VolElemGP(ipg,iElem) ;
                    BmatRIstT(:,iini:ifin) = BmatRIstTW;
                    BmatRIst_W{(m_NumElem(iElem)-1)*npg+ipg} = BmatRIT';
                    
                    %IntPhiTGen = IntPhiTGen + BmatRIstT(:,iini:ifin) ;
                    if iSet==Set_BANDAS
                        BmatRIT = PHI_DIS(iPG_ID,:)'*e_Dat_iSet.m_VolElemGP(ipg,iElem) ;
                        IntPhiTGen_DIS = IntPhiTGen_DIS + BmatRIT(:,[1 2 4]) ;
                        COMP_STR_DIS = [COMP_STR_DIS iini:ifin];
                        
                        % Convexification matrix - Discontinuous domain
                        K_CONVEX_DIS = K_CONVEX_DIS + PHI_DIS(iPG_ID,:)'*e_Dat_iSet.m_VolElemGP(ipg,iElem)*PHI_DIS(iPG_ID,:);
                        
                    else
                        BmatRIT = PHI_REG(iPG_ID,:)'*e_Dat_iSet.m_VolElemGP(ipg,iElem) ;                        
                        IntPhiTGen_REG = IntPhiTGen_REG + BmatRIT(:,[1 2 4]) ;
                        COMP_STR_REG = [COMP_STR_REG iini:ifin];
                        
                        % Convexification matrix - Regular domain
                        K_CONVEX_REG = K_CONVEX_REG + PHI_REG(iPG_ID,:)'*e_Dat_iSet.m_VolElemGP(ipg,iElem)*PHI_REG(iPG_ID,:);
                        
                    end
                    % CAVEAT: There is duplication of information. Take this
                    % into consideration in future, optimized versions
                    
                    %MATRICES FOR SOLVING WEAK PROBLEM
                    KG_WEAK(dof,dof)=KG_WEAK(dof,dof)+ e_BT(:,:,ipg)'*e_BT(:,:,ipg)*m_pesoPG(ipg,iElem);
                    
                end
            otherwise
            error('Element type not available')                
        end
    end
end

% Convexification matrices (for singular and regular domains)
%K_CONVEX_DIS = e_VG_ME.K_PARAM*K_CONVEX_DIS;
%K_CONVEX_REG = e_VG_ME.K_PARAM*K_CONVEX_REG;
K_CONVEX_DIS = e_VG_ME.K_PARAM*eye(size(ModoPHI_REG,2));
K_CONVEX_REG = e_VG_ME.K_PARAM*eye(size(ModoPHI_DIS,2));

%IntPhiTGen = e_VG_ME.E_MATRIX*[IntPhiTGen_REG;IntPhiTGen_DIS];
IntPhiTGen_REG_PREV = IntPhiTGen_REG;
IntPhiTGen_DIS_PREV = IntPhiTGen_DIS;

IntPhiTGen_REG = e_VG_ME.E_MATRIX*IntPhiTGen_REG;
IntPhiTGen_DIS = e_VG_ME.E_MATRIX*IntPhiTGen_DIS;

%MATRICES FOR SOLVING WEAK PROBLEM (FROM STRAINS TO DISPLACEMENTS)
ValRM = true(e_VG_ME.ndoft,1);
ValRM(e_VG_ME.doff_RM,1) = false;

Fint_WEAK = BmatGEN_W'*[ModoPHI_REG ModoPHI_DIS];

MAT_WEAK = KG_WEAK(ValRM,ValRM)\Fint_WEAK(ValRM,:);

% **************************
% EXPANSION MATRIZ SNAPSHOTS
% **************************
% JAHO_CRITICISM
% --------------
% A)
% Should not the SVD be carried over the B-matrix affected by the
% corresponding weights????
% Answer: Yes, it should !!!!
% B)
% Should not be considered the ENTIRE MATRIX (according to the expression below,
% only indicies indices_sigma_SC are considered )
% Answer to B): The third component in the STRAIN VECTOR is zero, and thus, so is
% the third row of B. Accordingly, it is not necessary to consider such a
% row (for planes stress).

% SVD of Bmat (for construction of the expanded basis) as [PHI U_B]
% U_B basis
% e_VG.BASIS_EXPANSION = DIRECT_B  NO LONGER AVAILABLE (NO ACCEPTABLE RESULTS)
%dbstop('261')

% Definition of sets (Energy decomposition)
% EnergySETS = 'DECOMP_01'; %JLM
switch EnergySETS
    case 'NO_DECOMP'
        SetCOMP{1} = any(SET_ENER_COMP,2);   % BANDAS Y MATRIZ
        
    case 'DECOMP_01'
        %Set_BANDAS = 1;
        SetCOMP{1} = any(SET_ENER_COMP(:,Set_BANDAS+1:nSet),2); % MATRIZ
        SetCOMP{2} = SET_ENER_COMP(:,Set_BANDAS);               % BANDAS
        
        % For ponder_factors (Q factor assessment)
        PGS = false(nGTOT,1);
        inBands = find(ponder_factors_SET(:,Set_BANDAS)~=0);
        PGS(inBands,1) = true;
        
        %Rearrangement of gauss weights at the small scale
        ponder_factors = {ponder_factors(~PGS) ponder_factors(PGS)};
        
        % Rotation matrix for stresses (local tension reconstruction)
        %T_BANDAS = T_GEN(SetCOMP{2},SetCOMP{2}); clear T_GEN
        
    otherwise
        error('Stress sets not defined!');
end

% ***************************************
% LECTURA DE TENSIONES - SVD DE TENSIONES
% ***************************************
disp('SVD of stress snapshot matrix')
if ~isempty(BASIS_DATA)
    switch EnergySETS
        case 'NO_DECOMP'
            %load(BASIS_DATA,'PHI_STRESS','SingleV_STRESS')
            %load(BASIS_DATA,'PHI_STRESS_REG','SingleV_STRESS')
            load(BASIS_DATA,'PHI_ENER_REG','SingleV_ENER')
            %U = PHI_STRESS ;
            %rankS = size(PHI_STRESS,2) ; clear PHI_STRESS ;
        case 'DECOMP_01'
            %load(BASIS_DATA,'PHI_STRESS_REG','PHI_STRESS_DIS','SingleV_STRESS')
            
            %load(BASIS_DATA,'PHI_ENER_REG','PHI_ENER_DIS','SingleV_ENER')
            
            load(BASIS_DATA,'PHI_ENER_REG','PHI_ENER_DIS','Sing_Val_ELAS_REG','Sing_Val_INELAS_REG','Sing_Val_ELAS_DIS','Sing_Val_INELAS_DIS')
    end
else
    %aaa = exist(SVD_STRESS_ws);
    %if load_SVD_STRESS == 0 || aaa == 0
    %    % dbstop('215')
    %    [U,S,V] = svd(SIGMA_SNAP_SC,0);
    %    rankS = rank(SIGMA_SNAP_SC);
    %    save(SVD_STRESS_ws,'U','S','V','rankS')
    %else
    %    load(SVD_STRESS_ws,'U','rankS');
    %end
    error('Decomposition not defined yet!!')
end

% Stabilization matrix, Consistent expansion for HROM Problem
%EXPAN_MATRIX = 'B_GEN_DECOMP_DOMAIN';
EXPAN_MATRIX = 'DIRECTLY_MODES';
%EXPAN_MATRIX = 'MODIFIED_TEST_FUNC_FOR_STRESSES';
switch EXPAN_MATRIX
    case 'PHI_STRAIN'
        EXP_MATRIX{1} = BmatRIstT';
        nmodes_disp = size(BmatRIstT',2);
    case 'B_GEN_STD'
        EXP_MATRIX{1} = BmatGEN ;
    case 'B_GEN_DECOMP_DOMAIN'
        EXP_MATRIX{1} = BmatGEN(SetCOMP{1},:); % COMPONENTS OUTSIDE THE BANDS
        EXP_MATRIX{2} = BmatGEN(SetCOMP{2},:); % COMPONENTS IN THE BANDS
        nmodes_disp{1} = size(BmatRIstT',2);
        nmodes_disp{2} = size(BmatRIstT',2);
    case 'DIRECTLY_MODES'
        %EXP_MATRIX{1} = ModoPHI_REG(SetCOMP{1},:); % COMPONENTS OUTSIDE THE BANDS
        %EXP_MATRIX{2} = ModoPHI_DIS(SetCOMP{2},:); % COMPONENTS IN THE BANDS 
        BmatRIstT_MOD=BmatRIstT';
        EXP_MATRIX{1} = BmatRIstT_MOD(SetCOMP{1},1:e_VG_ME.nModesEPS_bar); % COMPONENTS OUTSIDE THE BANDS
        EXP_MATRIX{2} = BmatRIstT_MOD(SetCOMP{2},e_VG_ME.nModesEPS_bar+1:end); % COMPONENTS IN THE BANDS   
        nmodes_disp{1} = size(ModoPHI_REG,2);
        nmodes_disp{2} = size(ModoPHI_DIS,2);
    otherwise
        error('Expansion not implemented, check this argument!');
end

% Treatment of stabilization matrix prior the expansion
U_B = cell(length(SetCOMP),1);

for iSET_STRESS = 1: length(SetCOMP)
    % Tolerance for divisors in normalization of snapshots in B
    tolDIV = 10e-10;
    switch BASIS_EXPANSION
        case 'DIRECT_U'
            aaa = exist(SVD_B_ws);
            if load_SVD_B == 0 || aaa == 0 ;
                %dbstop('201')
                disp('SVD of B-matrix')
                [U_B{iSET_STRESS},S,V] = svd(EXP_MATRIX{iSET_STRESS},0);
                save(SVD_B_ws,'U_B','S','V')
            else
                load(SVD_B_ws,'U_B','S','V');
            end
        case 'DIRECT_B'
            %if ALTERNATE_BASIS_STRESS_B>1 && NBASIS_STRESS_PREV==size(BmatRIstT',2) && NBASIS_STRESS_PREV<=ALTERNATE_BASIS_STRESS_B
            if ALTERNATE_BASIS_STRESS_B>1 && NBASIS_ENER_PREV==size(EXP_MATRIX{iSET_STRESS},2) && NBASIS_ENER_PREV<=ALTERNATE_BASIS_STRESS_B
                % See below, when ordering the basis vectors of the "expanded" basis matrix PHI
                % No special treatment
                aaa = exist(SVD_B_ws);
                if load_SVD_B == 0 || aaa == 0 ;
                    %dbstop('201')
                    disp('SVD of B-matrix')
                    [U_B{iSET_STRESS},S,V] = svd(EXP_MATRIX{iSET_STRESS},0);
                    save(SVD_B_ws,'U_B','S','V')
                else
                    load(SVD_B_ws,'U_B','S','V');
                end
            else
                % %dbstop('274')
                % U_B = BmatRIstT' ;
                % % Normalization .............................
                % nnnnorma = sqrt(sum(U_B.^2)) ;
                % for icoll = 1:length(nnnnorma)
                %     U_B(:,icoll) = U_B(:,icoll)/nnnnorma(icoll);
                % end
                DIV_FACT = sqrt(sum(EXP_MATRIX{iSET_STRESS}.^2));
                % We can impose directly 'DIV_FACT==0' because is not
                % probable find a small values in this norm!!
                NDCOMP=find(DIV_FACT<=tolDIV);  DIV_FACT(NDCOMP) = 1;
                
                U_B{iSET_STRESS} = bsxfun(@rdivide,EXP_MATRIX{iSET_STRESS},DIV_FACT) ;
            end
        case 'NO_EXPANSION'
            %dbstop('274')
            % U_B = BmatRIstT' ;
            % Normalization .............................
            % nnnnorma = sqrt(sum(U_B.^2)) ;
            % for icoll = 1:length(nnnnorma)
            %     U_B(:,icoll) = U_B(:,icoll)/nnnnorma(icoll);
            % end
            DIV_FACT = sqrt(sum(EXP_MATRIX{iSET_STRESS}.^2));
            % We can impose directly 'DIV_FACT==0' because is not
            % probable find a small values in this norm!!
            NDCOMP=find(DIV_FACT<=tolDIV);  DIV_FACT(NDCOMP) = 1;
            
            U_B{iSET_STRESS} = bsxfun(@rdivide,EXP_MATRIX{iSET_STRESS},DIV_FACT) ;
            
            % Normalization is introduced to avoid large condition numbers
            % We attempt to solve: Bmat_red'*(I-Rgappy)*SIGMA_RED
        case 'SIMPLIFIED_METHOD'
            % % We attempt to solve: Bmat_red'*SIGMA_RED  (5-Sept-2012, JAHO)
            % U_B = BmatRIstT' ;
            % % Normalization ..
            % nnnnorma = sqrt(sum(U_B.^2)) ;
            % for icoll = 1:length(nnnnorma)
            %     U_B(:,icoll) = U_B(:,icoll)/nnnnorma(icoll);
            % end
            DIV_FACT = sqrt(sum(EXP_MATRIX{iSET_STRESS}.^2));
            % We can impose directly 'DIV_FACT==0' because is not
            % probable find a small values in this norm!!
            NDCOMP=find(DIV_FACT<=tolDIV);  DIV_FACT(NDCOMP) = 1;
            
            U_B{iSET_STRESS} = bsxfun(@rdivide,EXP_MATRIX{iSET_STRESS},DIV_FACT) ;
            % JAHO, this strategy has proved unefficient; delete it when
            % "cleaning" the code
        case 'NO_NORMALIZATION'
            U_B{iSET_STRESS} = EXP_MATRIX{iSET_STRESS};
        otherwise
            error('BASIS_EXPANTION method not implemented')
    end
    
end

% Energy components on Gauss Points in the domain
TOT_ENER1 = repmat(1:e_VG_ME.nElem,npg,1);
TOT_ENER2 = repmat([1:npg]',e_VG_ME.nElem,1);

TOT_ENER = [TOT_ENER1(:) TOT_ENER2(:)];

% Numero de componentes total de tension en la microcelda
NumCOMP_ENER = [1:nGTOT]';

% ====================================
% LOOP OVER THE DOMAINS - REGULAR PART 
% ====================================
numat_B        = cell(length(SetCOMP),1);
npoints_energy = cell(length(SetCOMP),1);
NumCOMP_DOM    = cell(length(SetCOMP),1);
nDOF_ENER_SC   = cell(length(SetCOMP),1);
PHI            = cell(length(SetCOMP),1);
nENER_SET      = cell(length(SetCOMP),1);
nCOMP_TENS     = cell(length(SetCOMP),1);

CUB_DATA_GEN   = cell(length(SetCOMP),1);

INDICES_ORIGINAL_SETS = cell(length(SetCOMP),1);

nDOFEnergy_SC_TOT_SET = false(nGTOT,length(SetCOMP)); 

%nDOFsigma_SC_TOT = false(nGTOT*ntens,1);
%nmodes_disp_INIT = nmodes_disp;

% Loop over the sets of stresses at the microscale
for iSetCOMP=1:length(SetCOMP)

    switch iSetCOMP
        case 1 % REGULAR (CONTINUOUS) DOMAIN
            
            NAME = 'REGULAR_DOMAIN';            
            
            % Basta con estabilizar solo un subdominio, no se requieren de
            % los dos subdominios estabilizados, ya que el rango de la
            % matriz al final aplicara para ambos por igual
            
            % SVD de la matriz B de cada sub_dominio            
            %U_B = [];  nmodes_disp = [];

            % SVD de la matriz B de cada sub_dominio
            % [U_B,S_B,V_B] = svd(BmatGEN(SetCOMP{iSetCOMP},:),0);            
            %[EXPMAT,S_B,V_B] = svd(U_B{iSetCOMP},0);
            %nmodes_disp = 32;
            %EXPMAT = U_B{iSetCOMP};
            %EXPMAT = [];
            
            Sing_Val=[sqrt(Sing_Val_ELAS_REG); sqrt(Sing_Val_INELAS_REG)];
            
            U = PHI_ENER_REG;
            
            %U = bsxfun(@times,U,Sing_Val');
            
            %U_TRUNC = U(:,1:NBASIS_ENER_PREV_REG);
            %U_TRUNC = 10000*U(:,1:NBASIS_ENER_PREV_REG);
            
            U_TRUNC = U(:,1:12);
            
            ModosEPS = ModoPHI_REG(COMP_STR_REG,:) ;
            
            %Number of gauss points
            NBASIS_ENER_PREV =  NBASIS_ENER_PREV_REG ;

            % Integral for cubature method            
            GLOB_U=sum(bsxfun(@times,U_TRUNC,ponder_factors{iSetCOMP}),1) ;

            % Integral of strain modes over the whole domain            
            GLOB_EPS = IntPhiTGen_REG_PREV' ;
            GLOB_EPS = GLOB_EPS(:) ;
           
            % Num of Elements in Regular domain
            ELEM = sort(ElemREG) ;
            
            nmodesU = 50;            
            
        case 2 % REGULAR DOMAIN
            
            NAME = 'DISCONTINUOUS_DOMAIN';

            % SVD de la matriz B de cada sub_dominio
            % [U_B,S_B,V_B] = svd(BmatGEN(SetCOMP{iSetCOMP},:),0);
            %[EXPMAT,S_B,V_B] = svd(U_B{iSetCOMP},0);
            %nmodes_disp = 23;
            %EXPMAT = U_B{iSetCOMP};
            %EXPMAT = [];            
            
            Sing_Val=[sqrt(Sing_Val_ELAS_DIS); sqrt(Sing_Val_INELAS_DIS)];
  
            U = PHI_ENER_DIS;
            
            %U = bsxfun(@times,U(:,length(Sing_Val)),Sing_Val');
            
            
            %U_TRUNC = U(:,1:NBASIS_ENER_PREV_DIS);
            %U_TRUNC = 10000*U(:,1:NBASIS_ENER_PREV_DIS);
            %U_TRUNC = U(:,1:15);
            U_TRUNC = U;
            
            ModosEPS = ModoPHI_DIS(COMP_STR_DIS,:) ;            
            
            %Number of gauss points
            NBASIS_ENER_PREV =  NBASIS_ENER_PREV_DIS ;
            
            % Integral for cubature method
            GLOB_U=sum(bsxfun(@times,U_TRUNC,ponder_factors{iSetCOMP}),1) ;
            
            % Integral of strain modes over the whole domain
            GLOB_EPS = IntPhiTGen_DIS_PREV' ;
            GLOB_EPS = GLOB_EPS(:) ;            
            
            % Num of Elements in Regular domain
            ELEM = sort(ElemDIS) ;
            
            nmodesU = 200;

    end
    
    EXPMAT = eye(size(U,1)) - U(:,1:NBASIS_ENER_PREV)*U(:,1:NBASIS_ENER_PREV)';                

    %switch  e_VG.BASIS_EXPANSION
    % PHI = ConsExpMatrix(ALTERNATE_BASIS_STRESS_B,NBASIS_STRESS_PREV,U_B,U);
    PHI{iSetCOMP} = ConsExpMatrix(ALTERNATE_BASIS_STRESS_B,NBASIS_ENER_PREV,EXPMAT(:,1:nmodes_disp{iSetCOMP}),U);
    %PHI{iSetCOMP} = ConsExpMatrix(ALTERNATE_BASIS_STRESS_B,NBASIS_ENER_PREV,EXPMAT(:,1:NBASIS_ENER_PREV),U);
        
    % Determining the number of sample points  (npoints_gappy)
    [npoints_gappy,PHI_AMP,STUDY_GREEDY_LOC,ndof_sigma_SC_PRE,npoints_sigma_PRE,numat_B_PRE] = ...
        DetermineNpoingGAUSS(e_VG_ME,PHI{iSetCOMP},EXPMAT,nGTOT,U,NGAUSS_GAPPY_loc);
    
    % This is the "expanded" basis   [U_B PHI_SIGMA] (denominations are somehow "misleading")
    % In latex,  \BasisSp = PHI
       
    % SELECTION OF GAUSS POINTS
    NAME_LOAD_GREEDY = ['DATA_JAHO/',file(1:end-4),'_',GREEDY_ALGOR,'_',...
        num2str(npoints_gappy),'_',NAME,'_',num2str(NBASIS_ENER_PREV),'_',num2str(size(BmatRIstT',2))];    
    
    switch SNAPSHOTS_ORIGIN
        case 'ROMI'
            NAME_LOAD_GREEDY = [NAME_LOAD_GREEDY,'.mat'];
        case 'HF'
            NAME_LOAD_GREEDY = [NAME_LOAD_GREEDY,'_HF.mat'];
    end
    
    %e_VG.NAME_LOAD_GREEDY = NAME_LOAD_GREEDY ;
    if ~exist(NAME_LOAD_GREEDY)
        LOAD_RESULT_GREEDY = 0;
    end
    
    %LOAD_RESULT_GREEDY = 1;
    % ESTE ANALISIS FUNCIONA PARA MALLADOS CONFORMADOS POR ELEMENTOS DEL MISMO
    % NUMERO DE PUNTOS DE GAUSS!!!
    
    % ACTUALIZAR  ###################
    npg = 4 ;
    % ACTUALIZAR  ###################
    
    puntos_adicionales = 0 ;
    interpolation_method = 0;
    
    NAME_SAVE_GREEDY = NAME_LOAD_GREEDY;
    if  STUDY_GREEDY_LOC  == 1  ;
        if  LOAD_RESULT_GREEDY == 0 ;
            %GREEDY_ALGOR = e_VG.GREEDY_ALGOR ;
            
            if isempty(SET_GPS{iSetCOMP})              
                
                if (interpolation_method == 1)
                    
                    [~,~,cond_NUMB_WITHOUT,ndof_ener_SC,cond_NUMB_WITH,WORST_COMB,INDICES_ORIGINAL_SETS{iSetCOMP},...
                        INDICES_ENER_RECONST,INDICES_GAUSS_RECONST,SELECT_ELEMENT_STRESS]...
                        = greedy_algorithms(PHI_AMP,GREEDY_ALGOR,ntens,npg,'ng',ntens,'puntos_adicionales',...
                        puntos_adicionales,'CHATUR_NONREPEATED',CHATUR_NONREPEATED,'INDICES_GIVEN_LOC',INDICES_GIVEN_LOC,...
                        'PARAMETER_GREEDY',PARAMETER_GREEDY,'NMODES_stress',NBASIS_ENER_PREV,'NMODE_disp',nmodes_disp,...
                        'U_left',U,'S_singular',SingleV_ENER,'BASIS_EXPANSION',BASIS_EXPANSION);
                else
                                        
                    EPS_MAT = cell2mat(arrayfun(@(x)reshape(ModosEPS(:,x),ntens,length(ponder_factors{iSetCOMP}),[]),1:nmodes_disp{iSetCOMP},'UniformOutput',false)');
                    
                    compEPS = true(size(EPS_MAT,1),1); compEPS(3:ntens:end,1)=false;
                    EPS_MAT=EPS_MAT(compEPS,:); % without out-of-plane component (plane strain)
                    
                    %G = [U_TRUNC';EPS_MAT] ;
                    %b = [GLOB_U'; GLOB_EPS] ;
                    
                    %G = [U_TRUNC';EPS_MAT; ones(1,length(ponder_factors{iSetCOMP}))] ;
                    %b = [GLOB_U'; GLOB_EPS; sum(ponder_factors{iSetCOMP})] ;
                    
                    G = [U_TRUNC'; ones(1,length(ponder_factors{iSetCOMP}))] ;
                    b = [GLOB_U'; sum(ponder_factors{iSetCOMP})] ;                                        
                    
                    %G = [U_TRUNC';ones(1,length(ponder_factors{iSetCOMP}))] ;
                    %b = [GLOB_U'; sum(ponder_factors{iSetCOMP})] ;
                    
                    Jmin = G ;  bMIN = b ;
                    
                    DATA_CUBATURE.VOLUME_CONSTRAINT = 0;
                    DATA_CUBATURE.RELOAD_POINTS = 0;
                    DATA_CUBATURE.PLOT_WEIGHTS = 1;
                    DATA_CUBATURE.INTENSITY_ELASTIC_FACTOR_FORCE = 10;
                    DATA_CUBATURE.INTENSITY_ELASTIC_FACTOR_DISPL = 10;
                    DATA_CUBATURE.WEIGHTED_MODES= 1;
                    DATA_CUBATURE.INCLUDE_SVDISPL = 1;
                    DATA_CUBATURE.REMOVE_ZERO_WEIGHTS = 1;
                    DATA_CUBATURE.RELOAD_Jmin = 0;
                    DATA_CUBATURE.value_RESIDini = 0;
                    
                    setPoints = []; % npoints = [];
                    BOUND_WEIGHTS_m_mg = 0;
                    
                    %nmodesUelas = 3;
                    %TYPE_SVD = 'ELASTIC_INELASTIC_SVD' ;
                                       
                    % numero total de modos (energia mas deformacion)
                    %nmodesU = NBASIS_ENER_PREV + nmodes_disp{iSetCOMP} ;
                    
                    
                    %%%%%%%%nmodesU = NBASIS_ENER_PREV + 1 ;
                    
                    %[setPoints,weights,TIME_JMIN,TIME_POINTS] =  PointsCubature_2DFE(nmodesU,npoints,setPoints,...
                    %    BOUND_WEIGHTS_m_mg,DATA_CUBATURE,length(ponder_factors{iSetCOMP}),ponder_factors{iSetCOMP},G,b);

                    UNSTRESSED_POINTS = [] ;
                    
                    %[setPoints,weights,norm_resid,TIME_POINTS,setPointsALL,Jmin,bMIN] = ...
                    %    GreedyCubature_2DFE(bMIN,Jmin,ngausG,nmodesU,BOUND_WEIGHTS_m_mg,DATA_CUBATURE,UNSTRESSED_POINTS,...
                    %    weigE,setPoints,TIME_POINTS);
                    
                    [setPoints,weights,norm_resid,TIME_POINTS,setPointsALL,Jmin,bMIN] = GreedyCubature_2DFE(bMIN,Jmin,...
                        nmodesU,nmodesU,BOUND_WEIGHTS_m_mg,DATA_CUBATURE,UNSTRESSED_POINTS,...
                        ponder_factors{iSetCOMP},setPoints);                    
                    
                    % Imposing GP's in an Ad-hoc way
                    %ndof_ener_SC = IndicesGPsWCOMP_CUB(setPoints,e_DatSet.e_DatMat.e_DatSet(iSetCOMP).m_NumElem,npg,ntens) ;
                    
                    %ndof_ener_SC = IndicesGPsWCOMP_CUB(setPoints,ELEM,npg,ntens) ;
                    %ndof_ener_SC = sort(setPoints,ELEM,npg,ntens) ;
                    
                    % sorting integration points and weights
                    CUB_DATA = [setPoints' weights] ;
                    CUB_DATA = sortrows(CUB_DATA,1) ;
                    
                    ndof_ener_SC = CUB_DATA(:,1) ;
                    
                    % Storing points & weights for cubature method
                    CUB_DATA_GEN{iSetCOMP} = CUB_DATA ;
                    
                end


            else
                % Imposing GP's in an Ad-hoc way
                ndof_ener_SC = IndicesGPsWCOMP(SET_GPS{iSetCOMP},e_DatSet.e_DatMat.e_DatSet(Set_BANDAS).m_NumElem,npg,ntens) ;
            end
            
            NumCOMP_DOM{iSetCOMP}  = NumCOMP_ENER(SetCOMP{iSetCOMP});
            nDOF_ENER_SC{iSetCOMP} = ndof_ener_SC;
            
            % COMPONENTES DE TENSION SELECCIONADAS POR SET
            nDOFEnergy_SC_TOT_SET(NumCOMP_DOM{iSetCOMP}(ndof_ener_SC),iSetCOMP) = true;
            TMP = TOT_ENER(nDOFEnergy_SC_TOT_SET(:,iSetCOMP),:);
            %npoints_energy{iSetCOMP} = TMP(1:ntens:end,:) ;
            npoints_energy{iSetCOMP} = TMP;
            
            % ID - PHI MATRICES (the same ndof_ener_SC because of the use of energy per GP in 2nd reduction
            numat_B{iSetCOMP} = (npoints_energy{iSetCOMP}(:,1)-1)*npg+npoints_energy{iSetCOMP}(:,2) ;
            
            % For computing offline reconstruction operator (Energy Version)
            TMP=arrayfun(@(x)[(ndof_ener_SC(x)-1)*ntens+1:ndof_ener_SC(x)*ntens],1:length(ndof_ener_SC),'UniformOutput',false);
            
            %strain_comp{iSetCOMP} = sort(additional_rows(ntens,numat_B{iSetCOMP},'double'));
            strain_comp{iSetCOMP} = sort(cell2mat(TMP));
            
            %Sets of integration Gauss points
            if iSetCOMP == 1
                nENER_SET{iSetCOMP} = 1:length(npoints_energy{iSetCOMP});
                %nCOMP_ENER{iSetCOMP} = 1:length(npoints_energy{iSetCOMP});
                nCOMP_TENS{iSetCOMP} = 1:length(strain_comp{iSetCOMP});
            else
                nENER_SET{iSetCOMP} = [1:length(npoints_energy{iSetCOMP})] + repmat(length(npoints_energy{iSetCOMP-1}),1,length(npoints_energy{iSetCOMP}));
                %nCOMP_ENER{iSetCOMP} = [1:length(npoints_energy{iSetCOMP})] + repmat(length(npoints_energy{iSetCOMP-1}),1,length(npoints_energy{iSetCOMP}));
                nCOMP_TENS{iSetCOMP} = [1:length(strain_comp{iSetCOMP})] + repmat(length(strain_comp{iSetCOMP-1}),1,length(strain_comp{iSetCOMP}));
            end
            
            npoints_energy_SET = npoints_energy{iSetCOMP};
            numat_B_SET = numat_B{iSetCOMP};
            
            save(NAME_SAVE_GREEDY,'npoints_energy_SET','numat_B_SET');
            label_leg = NAME_SAVE_GREEDY  ;

        else
            load(NAME_LOAD_GREEDY); % Load previously computed set of indices
            % e_VG.GREEDY_ALGOR = GREEDY_ALGOR ;
        end
        % To ensure that they are "sorted" ( it is somehow redundant ...)
        % dbstop('261')
    end
       
end

% COMPONENTES DE TENSION SELECCIONADAS TOTALES
nDOF_ENER_SC_TOT=any(nDOFEnergy_SC_TOT_SET,2);

%setElem_red = npoints_sigma;
setElem_red = cell2mat(npoints_energy);
numat_B = cell2mat(numat_B);

% FLAG IN SELECTED GAUSS POINTS (FOR POSTRPROCES PURPOSES)
GaussFLAG=zeros(npg,e_VG_ME.nElem);
GaussFLAG(numat_B)=1;
GaussFLAG = [(1:e_VG_ME.nElem);GaussFLAG];

%ROM_II.nmodesU = nmodesU ;                 % indices of selected rows
ROM_II.nModesEPS_TOT = e_VG.nModesEPS_TOT ; % indices of selected rows

% % % ROM_II.ndof_sigma_SC = ndof_sigma_SC ;   % indices of selected rows
%ROM_II.npoints_sigma = npoints_sigma ;   % indicies of gauss points associated to the rows e_VG.ndof_sigma_SC
% % % ROM_II.npoints_sigma =  setElem_red;      % indices of gauss points associated to the rows e_VG.ndof_sigma_SC

% % % ROM_II.INDICES_ORIGINAL = INDICES_ORIGINAL;

% % % ROM_II.SELECTED_SAMPLE_POINTS = numat_B ;

ROM_II.GaussFLAG = GaussFLAG';

switch BASIS_EXPANSION
    case 'SIMPLIFIED_METHOD'
        
        % JAHO, this strategy has proved unefficient; delete it when
        % "cleaning" the code
        disp(['Number of selected points (EQUILIBRIUM POINTS) = ',num2str(length(numat_B))])
        disp(['SELECTED POINTS = ',num2str(numat_B)])
        disp(['Number of selected points (STRESS POINTS) = ',num2str(length(INDICES_GAUSS_RECONST))])
        try
            disp(['SELECTED POINTS = ',num2str(INDICES_GAUSS_RECONST)])
        catch
            disp(['SELECTED POINTS = ',num2str(INDICES_GAUSS_RECONST')])
        end
    otherwise
        try
            disp(['SELECTED POINTS = ',num2str(numat_B)])
        catch
            disp(['SELECTED POINTS = ',num2str(numat_B')])
        end
        disp(['Number of selected points = ',num2str(length(numat_B))])
end

%%%%
% APPROXIMATION ERROR ANALYSIS
NFIG_ERROR_1 = 70 ;

% REVISAR!!!!!!!!!!!!!!!!
%e_VG.ponder_factors = ponder_factors ;
% % switch COMPARE_WITH_SNAPSHOTS
% %     case 'YES'
% %         switch APPROXIMATION_ERROR_ANALYSIS
% %             case 'YES'
% %                 [Hplot1LOC,legendLOC]=study_error1(PHI,SIGMA_SNAP_SC,ndof_sigma_SC,e_VG,TYPE_NORM);
% %         end
% % end

% This is the "gappy" basis matrix (ndof_sigma_SC) corresponds to the selected indices
% In latex, G = \BasisSg, and G2 = \Mgappy
%dbstop('424')
switch  BASIS_EXPANSION
    case {'NO_EXPANSION','NO_NORMALIZATION'}
        for iSetCOMP = 1:length(SetCOMP)
            switch iSetCOMP
                case 1 % SINGULAR DOMAIN
                    NBASIS_ENER_PREV = NBASIS_ENER_PREV_REG ;
                case 2 % REGULAR DOMAIN
                    NBASIS_ENER_PREV = NBASIS_ENER_PREV_DIS ;
            end
            PHI_G=PHI{iSetCOMP}(nDOF_ENER_SC{iSetCOMP},1:NBASIS_ENER_PREV);
            M=PHI_G'*PHI_G;
            % CONDITION NUMBRE
            condM = cond(M) ;
            %condM = condM;
            fprintf('Condition number for M=PHI_G_T*PHI_G for Set No %i : %d \n',iSetCOMP,condM)
            fprintf('\n')
            %fprintf('\n')
            if  cond(M) > 100
                warning('condition number of M very high' );
            end
        end
    otherwise
        PHI_G=PHI(ndof_ener_SC,:);
        M=PHI_G'*PHI_G;
        
        % CONDITION NUMBER
        condM = cond(M) ;
        %condM = condM;
        fprintf('Condition number for M=PHI_G_T*PHI_G: %d \n',condM)
        fprintf('\n')
        fprintf('\n')
        if  cond(M) > 100
            warning('condition number of M very high' );
        end
        
end

% MATRICES B DE LOS PG'S SELECCIONADOS
% The ROM B-matrix corresponding to each sample point is stored in
% BmatRIst_W_red (as a cell array)
ROM_II.BmatRIst_W_red = BmatRIst_W(numat_B(:));

AENER      = cell(length(PHI),1);
B_II_sigma = cell(length(PHI),1);

% FOR CONVEXIFICATION PURPOSES
%B_II_diff_c_q = cell(length(PHI),1);

Q          = cell(length(SetCOMP),1);

switch  BASIS_EXPANSION
    case {'SIMPLIFIED_METHOD'}
        % JAHO, this strategy has proved unefficient; delete it when
        % "cleaning" the code
        PHI_G=PHI(INDICES_ENER_RECONST,1:NBASIS_ENER_PREV);
        M=PHI_G'*PHI_G;
        AENER=  PHI(:,1:NBASIS_ENER_PREV)*(M\PHI_G') ;
        
        % REduced Bmat operator --->  Bgappy'
        BmatRIst = BmatRIstT' ;
        Bgappy = BmatRIst(ndof_ener_SC,:) ;
        B_II_sigma= Bgappy'  ;
        
    case {'NO_EXPANSION','NO_NORMALIZATION'}
        %PHI_G=PHI(ndof_sigma_SC,1:NBASIS_STRESS_PREV);
        %M=PHI_G'*PHI_G;
        %ASIGMA=  PHI(:,1:NBASIS_STRESS_PREV)*(M\PHI_G') ;
        
        % In this case, the U_B modes are not considered
        % for reconstruction purposes (on the grounds that the associated fourier coefficients
        % becomes negligible )
        %
        % REduced Bmat operator --->  Bgappy'*(I-Rgappy)
        BmatRIst = BmatRIstT' ;        
        for iSetCOMP = 1:length(SetCOMP)
            switch iSetCOMP
                case 1 % REGULAR DOMAIN
                    NBASIS_ENER_PREV = NBASIS_ENER_PREV_REG ;
                    %nCOL = 1:length(nDOF_SIGMA_SC{1});
                    nCOL = 1:size(ModoPHI_REG,2);
                    
                    nMODE = size(ModoPHI_REG,2);
                    
                case 2 % DISCONTINUOUS DOMAIN (BANDS)
                    NBASIS_ENER_PREV = NBASIS_ENER_PREV_DIS ;
                    %nCOL = length(nDOF_SIGMA_SC{1})+1:SIGMA_COMP;
                    nCOL = size(ModoPHI_REG,2)+1:size(ModoPHI_REG,2)+size(ModoPHI_DIS,2);
                    
                    nMODE = size(ModoPHI_DIS,2);
            end
            
            PHI_G=PHI{iSetCOMP}(nDOF_ENER_SC{iSetCOMP},:);
            M=PHI_G'*PHI_G;
            AENER{iSetCOMP} = PHI{iSetCOMP}*(M\PHI_G') ;
            
            %%% Bgappy = BmatRIst(ndof_sigma_SC,:) ;
            %Bgappy = BmatRIst(nDOFsigma_SC_TOT_SET(:,iSetCOMP),:) ;
            Bgappy = BmatRIst(nDOFEnergy_SC_TOT_SET(:,iSetCOMP),nCOL) ;
            
            %%%% ReconsGappy = ASIGMA(ndof_sigma_SC,:);
            ReconsGappy = AENER{iSetCOMP}(nDOF_ENER_SC{iSetCOMP},:);
            
            Ident = eye(size(ReconsGappy)) ;
            
            %B_II_sigma{iSetCOMP}= Bgappy'*(Ident-ReconsGappy) ;

            % Be aware when upgrading the code at this part of the code
            if interpolation_method==1
                
                Q_SET = zeros(1,size(AENER{iSetCOMP},2));
                for ipoint = 1:length(ponder_factors{iSetCOMP})
                    Q_SET = Q_SET + AENER{iSetCOMP}(ipoint,:)*ponder_factors{iSetCOMP}(ipoint);
                end
                
                %Q_SET = ones(1,size(AENER{iSetCOMP},2));
                % Pass from scalar integration weigth to an array of new
                % quadrature weigth
                Q_SET_TMP = repmat(Q_SET,ntens,1);
                Q_SET2 = Q_SET_TMP(:) ;
                W2 = diag(Q_SET2);
                
                % PARA PESOS PARTE NO ADMISIBLE
                %Q_SET_TMP2 = sum(Q_SET) ;
                %W3 = Q_SET_TMP2*eye(nMODE) ;
                
            else
                
                % Reduced integration method
                Q_SET = CUB_DATA_GEN{iSetCOMP}(:,2) ;
                Q_SET_TMP = repmat(Q_SET',ntens,1);
                Q_SET2 = Q_SET_TMP(:) ;
                W2 = diag(Q_SET2) ;
                
            end
            
            %Q2 = Q_SET/e_DatSet.e_DatMat.omegaMicro;
            %ROM_II.Q_oper = Q/(sum(sum(bsxfun(@times,e_VG_ME.wg,m_DetJT))));
            
            % If is possible to use the bsxfun function on the assessment of R, we can
            % use the f_HomogArea function for the assessment of Q_oper!!!!
            
            %B_II_sigma{iSetCOMP}= BmatRIst(SetCOMP{iSetCOMP},nCOL)'*AENER{iSetCOMP} ;
            
            % In this case, the PHI matrix used to compute the internal forces is
            % the reduced one instad the total one (2nd reduction based on energy)
            
            % for 2nd reduction in stresses
            %B_II_sigma{iSetCOMP}= BmatRIst(SetCOMP{iSetCOMP}(strain_comp{iSetCOMP}),nCOL)'*W2 ;
            
            if iSetCOMP == 1 % REGULAR DOMAIN
                
                % PHI_REG or PHI_DIS without weigth!!!!  the weights are stored in W2 (for 2nd reduction in energy)
                B_II_sigma{iSetCOMP}= ModoPHI_REG(COMP_STR_REG(strain_comp{iSetCOMP}),:)'*W2 ;
                
                %B_II_diff_c_q{iSetCOMP} = e_VG_ME.K_PARAM*W3 ;
                
                %WCELL = arrayfun(@(x)diag(repmat(ponder_factors{iSetCOMP}(x),1,ntens)),1:length(ponder_factors{iSetCOMP}),'UniformOutput',false);
                WCELL = arrayfun(@(x)diag(Q_SET_TMP(:,x)),1:size(Q_SET_TMP,2),'UniformOutput',false);
                WMAT = cell2mat(WCELL);
                
                Q{iSetCOMP} = WMAT./e_DatSet.e_DatMat.omegaMicro ; 
                
                %IntPhiTGen_REG = WMAT*ModoPHI_REG(COMP_STR_REG(strain_comp{iSetCOMP}),:);
                %IntPhiTGen_REG = e_VG_ME.E_MATRIX*IntPhiTGen_REG([1 2 4],:)';
                
                % For computing homogenized macro-stress
                %ROM_II.nCOMP_TENS{iSetCOMP}=COMP_STR_REG(strain_comp{iSetCOMP})';
                
            elseif iSetCOMP == 2 % DISCOTINUOUS DOMAIN
                
                % PHI_REG or PHI_DIS without weigth!!!!  the weights are stored in W2 (for 2nd reduction in energy)
                B_II_sigma{iSetCOMP}= ModoPHI_DIS(COMP_STR_DIS(strain_comp{iSetCOMP}),:)'*W2 ;                
                
                %WCELL = arrayfun(@(x)diag(repmat(ponder_factors{iSetCOMP}(x),1,ntens)),1:length(ponder_factors{iSetCOMP}),'UniformOutput',false);
                WCELL = arrayfun(@(x)diag(Q_SET_TMP(:,x)),1:size(Q_SET_TMP,2),'UniformOutput',false);
                WMAT = cell2mat(WCELL);
                
                Q{iSetCOMP} = WMAT./e_DatSet.e_DatMat.omegaMicro ;   
                
                %B_II_diff_c_q{iSetCOMP} = e_VG_ME.K_PARAM*W3 ;                
                
                %IntPhiTGen_DIS = WMAT*ModoPHI_DIS(COMP_STR_DIS(strain_comp{iSetCOMP}),:);
                %IntPhiTGen_DIS = e_VG_ME.E_MATRIX*IntPhiTGen_DIS([1 2 4],:)';
                
                % For computing homogenized macro-stress
                %ROM_II.nCOMP_TENS{iSetCOMP}=COMP_STR_DIS(strain_comp{iSetCOMP})';
            end
            
        end
        
        %    N_gappy = B_II_sigma*Bgappy ; % N_gappy = Bgappy'*(Ident - ReconsGappy)*Bgappy
        % The "NO_EXPANSION" approach is predicated on the
        % existence of the inverse of N_gappy. A poorly conditioned Ngappy
        % is likely to produce ill-conditioned tangent stiffness matrices.
    otherwise
        PHI_G=PHI(ndof_ener_SC,1:NBASIS_ENER_PREV);
        M=PHI_G'*PHI_G;
        AENER=  PHI(:,1:NBASIS_ENER_PREV)*(M\PHI_G') ;
        B_II_sigma= BmatRIstT*AENER ;
end

ROM_II.nGTOT = nGTOT ;
ROM_II.nTENS_SET = nENER_SET;

% For computing homogenized stress tensor
ROM_II.nCOMP_TENS = nCOMP_TENS;

%ROM_II.IntPhiTGen = IntPhiTGen;
ROM_II.IntPhiTGen_REG = IntPhiTGen_REG;
ROM_II.IntPhiTGen_DIS = IntPhiTGen_DIS;

% Convexification matrices
% % % ROM_II.K_CONVEX_REG = K_CONVEX_REG;
% % % ROM_II.K_CONVEX_DIS = K_CONVEX_DIS;

% % % ROM_II.BmatRIstT_REG = BmatRIstT(1:nModesEPS_bar,SetCOMP{1}) ;
% % % ROM_II.BmatRIstT_DIS = BmatRIstT(nModesEPS_bar+1:nModesEPS_bar+nModesEPS_bar2,SetCOMP{2}) ;

ROM_II.NPG_DOM = [length(ponder_factors{1}) length(ponder_factors{2})] ;

% % % ROM_II.ASIGMA = AENER;

% % % ROM_II.B_II_sigma = B_II_sigma;
% % % ROM_II.PHI_ASIGMA = BmatRIstT*ASIGMA ;
%ROM_II.Q_ROMII = BmatRIstT*ASIGMA;

%%% JAHO, NEW RECONSTRUCTION OPERATOR  : THIS STRATEGY HAS PROVED UNEFFECTIVE!!!!
% -------------------------------------------------------------------------------
% 18-May-2012: So far, this operator will only affect postprocess results. Instead of defining PHI as an expanded
% matrix containing the SVD stress basis vectors themselves and the U_B
% vectors, we discard U_B (on the grounds that it plays no role in re-constructing or interpolating
% the stress field) and replace it by the next nU modes . The new basis
% matrix will be called  PHI_RECONS
%dbstop('633')
% PHI_RECONS  =[U(:,1:size(PHI,2)) ];
%
% % Now we select the rows corresponding to indices ndof_sigma_SC
% PHI_G_RECONST=PHI_RECONS(ndof_sigma_SC,:);
% % and calculate M_RECONS
% M_RECONS=PHI_G_RECONST'*PHI_G_RECONST;
% % on which the existence of the reconstruction operator is predicated
% %
% condM_RECONS = cond(M_RECONS) ;
% disp(['Conditioning number for  PHI_RECONST''*PHI_RECONST =',num2str(condM_RECONS)])
%
% if condM_RECONS > 1e10  % Threshold for M_RECONS to be deemed as ill-conditioned
%     warning('Too large a condition number for PHI_RECONST')
%     pause(2)
% end
% ASIGMA_RECONS=  PHI_RECONS*(M_RECONS\PHI_G_RECONST') ;
% e_VG.ASIGMA_RECONS = ASIGMA_RECONS

% NEW STRATEGY
if PARTIAL_RECONSTR == 1
    switch GREEDY_ALGOR
        case 'HIERARCHICAL_cycle'
        otherwise
            error(['Option PARTIAL_RECONSTR = 1 is, so far, only valid for GREEDY_ALGOR_GLO_CHOICE =HIERARCHICAL_cycle'])            
    end
    
    %   Only the first "imod_stress" gauss points are employed in the reconstruction
    %   of the stress field.
    
    % STEP 1) BASIS MATRIX. ONLY STRES BASES ARE CONSIDERED FOR RECONSTRUCTION PURPOSES
    PHI_RECONS  =[U(:,1:NBASIS_ENER_PREV) ] ;
    % STEP 2) Now we select the   rows corresponding to the first
    % "imod_stress*e_VG.ntens"  indices
    ndof_sigma_SC_stresses = ndof_ener_SC(1:NBASIS_ENER_PREV*ntens) ;
    PHI_G_RECONST=PHI_RECONS(ndof_sigma_SC_stresses,:);
    % STEP 3) and calculate M_RECONS
    M_RECONS=PHI_G_RECONST'*PHI_G_RECONST;
    % STEP 4) The reconstruction operator is calculated as
    ASIGMA_RECONS=  PHI_RECONS*(M_RECONS\PHI_G_RECONST') ;
    %
    % % % ROM_II.ASIGMA_RECONS = ASIGMA_RECONS;
end

%ndof_sigma_red = length(ndof_sigma_SC);

% % REDUCED CONSTITUTIVE OPERATOR
% switch BASIS_EXPANSION
%     case {'DIRECT_U','SIMPLIFIED_METHOD','NO_EXPANSION','NO_NORMALIZATION'}
%         
%         % NOT IMPLEMENTED YET
%         % -------------------
%         %   ROM_II.R_oper = []; %sum(R,3);
%         %   ROM_II.Q_oper =[] ;
%         % JAHO
%         % otherwise
%         % **********************
%         % Constitutive operators
%         % **********************
%         % In this case, the matrix BmatRIstT now contains the ponder factor
%         %%%% Q = zeros(e_VG_ME.ntens,ndof_sigma_red);
%         Q = cell(length(SetCOMP),1);
%         
%         for iSetCOMP=1:length(SetCOMP)
%             % Be aware when upgrading the code at this point of the
%             % code
%             %bsxfun(@times,repmat(ASIGMA,)
%             Q_SET = zeros(ntens,size(AENER{iSetCOMP},2));
%             for ipoint = 1:length(ponder_factors{iSetCOMP})
%                  Q_SET = Q_SET + AENER{iSetCOMP}((ipoint-1)*ntens+1:ipoint*ntens,:)*ponder_factors{iSetCOMP}(ipoint);
%             end
%             Q{iSetCOMP} = Q_SET/e_DatSet.e_DatMat.omegaMicro;
%             %ROM_II.Q_oper = Q/(sum(sum(bsxfun(@times,e_VG_ME.wg,m_DetJT))));
%             
%             % If is possible to use the bsxfun function on the assessment of R, we can
%             % use the f_HomogArea function for the assessment of Q_oper!!!!
%         end
% 
%         ROM_II.Q_oper = Q ;
%         
%         %R = BmatRIstT*ASIGMA;        
%         R = B_II_sigma;  % JAHO: I am not sure if this modification is correct .... 28-May-2012.
%         ROM_II.R_oper = R; %sum(R,3);        
% end

R = B_II_sigma;  % JAHO: I am not sure if this modification is correct .... 28-May-2012.
ROM_II.Q_oper = Q;
ROM_II.R_oper = R;

% When weights are incorporated into the STRESS VECTOR,
% we have to multiply the gappy stresses returned by the
% constitutive equation subroutine by the corresponding weights.
% In the following, the "reduced" weighing vector is defined

ROM_II.B_II_sigma_W = B_II_sigma;
%ROM_II.B_II_diff_c_q = B_II_diff_c_q ;

%ROM_II.W_STRESS = (sum(ponder_factors{1})+sum(ponder_factors{2}))/e_DatSet.e_DatMat.omegaMicro;
% % % ROM_II.W_STRESS = (sum(ponder_factors{1})+sum(ponder_factors{2}))/e_DatSet.e_DatMat.omegaMicro;

% Matrix for recovering displacements from strains
ROM_II.ValRM = ValRM;
ROM_II.MAT_WEAK = MAT_WEAK;

%e_VG.ASIGMA_W =ASIGMA;
% dbstop('310')

switch WHERE_WEIGHTS
    case 'WinB'
        % With no weights in S (they are already in B_II_sigma)
        % e_VG.B_II_sigma_W = B_II_sigma;
        
    case 'WinS'
        % dbstop('350')
        warning(['Option ',WHERE_WEIGHTS, 'contains some sort of error, since it furnishes non-sensical results']);
        ppff=  repmat(ponder_factors(setElem_red),[1,ntens]);
        ppff = reshape(ppff',[ntens*npg*length(setElem_red),1]);
        ROM_II.ponder_factors_RED = ppff ;
        diag_ponder = diag(ppff);
        % e_VG.B_II_sigma_W = bsxfun(@times,B_II_sigma,ppff');  % THIS ONE INCLUDE THE WEIGHTS !!!
        ROM_II.B_II_sigma_W =  B_II_sigma*diag_ponder; % = diag(ppff);
end

% Fix the set of the selected GP's
ROM_ELEM_SET = zeros(length(setElem_red),1) ;
ROM_ksb = zeros(length(setElem_red),1) ;
for iSet=1:nSet
    %SET_TMP = e_DatSet(ME_Set).e_DatMat.e_DatSet(iSet).m_NumElem;
    [iD,~] = ismember(setElem_red(:,1),e_DatSet.e_DatMat.e_DatSet(iSet).m_NumElem) ;
    ROM_ELEM_SET(iD) = iSet ;
    if e_DatSet.e_DatMat.e_DatSet(iSet).e_DatElem.eltype==31 && (~isempty(iD))
        element_set=setElem_red(iD,1) ;
        AUX_ksb=[];
        for iKsb=1:sum(iD)
            index=find(e_DatSet.e_DatMat.e_DatSet(iSet).m_NumElem==element_set(iKsb));
            AUX_ksb(iKsb,1)=e_DatSet.e_DatMat.e_DatSet(iSet).e_DatElem.ksb(index);
        end
        ROM_ksb(iD)=AUX_ksb ;
    end
end

diffEL=find(ROM_ksb==0);
ROM_ksb(diffEL)=1;

ROM_II.ROM_ELEM_SET=ROM_ELEM_SET ;
ROM_II.ROM_ksb=ROM_ksb ;

% e_DatSet for ROM SET OF ELEMENTS
% ********************************
e_DatSet_ROM.nElem   = length(setElem_red);
e_DatSet_ROM.sihvare = e_DatSet.e_DatMat.e_DatSet(1).sihvare/npg ;
e_DatSet_ROM.siavare = e_DatSet.e_DatMat.e_DatSet(1).siavare/npg ;
e_DatSet_ROM.sitvare = e_DatSet.e_DatMat.e_DatSet(1).sitvare/npg ;
e_DatSet_ROM.dofpe   = e_DatSet.e_DatMat.e_DatSet(1).e_DatElem.dofpe;

e_DatSet_ROM.siavarpg = e_DatSet.e_DatMat.e_DatSet(1).e_DatMat.siavarpg;
e_DatSet_ROM.sihvarpg = e_DatSet.e_DatMat.e_DatSet(1).e_DatMat.sihvarpg;

e_DatSet_ROM.e_DatElem.nVarHistElem = e_DatSet.e_DatMat.e_DatSet(1).e_DatElem.nVarHistElem ;
e_DatSet_ROM.e_DatElem.nVarAuxElem  = e_DatSet.e_DatMat.e_DatSet(1).e_DatElem.nVarAuxElem ;
e_DatSet_ROM.e_DatElem.eltype       = e_DatSet.e_DatMat.e_DatSet(1).e_DatElem.eltype ;
e_DatSet_ROM.e_DatElem.pointersVAE.p_kSD = e_DatSet.e_DatMat.e_DatSet(1).e_DatElem.pointersVAE.p_kSD;
e_DatSet_ROM.e_DatElem.ksb          = ROM_ksb;

e_DatSet_ROM.e_DatMat.conshyp  = e_DatSet.e_DatMat.e_DatSet(1).e_DatMat.conshyp ;
e_DatSet_ROM.e_DatMat.esImplex = e_DatSet.e_DatMat.e_DatSet(1).e_DatMat.esImplex ;
e_DatSet_ROM.e_DatMat.tit      = e_DatSet.e_DatMat.e_DatSet(1).e_DatMat.tit ;

ROM_II.e_DatSet_ROM = e_DatSet_ROM ;

% THIS STRUCTURE (ROM_II) IS ASSIGNED FROM NOW ON OUTSIDE THESE ROUTINE
% ===========================================================
% % % % Adding the hrom data in general multiscale data structure.
% % % e_DatSet.e_DatMat.ROM_II = ROM_II ;
% ===========================================================

% Modify variables in e_DatSet for ROMII purposes
% e_DatSet(ME_Set).e_DatMat.e_VG.nSet = 1 ; % Corresponding to one (unique) set of elements in the small scale.

% JAHO, 8-Jun-2012
% RECONSTRUCTION OPERATOR FOR INTERNAL VARIABLES
disp('Reconstruction operator for Internal Variables ')
% STEP 1: Retrieving internal variables bases
% -------------------------------------------
RECONST_OPER_INTV_GLO = [] ;
if  ~isempty(LEVEL_TRUN_INTERN)
    load(BASIS_DATA,'PHI_intern_GLO','DATA_VARPLOT')
    for iintvar = 1:length(PHI_intern_GLO)
        PHI_RED = PHI_intern_GLO{iintvar}(numat_B,1:LEVEL_TRUN_INTERN) ;
        M_INT_RED = PHI_RED'*PHI_RED ;
        condMM = cond(M_INT_RED) ;
        disp('-------------------------------------------------------------')
        disp(['RECONSTR. of INTERNAL VARIABLE -->',DATA_VARPLOT(iintvar).NAME])
        disp(['Condition number of Mgappy = ',num2str(condMM)]) ;
        disp('-------------------------------------------------------------')
        COEFF = (M_INT_RED\PHI_RED') ;
        RECONST_OPER_INTV_GLO{iintvar} = PHI_intern_GLO{iintvar}(:,1:LEVEL_TRUN_INTERN)*COEFF ;
    end
end

