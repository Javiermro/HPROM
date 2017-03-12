function [PHI_GEN,IntPhiTGen,MAT_WEAK,ValRM] = HROM_OFFLINE_PROC_ROMI(e_DatSet,ModoPHI_REG,ModoPHI_DIS,e_VG,isMICRO)

% **********************************
% Definicion del set multiescala 
% Nota: en este caso inicialmente solo se considerara un posible tipo de microcelda..
% **********************************

%ME_Set=find(conshyp_GEN==53);
nGTOT=0;
if isMICRO %~e_VG.esME % MONOSCALE CASE
    e_VG_ME = e_VG;
    nSet = e_VG.nSet;
    for iSet = 1:nSet
        nGTOT = nGTOT + e_DatSet(iSet).e_DatElem.npg*e_DatSet(iSet).nElem;
    end
else
    e_VG_ME = e_DatSet.e_DatMat.e_VG;
    nSet = e_DatSet.e_DatMat.e_VG.nSet;
    for iSet = 1:nSet
        nGTOT = nGTOT + e_DatSet.e_DatMat.e_DatSet(iSet).e_DatElem.npg*...
            e_DatSet.e_DatMat.e_DatSet(iSet).nElem;
    end
end

ntens = e_VG_ME.ntens ;

BmatRIstT  = zeros(e_VG_ME.nModesEPS_TOT,ntens*nGTOT);
BmatRIst_W = cell(nGTOT,1);             % (% BmatRIst, not affected by the WEIGHTS)
IntPhiTGen  = zeros(e_VG_ME.nModesEPS_TOT,ntens-1);

% General B Matrix (for retriving displacements from strains)
BmatGEN_W = zeros(ntens*nGTOT,e_VG_ME.ndoft);

KG_WEAK = zeros(e_VG_ME.ndoft);

ponder_factors = zeros(nGTOT,1) ;
for iSet = 1:nSet
    
    if isMICRO %~e_VG.esME;  
        e_Dat_iSet = e_DatSet(iSet); 
    else
        e_Dat_iSet = e_DatSet.e_DatMat.e_DatSet(iSet);
    end
    
    nElem = e_Dat_iSet.nElem;
    eltype = e_Dat_iSet.e_DatElem.eltype;
    npg = e_Dat_iSet.e_DatElem.npg;
    
    % B Matrix by set
    %m_BT_SET = e_DatSet.e_DatMat.e_DatSet(iSet).m_BT;
    m_BT_SET = e_Dat_iSet.m_BT;
    
    % Element id
    %m_NumElem =  e_DatSet.e_DatMat.e_DatSet(iSet).m_NumElem ;
    m_NumElem = e_Dat_iSet.m_NumElem ;
    
    %for integration of the weak problem (recover displacements from strains)
    %m_DetJT = e_DatSet.e_DatMat.e_DatSet(iSet).m_DetJT;
    m_DetJT = e_Dat_iSet.m_DetJT;
    
    %wg = e_DatSet.e_DatMat.e_DatSet(iSet).e_DatElem.wg;
    wg = e_Dat_iSet.e_DatElem.wg;
    
    m_pesoPG = bsxfun(@times,m_DetJT,wg);    
    
    % Loop over the elements in each Set
    for iElem=1:nElem
        dof = f_DofElem(e_Dat_iSet.conec(iElem,:),e_VG_ME.ndn);
        ponder_factors((m_NumElem(iElem)-1)*npg+1:m_NumElem(iElem)*npg) = e_Dat_iSet.m_VolElemGP(:,iElem);
        
        % Set sobre el tipo de elementos en la microescala
        switch eltype
            case 2
            case {4,8,31,32,108}
                e_PGsID = (m_NumElem(iElem)-1)*ntens*npg+1:m_NumElem(iElem)*ntens*npg ;                
                %PHI_GEN = [ModoPHI_REG(e_PGsID,:) ModoPHI_DIS(e_PGsID,:) ModoPHI_ELAS2(e_PGsID,:)] ;
                if e_VG.nModesEPS_bar2 == 0 %No Domain decomposition
                    PHI_GEN = [ModoPHI_REG(e_PGsID,:)] ;
                else
                    PHI_GEN = [ModoPHI_REG(e_PGsID,:) ModoPHI_DIS(e_PGsID,:)] ;
                end                
                e_BT = m_BT_SET(:,:,:,iElem);                
                for ipg = 1:npg
                    ilocp = (m_NumElem(iElem)-1)*npg+ipg ;

                    % Definint the ilocp entry of BmatRIstT
                    iini = (ilocp-1)*ntens +1;
                    ifin = ilocp*ntens ;
                    
                    iPG_ID = (ipg-1)*ntens+1:ipg*ntens ;                    
                    BmatRIT = PHI_GEN(iPG_ID,:)' ;
                    
                    BmatGEN_W(e_PGsID(iPG_ID),dof) = e_BT(:,:,ipg)*m_pesoPG(ipg,iElem);                    
                    
                    %  WHERE_WEIGHTS = 'WinB' ;
                    % WinB: Weighting factors are included in B-matrix (default)
                    % WinS: Weighting factors are included in SNAPSHOTS
                    BmatRIstTW = BmatRIT*e_Dat_iSet.m_VolElemGP(ipg,iElem) ;
                    BmatRIstT(:,iini:ifin) = BmatRIstTW;
                    BmatRIst_W{(m_NumElem(iElem)-1)*npg+ipg} = BmatRIT';
                    
                    if eltype==108
                        IntPhiTGen = IntPhiTGen + BmatRIstTW(:,[1 2 4 5]) ;
                    else
                        %IntPhiTGen = IntPhiTGen + BmatRIstT(:,iini:ifin) ;
                        IntPhiTGen = IntPhiTGen + BmatRIstTW(:,[1 2 4]) ;
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

IntPhiTGen = e_VG_ME.E_MATRIX*IntPhiTGen;

if e_VG.nModesEPS_bar2 == 0 %No Domain decomposition
    PHI_GEN = [ModoPHI_REG] ;
else
    PHI_GEN = [ModoPHI_REG ModoPHI_DIS] ;
end  


%MATRICES FOR SOLVING WEAK PROBLEM (FROM STRAINS TO DISPLACEMENTS)
ValRM = true(e_VG_ME.ndoft,1);
ValRM(e_VG_ME.doff_RM,1) = false;

Fint_WEAK = BmatGEN_W'*PHI_GEN;

MAT_WEAK = KG_WEAK(ValRM,ValRM)\Fint_WEAK(ValRM,:);

