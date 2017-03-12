
% HROM - ONLINE
% *************
c_TensorTang = cell(1); 

% Variables globales
nElem = e_VG.nElem;
dofpe = e_VG.dofpe;
ntens = e_VG.ntens;
DefMacro = DefMacro{1} ;

ROM_II = e_VG.ROM_II ;
ROM_ELEM_SET = ROM_II.ROM_ELEM_SET ;
%POINTS_SIGMA =  ROM_II.npoints_sigma ;
nGTOT = ROM_II.nGTOT ;
nTENS_SET = ROM_II.nTENS_SET;

e_DatSet_ROM = ROM_II.e_DatSet_ROM ;
siavarpg     = e_DatSet_ROM.siavarpg;
sihvarpg     = e_DatSet_ROM.sihvarpg;

% MATRICES PRECOMPUTADAS
BmatRIst_W_red = ROM_II.BmatRIst_W_red ;
%IntPhiTGen = ROM_II.IntPhiTGen ;
IntPhiTGen_REG = ROM_II.IntPhiTGen_REG ;
IntPhiTGen_DIS = ROM_II.IntPhiTGen_DIS ;

% Convexification matrices
%K_CONVEX_REG = ROM_II.K_CONVEX_REG ;
%K_CONVEX_DIS = ROM_II.K_CONVEX_DIS ;

% Recuperacion de variables
% (evita llamar desde las estructura de los sets en el bucle de los elementos, ver que solo recupera punteros)
sigma_old        = e_VarEst_old.sigma;
eps_old          = e_VarEst_old.eps;
hvar_old         = e_VarEst_old.hvar;
m_VarHistElemOld = e_VarEst_old.VarHistElem;
aux_var          = e_VarAux.VarAuxGP;

% VARIABLES NUEVA REDUCCION (EN DEFORMACIONES)
nModesEPS_bar = e_VG.nModesEPS_bar ;
nModesEPS_bar2 = e_VG.nModesEPS_bar2 ;
nModesEPS_TOT = e_VG.nModesEPS_TOT; 

%%% Pointers for KT and Fint
p_bar  = 1:nModesEPS_bar ;
p_bar2 = nModesEPS_bar+1:nModesEPS_bar+nModesEPS_bar2 ;

% INICIALIZACIONES
sigma_new        = e_VarEst_new.sigma;
eps_new          = e_VarEst_new.eps;
hvar_new         = e_VarEst_new.hvar;
eps_fluct        = e_VarEst_new.eps_fluct;
m_VarHistElemNew = e_VarEst_new.VarHistElem;

m_VarAuxElem = e_VarAux.VarAuxElem;
conshyp = e_DatSet_ROM.e_DatMat.conshyp;

switch conshyp
%     case 11
%         p_kSD  = ROM_II.e_DatSet_ROM.e_DatElem.pointersVAE.p_kSD;
%         leq    = m_VarAuxElem(p_kSD,:);
    case {100,110}
%         %Matriz 2D de F0
%         uElem = mean(kin_Var,2) ;
%         m_F0 = DefMacro + uElem ;
%         m_F2D = [m_F0(1),m_F0(4);m_F0(5),m_F0(2)];
%         Fzz = m_F0(3);
        tipoFBar = 1;
%         switch tipoFBar
%            case 1  %Plane deformation
%               %Se utiliza lo propuesto por Souza, de no modificar la componente Fzz del tensor.
%               J0 = det(m_F2D);
%               pot = 2;
%               m_FBar = zeros(ntens,1);
%            case 2   %Plane stress ó en el caso que se considere modificar también la Fzz del tensor.
%               J0 = Fzz*det(m_F2D);
%               pot = 3;
%            otherwise
%               error('Matrices Elementales FBar_q1: Tipo de cálculo del FBar no implementado.')
%         end        
end

sigma_new_ROM        = cell(nElem,1); % nElem = nElem_ROM_II

Cconst_PHI_atgauss_ROM      = cell(nElem,1); % nElem = nElem_ROM_II
Cconst_PHI_atgauss_ROM_IMPL = cell(nElem,1); % nElem = nElem_ROM_II

if e_DatSet_ROM.e_DatMat.esImplex
   %m_TensorTang = zeros(ntens,ntens,2,nElem);
   m_TensorTang = cell(2*nElem,1);
else
   %m_TensorTang = zeros(ntens,ntens,nElem);
   m_TensorTang = cell(nElem,1);
end

%eta_voigt = e_VG.E_MATRIX*[eta(1) eta(2) 0 eta(3)]';
%eta_voigt = e_VG.E_MATRIX*eta;


for iROM_GP=1:nElem

    e_DatMatSet  = e_DatSet(ROM_ELEM_SET(iROM_GP)).e_DatMat;
    e_DatElemSet = e_DatSet(ROM_ELEM_SET(iROM_GP)).e_DatElem;
    
    ModoPHI_EPS = BmatRIst_W_red{iROM_GP};
    
%     if ROM_ELEM_SET(iROM_GP)==1; p_COMP=p_bar2; else p_COMP=p_bar; end
%     if ROM_ELEM_SET(iROM_GP)==1; p_COMP=p_bar; else p_COMP=p_bar2; end %JLM
    p_COMP=p_bar;
    if ROM_ELEM_SET(iROM_GP)==3; p_COMP=p_bar2; end  % SOLO PARA RVE_Struct_Hole_01 JLM
    
    kin_Var = ModoPHI_EPS(:,p_COMP)*EPS_ROMI(p_COMP);
    
    % Deformacion fluctuante
    eps_fluct(:,iROM_GP) = kin_Var;

    % Deformacion aplicada a cada punto de Gauss
    eps_new(:,iROM_GP) = DefMacro(:,iROM_GP) + eps_fluct(:,iROM_GP);


   
   
   switch conshyp    
       case {100,110}
           % Gradiente de deformación aplicado a cada punto de Gauss.   
           m_F = eps_new(:,iROM_GP);          
           %Matriz 2D de F
           m_F2D = [m_F(1),m_F(4);m_F(5),m_F(2)];
           Fzz = m_F(3);
          switch tipoFBar
              case 1  %Plane deformation
                 %Se utiliza lo propuesto por Souza, de no modificar la componente Fzz del tensor.         
                 % Determinante del gradiente de deformación.
                 J = det(m_F2D);
                 %Gradiente de deformación FBar
%                  facFBar = (J0/J)^(1/pot);
                 facFBar = 1 ; % esto debe calcularse con J0, lo dejo asi a modo de prueba JLM
                 m_FBar([1;2;4;5]) = facFBar*m_F([1;2;4;5]);
                 %Por si acaso se mantiene el valor que componente zz, aunque debería ser 1.
                 m_FBar(3) = Fzz;
              case 2   %Plane stress ó en el caso que se considere modificar también la Fzz del tensor.
                 % Determinante del gradiente de deformación.
                 J = Fzz*det(m_F2D);            
                 % Gradiente de deformación FBar
%                  facFBar = (J0/J)^(1/pot);
                 facFBar = 1 ; % esto debe calcularse con J0, lo dejo asi a modo de prueba JLM
                 m_FBar = facFBar*m_F;
           end
   end
   
    switch conshyp
        case 2  % J2 Plasticity
            [ct,sigma_new(:,iROM_GP),eps_new(:,iROM_GP),hvar_new(:,iROM_GP),aux_var(:,iROM_GP)] = rmap_plasJ2(...
                eps_new(:,iROM_GP),hvar_old(:,iROM_GP),e_DatMatSet,e_VG);
        case 11 % Softening Regularized damage
            Reg_PAR.ksd = leq(iROM_GP);
            Reg_PAR.leq = leq(iROM_GP);
            [ct,sigma_new(:,iROM_GP),hvar_new(:,iROM_GP),aux_var(:,iROM_GP)] = rmap_damage_BCN(...
                eps_new(:,iROM_GP),hvar_old(:,iROM_GP),e_DatMatSet,e_DatElemSet,e_VG,Reg_PAR);
      case 100 %Elastic Material neo-Hookean JLM
         [ct,sigma_new(:,iROM_GP)] = f_RMap_NeoHook_DefPl(m_FBar,e_DatMatSet,e_VG);
         f_Const = @(m_F)f_RMap_NeoHook_DefPl(m_F,e_DatMatSet,e_VG);
      case 110 %Large deformations J2 Plasticity JLM
         [ct,sigma_new(:,iROM_GP),hvar_new(:,iROM_GP),aux_var(:,iROM_GP),strsg,dmatx] = ...
            f_RMapPlastJ2LD(m_FBar,hvar_old(:,iROM_GP),aux_var(:,iROM_GP),e_DatMatSet,e_VG);
         %Tensor tangente elástico de Henky en función de G y K.
         %G=e_DatMatSet.young/2/(1+e_DatMatSet.poiss);K=e_DatMatSet.young/3/(1-2*e_DatMatSet.poiss);d1=4/3*G+K;d2=K-2/3*G;d3=G;c=[d1,d2,0,d2;d2,d1,0,d2;0,0,d3,0;d2,d2,0,d1];
         %f_Const = @(m_FBar,m_F)f_RMapPlastJ2LD(m_FBar,hvar_old(:,iPG),aux_var(:,iPG),e_DatMatSet,e_VG);
      otherwise
         error('Matrices Elementales FBar_q1: Modelo constitutivo no definido.')
         
    end
    
    % storage of stresses for assessment of Fint
    %sigma_new_ROM{iROM_GP} = sigma_new(:,iROM_GP)+eta_voigt;
    sigma_new_ROM{iROM_GP} = sigma_new(:,iROM_GP);
    
    %sigma_new(:,iROM_GP) = sigma_new(:,iROM_GP)+eta_voigt; 
    sigma_new(:,iROM_GP) = sigma_new(:,iROM_GP);
    
    if e_DatSet_ROM.e_DatMat.esImplex
        % IMPLEX: En el caso de tener activado el implex se almacena las matrices constitutivas
        %tangentes implicitas para el anï¿½lisis de bifurcacion.
        %Por ahora se asume que el modelo constitutivo siempre devuelve una estructura cuando esta seleccionado
        %el implex.
        
        % m_TensorTang(:,:,1,iROM_GP) = ct.Implex;
        m_TensorTang{(iROM_GP-1)*2+1} = ct.Implex;
        
        %Se almacena los tensores tangentes constitutivo implicitos para anï¿½lisis de bifurcaciï¿½n como si fuera
        %PG adicionales, tantos como nPG. Se almacena en los indices (:,:,nPG+1:2*nPG).
        %Hay que tener cuidado cuando se hace un loop por la componente de PG sea hasta nPG, y no usando la
        %tercera dimension de m_TensorTang (size(m_TensorTang,3)).
        
        % m_TensorTang(:,:,2,iROM_GP) = ct.Impli;
        m_TensorTang{iROM_GP*2} = ct.Impli;
        
        Cconst_PHI_atgauss_ROM{iROM_GP} =  ct.Implex*ModoPHI_EPS(:,p_COMP);
        Cconst_PHI_atgauss_ROM_IMPL{iROM_GP} = ct.Impli*ModoPHI_EPS(:,p_COMP);
        
        %Se almacena para la homogenizaciï¿½n del tensor tangente
        % m_TensorTang{iROM_GP} = c_tilde;
    else
        % m_TensorTang(:,:,1,iROM_GP) = ct;
        m_TensorTang{iROM_GP} = ct;
        
        %Cconst_PHI_atgauss_ROM{iROM_GP} =  ct.Impli*ModoPHI_EPS(:,p_COMP);
        Cconst_PHI_atgauss_ROM{iROM_GP} =  ct*ModoPHI_EPS(:,p_COMP);
        Cconst_PHI_atgauss_ROM_IMPL{iROM_GP} =  Cconst_PHI_atgauss_ROM{iROM_GP};
        
        %Se almacena para la homogenizaciï¿½n del tensor tangente
        % m_TensorTang{iROM_GP} = cconst;
    end
    
end

B_II_sigma_W  = ROM_II.B_II_sigma_W ;
%B_II_diff_c_q = ROM_II.B_II_diff_c_q ;

%eta_matrix = repmat(eta_voigt,nElem,1);

% Derivada del multiplicador de Lagrange de forma matricial
A=zeros(ntens,ntens-1); A(1,1)=1; A(2,2)=1; A(end,end)=1;
A=e_VG.E_MATRIX*A;

% Calculo de fuerzas internas reducidas
%sigma_new_ROM = cell2mat(sigma_new_ROM);
%Fint = [ROM_II.B_II_sigma_W*sigma_new_ROM+eta_matrix;
%         IntPhiTGen'*EPS_ROMI];
sigma_new_ROM_1 = cell2mat(sigma_new_ROM(nTENS_SET{1}));
sigma_new_ROM_2 = cell2mat(sigma_new_ROM(nTENS_SET{2}));
%Fint = [B_II_sigma_W{1}*(sigma_new_ROM_1+eta_matrix(ROM_II.nCOMP_TENS{1}));
%        B_II_sigma_W{2}*(sigma_new_ROM_2+eta_matrix(ROM_II.nCOMP_TENS{2}));
%         IntPhiTGen_REG'*EPS_ROMI(p_bar)+ IntPhiTGen_DIS'*EPS_ROMI(p_bar2)];     

%Fint = [B_II_sigma_W{1}*sigma_new_ROM_1 + B_II_diff_c_q{1}*(EPS_ROMI(p_bar)-Q_ROMI(p_bar));
%        B_II_sigma_W{2}*sigma_new_ROM_2 + B_II_diff_c_q{2}*(EPS_ROMI(p_bar2)-Q_ROMI(p_bar2));
%        IntPhiTGen_REG'*EPS_ROMI(p_bar)+ IntPhiTGen_DIS'*EPS_ROMI(p_bar2);
%        -B_II_diff_c_q{1}*(EPS_ROMI(p_bar)-Q_ROMI(p_bar));
%        -B_II_diff_c_q{2}*(EPS_ROMI(p_bar2)-Q_ROMI(p_bar2))];
    
Fint = [B_II_sigma_W{1}*sigma_new_ROM_1 + IntPhiTGen_REG*eta ;
        B_II_sigma_W{2}*sigma_new_ROM_2 + IntPhiTGen_DIS*eta ;
        IntPhiTGen_REG'*EPS_ROMI(p_bar)+ IntPhiTGen_DIS'*EPS_ROMI(p_bar2)];
     
% Calculo de la matriz de rigidez reducida
%KT.Implex =  ROM_II.B_II_sigma_W*Cconst_B_atgauss_ROM ;
%Cconst_PHI_atgauss_ROM = cell2mat(Cconst_PHI_atgauss_ROM);
%KT.Implex = [ ROM_II.B_II_sigma_W*Cconst_PHI_atgauss_ROM  IntPhiTGen;
%               IntPhiTGen'   zeros(ntens-1) ];
Cconst_PHI_atgauss_ROM_1 = cell2mat(Cconst_PHI_atgauss_ROM(nTENS_SET{1}));           
Cconst_PHI_atgauss_ROM_2 = cell2mat(Cconst_PHI_atgauss_ROM(nTENS_SET{2}));  
%KT.Implex = [B_II_sigma_W{1}*Cconst_PHI_atgauss_ROM_1    zeros(nModesEPS_bar,nModesEPS_bar2)     IntPhiTGen_REG;
%               zeros(nModesEPS_bar2,nModesEPS_bar)     B_II_sigma_W{2}*Cconst_PHI_atgauss_ROM_2  IntPhiTGen_DIS;
%                         IntPhiTGen_REG'                           IntPhiTGen_DIS'               zeros(ntens-1) ];

%KT.Implex = [B_II_sigma_W{1}*Cconst_PHI_atgauss_ROM_1    zeros(nModesEPS_bar,nModesEPS_bar2)     B_II_sigma_W{1}*repmat(A,length(nTENS_SET{1}),1);
%               zeros(nModesEPS_bar2,nModesEPS_bar)     B_II_sigma_W{2}*Cconst_PHI_atgauss_ROM_2  B_II_sigma_W{2}*repmat(A,length(nTENS_SET{2}),1);
%                         IntPhiTGen_REG'                           IntPhiTGen_DIS'               zeros(ntens-1) ];

%KT.Implex = [B_II_sigma_W{1}*Cconst_PHI_atgauss_ROM_1           zeros(nModesEPS_bar,nModesEPS_bar2)     B_II_sigma_W{1}*repmat(A,length(nTENS_SET{1}),1);
%               zeros(nModesEPS_bar2,nModesEPS_bar)           B_II_sigma_W{2}*Cconst_PHI_atgauss_ROM_2   B_II_sigma_W{2}*repmat(A,length(nTENS_SET{2}),1);
%         (B_II_sigma_W{1}*repmat(A,length(nTENS_SET{1}),1))'     (B_II_sigma_W{2}*repmat(A,length(nTENS_SET{2}),1))'      zeros(ntens-1) ];

KT.Implex = [B_II_sigma_W{1}*Cconst_PHI_atgauss_ROM_1           zeros(nModesEPS_bar,nModesEPS_bar2)     IntPhiTGen_REG;
               zeros(nModesEPS_bar2,nModesEPS_bar)           B_II_sigma_W{2}*Cconst_PHI_atgauss_ROM_2   IntPhiTGen_DIS;
                         IntPhiTGen_REG'                                  IntPhiTGen_DIS'               zeros(ntens-1) ];


%KT.Implex = [B_II_sigma_W{1}*Cconst_PHI_atgauss_ROM_1    zeros(nModesEPS_bar,nModesEPS_bar2)     IntPhiTGen_REG;
%               zeros(nModesEPS_bar2,nModesEPS_bar)     B_II_sigma_W{2}*Cconst_PHI_atgauss_ROM_2  IntPhiTGen_DIS;
%                         IntPhiTGen_REG'                           IntPhiTGen_DIS'               zeros(ntens-1) ];


% CONVEXIFICATION METHOD W/ LAGRANGE MULTIPLIER RECONSTRUCTION
%KT.Implex = [B_II_sigma_W{1}*Cconst_PHI_atgauss_ROM_1+B_II_diff_c_q{1}    zeros(nModesEPS_bar,nModesEPS_bar2)           IntPhiTGen_REG                           -B_II_diff_c_q{1}                   zeros(nModesEPS_bar,nModesEPS_bar2);
%               zeros(nModesEPS_bar2,nModesEPS_bar)        B_II_sigma_W{2}*Cconst_PHI_atgauss_ROM_2+B_II_diff_c_q{2}     IntPhiTGen_DIS                 zeros(nModesEPS_bar2,nModesEPS_bar)                 -B_II_diff_c_q{2};
%                         IntPhiTGen_REG'                                        IntPhiTGen_DIS'                         zeros(ntens-1)                     zeros(ntens-1,nModesEPS_bar)                zeros(ntens-1,nModesEPS_bar2) ;
%                         -B_II_diff_c_q{1}                                zeros(nModesEPS_bar,nModesEPS_bar2)       zeros(nModesEPS_bar,ntens-1)                  B_II_diff_c_q{1}                  zeros(nModesEPS_bar,nModesEPS_bar2);         
%               zeros(nModesEPS_bar2,nModesEPS_bar)                         -B_II_diff_c_q{2}                        zeros(nModesEPS_bar2,ntens-1)      zeros(nModesEPS_bar2,nModesEPS_bar)                  B_II_diff_c_q{2}             ];

                     
%KT.Impli = ROM_II.B_II_sigma_W*Cconst_B_atgauss_ROM_IMPL ;
%Cconst_PHI_atgauss_ROM_IMPL = cell2mat(Cconst_PHI_atgauss_ROM_IMPL);
%KT.Impli = [ ROM_II.B_II_sigma_W*Cconst_PHI_atgauss_ROM_IMPL IntPhiTGen;
%               IntPhiTGen'   zeros(ntens-1) ];
Cconst_PHI_atgauss_ROM_IMPL_1 = cell2mat(Cconst_PHI_atgauss_ROM_IMPL(nTENS_SET{1}));
Cconst_PHI_atgauss_ROM_IMPL_2 = cell2mat(Cconst_PHI_atgauss_ROM_IMPL(nTENS_SET{2}));
%KT.Impli = [B_II_sigma_W{1}*Cconst_PHI_atgauss_ROM_IMPL_1        zeros(nModesEPS_bar,nModesEPS_bar2)     IntPhiTGen_REG;
%                  zeros(nModesEPS_bar2,nModesEPS_bar)     B_II_sigma_W{2}*Cconst_PHI_atgauss_ROM_IMPL_2  IntPhiTGen_DIS;
%                             IntPhiTGen_REG'                              IntPhiTGen_DIS'                zeros(ntens-1) ];

%KT.Impli = [B_II_sigma_W{1}*Cconst_PHI_atgauss_ROM_IMPL_1        zeros(nModesEPS_bar,nModesEPS_bar2)     B_II_sigma_W{1}*repmat(A,length(nTENS_SET{1}),1);
%                  zeros(nModesEPS_bar2,nModesEPS_bar)     B_II_sigma_W{2}*Cconst_PHI_atgauss_ROM_IMPL_2  B_II_sigma_W{2}*repmat(A,length(nTENS_SET{2}),1);
%         (B_II_sigma_W{1}*repmat(A,length(nTENS_SET{1}),1))'    (B_II_sigma_W{2}*repmat(A,length(nTENS_SET{2}),1))'           zeros(ntens-1) ];

KT.Impli = [B_II_sigma_W{1}*Cconst_PHI_atgauss_ROM_IMPL_1        zeros(nModesEPS_bar,nModesEPS_bar2)     IntPhiTGen_REG;
                  zeros(nModesEPS_bar2,nModesEPS_bar)     B_II_sigma_W{2}*Cconst_PHI_atgauss_ROM_IMPL_2  IntPhiTGen_DIS;
                          IntPhiTGen_REG'                                  IntPhiTGen_DIS'               zeros(ntens-1) ];


%KT.Impli = [B_II_sigma_W{1}*Cconst_PHI_atgauss_ROM_IMPL_1        zeros(nModesEPS_bar,nModesEPS_bar2)     IntPhiTGen_REG;
%                  zeros(nModesEPS_bar2,nModesEPS_bar)     B_II_sigma_W{2}*Cconst_PHI_atgauss_ROM_IMPL_2  IntPhiTGen_DIS;
%                             IntPhiTGen_REG'                              IntPhiTGen_DIS'                zeros(ntens-1) ];


% CONVEXIFICATION METHOD W/ LAGRANGE MULTIPLIER RECONSTRUCTION
%KT.Impli = [B_II_sigma_W{1}*Cconst_PHI_atgauss_ROM_IMPL_1+B_II_diff_c_q{1}                 zeros(nModesEPS_bar,nModesEPS_bar2)             IntPhiTGen_REG                           -B_II_diff_c_q{1}                   zeros(nModesEPS_bar,nModesEPS_bar2);
%               zeros(nModesEPS_bar2,nModesEPS_bar)                      B_II_sigma_W{2}*Cconst_PHI_atgauss_ROM_IMPL_2+B_II_diff_c_q{2}     IntPhiTGen_DIS                 zeros(nModesEPS_bar2,nModesEPS_bar)                 -B_II_diff_c_q{2};
%                         IntPhiTGen_REG'                                                           IntPhiTGen_DIS'                     zeros(ntens-1)                    zeros(ntens-1,nModesEPS_bar)             zeros(ntens-1,nModesEPS_bar2) ;
%                         -B_II_diff_c_q{1}                                                  zeros(nModesEPS_bar,nModesEPS_bar2)       zeros(nModesEPS_bar,ntens-1)                  B_II_diff_c_q{1}                   zeros(nModesEPS_bar,nModesEPS_bar2);         
%               zeros(nModesEPS_bar2,nModesEPS_bar)                                                 -B_II_diff_c_q{2}                  zeros(nModesEPS_bar2,ntens-1)       zeros(nModesEPS_bar2,nModesEPS_bar)                  B_II_diff_c_q{2}             ];


% Almacenamiento de los tensores constitutivos
c_TensorTang{1} = m_TensorTang;

% Actualizacion de las variables internas y de estado
e_VarEst_new.sigma = sigma_new; 
e_VarEst_new.eps   = eps_new;
e_VarEst_new.hvar  = hvar_new;
e_VarEst_new.eps_fluct = eps_fluct;
e_VarEst_new.VarHistElem = m_VarHistElemNew;
e_VarAux.VarAuxGP  = aux_var;
