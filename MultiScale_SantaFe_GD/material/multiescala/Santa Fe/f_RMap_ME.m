function [m_CTHomog,m_DsigmaHomog,hvar_newMacro] = f_RMap_ME(...
   m_IDefMacroReg,hvar_oldMacro,e_DatMatSetMacro,e_VGMacro)

    e_VGMicro     = e_DatMatSetMacro.e_VG;   
    e_VGMicro.INV_PROB = e_VGMacro.INV_PROB;
    e_DatSet      = e_DatMatSetMacro.e_DatSet;     
    e_VGMicro.MOD_TYPE = e_VGMacro.MOD_TYPE;  % Model Flag for microstructure
    
    switch e_VGMicro.MOD_TYPE            
        case 1 % ROM_I
            
            % VARIABLES GLOBALES
            ndoft         = e_VGMicro.ndoft;
            u             = hvar_oldMacro.u;
            e_VarEst_old  = hvar_oldMacro.e_VarEst;
            e_VarAux      = hvar_oldMacro.e_VarAux;
            
            EPS_ROMI      = hvar_oldMacro.EPS_ROMI;
            eta           = hvar_oldMacro.eta;
            Q_ROMI        = [];
            
            % INICIALIZACION DE VARIABLES
            Du_step_old   = zeros(ndoft,1);
            Du_step_new = zeros(ndoft,1);
            
            % Displacement Reduction matrix (Modes)
            %e_VGMicro.ModoPHI = e_VGMacro.ModoPHI;
            
            % Para seleccion de matrix en el solver
            e_VGMicro.esImplex = e_DatSet(1).e_DatMat.esImplex ;            
            
            % DEFORMACION MACRO TOTAL A APLICAR EN CADA PG
            c_DefMacro = arrayfun(@(x)repmat(m_IDefMacroReg,[1,x.e_DatElem.npg,x.nElem]),...
                e_DatSet,'UniformOutput',false);
            
        case 2  % ROM_II
            
            % VARIABLES GLOBALES
            ndoft         =  e_VGMicro.ndoft;
            
            % Hiperreduction Model
            u             = hvar_oldMacro.u;
            e_VarEst_old  = hvar_oldMacro.e_VarEst_ROMII;
            e_VarAux      = hvar_oldMacro.e_VarAux_ROMII;
            
            EPS_ROMI      = hvar_oldMacro.EPS_ROMI;
            eta           = hvar_oldMacro.eta;
            %Q_ROMI      = hvar_oldMacro.Q_ROMI;
            Q_ROMI      = [];
            
            % INICIALIZACION DE VARIABLES
            Du_step_old   = zeros(ndoft,1);
            Du_step_new = zeros(ndoft,1);
            
            e_VGMicro.ndoft  = ndoft;
            e_VGMicro.nElem  = e_DatMatSetMacro.ROM_II.e_DatSet_ROM.nElem ;
            e_VGMicro.dofpe  = e_DatMatSetMacro.ROM_II.e_DatSet_ROM.dofpe ;
            e_VGMicro.ROM_II = e_DatMatSetMacro.ROM_II ;
            
            e_DatSet_ROM = e_DatMatSetMacro.ROM_II.e_DatSet_ROM ;
            
            % Para seleccion de matrix en el solver
            e_VGMicro.esImplex = e_DatSet_ROM.e_DatMat.esImplex ;
            
            % DEFORMACION MACRO TOTAL A APLICAR EN CADA PG
            c_DefMacro = arrayfun(@(x)repmat(m_IDefMacroReg,[1,x.nElem]),...
                e_DatSet_ROM,'UniformOutput',false);
            
        case 0 % HF
            
            % VARIABLES GLOBALES
            ndoft         = e_VGMicro.ndoft;
            u             = hvar_oldMacro.u;
            e_VarEst_old  = hvar_oldMacro.e_VarEst;
            e_VarAux      = hvar_oldMacro.e_VarAux;
            
            % INICIALIZACION DE VARIABLES
            Du_step_old   = zeros(ndoft,1);
            Du_step_new = zeros(ndoft,1);
            
            EPS_ROMI = []; eta=[];  Q_ROMI = [];
            
            % DEFORMACION MACRO TOTAL A APLICAR EN CADA PG
            c_DefMacro = arrayfun(@(x)repmat(m_IDefMacroReg,[1,x.e_DatElem.npg,x.nElem]),...
                e_DatSet,'UniformOutput',false);

        otherwise % HF & TRAINING TEST
            
            % VARIABLES GLOBALES
            ndoft         = e_VGMicro.ndoft;
            u             = hvar_oldMacro.u;
            e_VarEst_old  = hvar_oldMacro.e_VarEst;
            e_VarAux      = hvar_oldMacro.e_VarAux;
            
            % INICIALIZACION DE VARIABLES
            Du_step_old   = zeros(ndoft,1);
            Du_step_new = zeros(ndoft,1);
            
            EPS_ROMI = []; eta=[];  Q_ROMI = [];
            
            % DEFORMACION MACRO TOTAL A APLICAR EN CADA PG
            c_DefMacro = arrayfun(@(x)repmat(m_IDefMacroReg,[1,x.e_DatElem.npg,x.nElem]),...
                e_DatSet,'UniformOutput',false);            
            
    end
   
   % VARIABLES GLOBALES
    nElem        = e_VGMicro.nElem;
    ntens        = e_VGMicro.ntens;
    nSet         = e_VGMicro.nSet;
   
   
%% %Se recupera variables micro JLM
%    xx = e_DatMatSetMacro.xx;
%    omegaMicro = e_DatMatSetMacro.omegaMicro;
%    e_DatSet = e_DatMatSetMacro.e_DatSet;
%    e_VG = e_DatMatSetMacro.e_VG;
%    esImplexMacro = e_DatMatSetMacro.esImplex;
%    e_VarEst_old = hvar_oldMacro.e_VarEst;
%    u = hvar_oldMacro.u;
%    c_GdlCond = hvar_oldMacro.c_GdlCond;
%    Fint         = hvar_oldMacro.Fint;
%    m_LinCond    = hvar_oldMacro.m_LinCond;
%    doff         = hvar_oldMacro.doff;
%    dofl         = hvar_oldMacro.dofl;
%    e_VarAux = hvar_oldMacro.e_VarAux;
%%

% Se recupera variables micro
   xx = e_DatMatSetMacro.xx;
   omegaMicro = e_DatMatSetMacro.omegaMicro;
   
   c_GdlCond = hvar_oldMacro.c_GdlCond;
   
   Fint = hvar_oldMacro.Fint;
   m_LinCond = hvar_oldMacro.m_LinCond;
   doff = hvar_oldMacro.doff;
   dofl = hvar_oldMacro.dofl;

   % INICIALIZACION DE VARIABLES
%    Du_step_old = zeros(ndoft,1); % JLM  % Vector de incrementos de desplazamientos en el paso de tiempo previo
   Fext = zeros(ndoft,1);
   
   %Por si es necesario el paso tiempo a nivel micro.
   e_VG.istep = e_VGMacro.istep;
   %Se guarda que se quiere que el modelo constitutivo se comporte elásticamente.
   e_VG.elast = e_VGMacro.elast;
   %Para impresión y debug se guarda la iteración macro.
   e_VG.iterMacro = e_VGMacro.iter;
   %Para imprensión se guarda el número de elemento macro y número de PG macro.
   e_VG.iElemNumMacro = e_VGMacro.iElemNum;
   e_VG.iPGMacro = e_VGMacro.iPG;
   %
%    e_VG.SharedFS = e_VGMacro.SharedFS;
   
   %Nombre interno de la celda unitaria 
   %(hacer esto en cada iteración puede ser medio lento, ver que hacer sino).
   %e_VG.fileCompleto = [e_VG.fileCompleto,'_EM',int2str(e_VGMacro.iElemNum),'PGM',int2str(e_VGMacro.iPG)];
   %Se guarda los datos del matlabpool (se utiliza para imprimir de pasos y iteraciones a nivel micro que no
   %convergieron).
   e_VG.nLab = e_VGMacro.nLab;
   e_VG.tipoMPool = e_VGMacro.tipoMPool;

   % DESPLAZAMIENTO IMPUESTO
   %No sería necesario estas operaciones porque vfix en todos los modelos clásicos de las
   %formulaciones multiescala son nulos. Además que habría que interpretar como realizar el delta
   %psi_value a nivel micro, ya debería se corresponder con el incremento de tiempo a nivel macro,
   %pero lo que implicaría que en cada paso de tiempo se esté resolviendo un RVE distinto.
   %vDeltaFixTemp = vfix;
   %vDeltaFixTemp(doffCondCte) = vfix(doffCondCte)*(psi_value - psi_value_old);
   %u(doff) = u(doff) + m_InvCRR*vDeltaFixTemp(doff);
  
   % Deformación macro aplicada en la Celda unitaria.
   %Se aplica en forma uniforme en todo dominio.
   %Por simplificidad se considera una deformación macro aplicada distinta por elemento, y no por PG.
   %(no es necesario considerar una estructura para esta variable).
   %m_DefMacro = repmat(eps_new,1,nElem);
   %Deformación macro por elemento y por punto de gauss, dividida en sets.
   %c_DefMacro = arrayfun(@(x)repmat(eps_new,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);
   %Deformación macro por elemento, dividad en sets.
%    c_DefMacro = arrayfun(@(x)repmat(m_IDefMacroReg,[1,x.nElem]),e_DatSet,'UniformOutput',false);%    %JLM
   
   %ticIDNewt = tic;
   % ESQUEMA DE NEWTON-RAPHSON
%    Du_step_new = zeros(ndoft,1);
%    [u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT] = newton_raphson(...
%       xx,m_LinCond,dofl,doff,u,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,e_VarEst_old,e_VarAux,e_DatSet,...
%       c_DefMacro,e_VG);
  
   [u,EPS_ROMI,eta,Q_ROMI,c_GdlCond,Fint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT] = newton_raphson_MICRO(...
       xx,m_LinCond,dofl,doff,u,EPS_ROMI,eta,Q_ROMI,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,e_VarEst_old,e_VarAux,e_DatSet,...
       c_DefMacro,e_VGMicro);
   
%    %Se guarda las fuerzas internas sólo para impresión de las gráficas X-Y
%    %Se está usando equilibrio en totales, por lo que no es necesario la siguiente línea.
%    %Fint = Fint+DFint;
%    %
%    %fprintf('Tiempo del Newton micro: %f\n',toc(ticIDNewt));   
%    
%    % TENSOR TANGENTE HOMOGENEIZADO
%    m_CTHomog = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,...
%         true(nElemTot,1),true(nElemTot,1),e_VG);
%    
%    %Se asume que no se realiza análisis de bifurcación con el tensor tangente constitutivo homogeneizado, por
%    %lo que en el caso ser implex, se devuelve nulo el tensor implícito homogeneizado.
%    if esImplexMacro
%       m_CTHomog = struct('Implex',m_CTHomog,'Impli',zeros(ntens,ntens));
%    end
%    
%    % CÁLCULO DE VARIABLES HOMOGENEIZADAS
%    m_DsigmaHomog = f_HomogArea({e_VarEst_new.sigma},ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VG);
%    %defHomog = f_HomogArea({e_VarEst_new.eps},ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VG);
%    defHomogFl = f_HomogArea({e_VarEst_new.eps_fluct},ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VG);
%    %Verificación de la media del desplazamiento
%    m_uMedioFl = f_MediaDespCU(u,omegaMicro,e_DatSet,e_VG);
%    %
%    fprintf('Elemento %d: PG %d: Norma de la deformación fluctuante media: %g\n',e_VGMacro.iElemNum,...
%       e_VGMacro.iPG,norm(defHomogFl))
%    fprintf('Elemento %d: PG %d: Norma del desplazamiento fluctuante medio: %g\n',e_VGMacro.iElemNum,...
%       e_VGMacro.iPG,norm(m_uMedioFl))
%    
%    hvar_newMacro = struct('u',u,'c_GdlCond',{c_GdlCond},'Fint',Fint,'e_VarEst',e_VarEst_new,...
%       'e_VarAux',e_VarAux,'m_LinCond',m_LinCond,'doff',doff,'dofl',dofl,'c_DefMacro',{c_DefMacro});
%    
% end
 switch e_VGMacro.MOD_TYPE
       
       case 2
           %m_DsigmaHomog = e_DatMatSetMacro.ROM_II.Q_oper*(e_VarEst_new.sigma(:)) ;
           nCOMP_TENS = e_DatMatSetMacro.ROM_II.nCOMP_TENS;
           %W_STRESS   = e_DatMatSetMacro.ROM_II.W_STRESS;
           sigma_GEN = e_VarEst_new.sigma(:) ;
           
           m_DsigmaHomog = e_DatMatSetMacro.ROM_II.Q_oper{1}*(sigma_GEN(nCOMP_TENS{1})) + ...
               e_DatMatSetMacro.ROM_II.Q_oper{2}*(sigma_GEN(nCOMP_TENS{2})) ;
           
       otherwise
           %Para considerar que en el caso del cuadrangulo Q1 se debe integrar las tensiones estabilizadas.
           c_Tens_old = cell(nSet,1);
           c_Tens_new = cell(nSet,1);
           for iSet = 1:nSet
               e_DatElemSet = e_DatSet(iSet).e_DatElem;
               eltype = e_DatElemSet.eltype;
               switch eltype
                   case {2,4,8,31,32,100,108}
                       c_Tens_old{iSet} = e_VarEst_old(iSet).sigma;
                       c_Tens_new{iSet} = e_VarEst_new(iSet).sigma;
                   case 20
                       c_Tens_old{iSet} = e_VarEst_old(iSet).VarHistElem;
                       c_Tens_new{iSet} = e_VarEst_new(iSet).VarHistElem;
                   otherwise
                       error(['Modelo Multiescala: Homogeneizaciï¿½n: Tensiones de homogeneizaciï¿½n: Modelo ',...
                           'constitutivo no definido.'])
               end
           end
           
           m_DsigmaHomog = f_HomogArea(c_Tens_new,ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VGMicro);
           
   end
   
   % TENSOR TANGENTE HOMOGENEIZADO
   %m_CTHomog = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,...
   %     true(ndoft,1),true(ndoft,1),e_VG);
   
   %Se asume que no se realiza anï¿½lisis de bifurcaciï¿½n con el tensor tangente constitutivo homogeneizado, por
   %lo que en el caso ser implex, se devuelve nulo el tensor implï¿½cito homogeneizado.
   %if esImplexMacro
   %   m_CTHomog = struct('Implex',m_CTHomog,'Impli',zeros(ntens,ntens));
   %end
   
   
   % Tensor Tangente Homogeneizado (IMPLEX)
   switch e_VGMacro.MOD_TYPE
       case 1
           % Tensor Tangente Homogeneizado
           m_CTHomog = f_ModTangHomogROM(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,...
               true(nElem,1),true(nElem,1),e_VGMicro);
           
           if e_VGMicro.esImplex
               MAT_LOC = inv(KT.Implex)*KT.Q_ROM ;
           else
               MAT_LOC = inv(KT.Impli)*KT.Q_ROM ;
           end           
           
           m_CTHomog = m_CTHomog-(1/omegaMicro)*KT.Q_ROM'*MAT_LOC;
           
       case 2 %
                    
           nTENS_SET = e_DatMatSetMacro.ROM_II.nTENS_SET;
           BmatRIst_W_red = e_DatMatSetMacro.ROM_II.BmatRIst_W_red;
           
           if e_VGMicro.esImplex
               c_CT_IMPLEX = c_CT{1}(1:2:end);
           else
               c_CT_IMPLEX = c_CT{1};
           end
           
           c_CT_SET_REG = c_CT_IMPLEX(nTENS_SET{1});
           m_TensorTang_mat_IMPLEX_REG = cell2mat(c_CT_SET_REG);
           
           c_CT_SET_DIS = c_CT_IMPLEX(nTENS_SET{2});
           m_TensorTang_mat_IMPLEX_DIS = cell2mat(c_CT_SET_DIS);
           
           Q = [e_DatMatSetMacro.ROM_II.R_oper{1}*m_TensorTang_mat_IMPLEX_REG;
               e_DatMatSetMacro.ROM_II.R_oper{2}*m_TensorTang_mat_IMPLEX_DIS;
               zeros(ntens-1,ntens)]; % JLM CAMBIE 4 POR ntens
           
           D_Taylor = e_DatMatSetMacro.ROM_II.Q_oper{1}*m_TensorTang_mat_IMPLEX_REG + ...
               e_DatMatSetMacro.ROM_II.Q_oper{2}*m_TensorTang_mat_IMPLEX_DIS;
           
           for iSet_COMP = 1:2 % CAMBIAR INDICE!!!!
               
               BmatRIst_W_red_SET = BmatRIst_W_red(nTENS_SET{iSet_COMP}) ;
               c_CT_SET = c_CT_IMPLEX(nTENS_SET{iSet_COMP});
               
               % Esto se puede quitar usando bases no generalizadas sino particulares para cada dominio
               if iSet_COMP==1;
                   p_COMP = 1:e_VGMicro.nModesEPS_bar;
               else
                   p_COMP=e_VGMicro.nModesEPS_bar+1:e_VGMicro.nModesEPS_bar+e_VGMicro.nModesEPS_bar2 ;
               end
               
               for iSetPROD = 1:length(nTENS_SET{iSet_COMP})
                   CB{iSet_COMP}{iSetPROD}=c_CT_SET{iSetPROD}*BmatRIst_W_red_SET{iSetPROD}(:,p_COMP);
               end
               
           end
           
           %calculo de sensibilidades con respecto a la deformacion macro
           Oper = -KT.Implex\Q ;
           
           D_Fluc = zeros(ntens);
           %eta_var = Oper(e_VGMicro.nModesEPS_TOT+1:end,:) ;
           %eta_var = repmat(e_VGMicro.E_MATRIX*[eta_var(1,:);eta_var(2,:);zeros(1,ntens);eta_var(3,:)],e_VGMicro.nElem,1);
           
           for iSet_COMP = 1:2  % CAMBIAR INDICE!!!!
               
               % Esto se puede quitar usando bases no generalizadas sino particulares para cada dominio
               if iSet_COMP==1;
                   p_COMP = 1:e_VGMicro.nModesEPS_bar;
               else
                   p_COMP=e_VGMicro.nModesEPS_bar+1:e_VGMicro.nModesEPS_bar+e_VGMicro.nModesEPS_bar2 ;
               end
               
               CB_mat = CB{iSet_COMP};
               CB_mat = cell2mat(CB_mat');
               
               %D_Fluc = D_Fluc + e_DatMatSetMacro.ROM_II.Q_oper{iSet_COMP}*(CB_mat*Oper(p_COMP,:)+eta_var(e_VGMicro.ROM_II.nCOMP_TENS{iSet_COMP},:));
               D_Fluc = D_Fluc + e_DatMatSetMacro.ROM_II.Q_oper{iSet_COMP}*(CB_mat*Oper(p_COMP,:));
           end
           
           %D_Fluc = D_Fluc - W_STRESS*eta_var(1:ntens,:);
           
           % TENSOR CONSTITUTIVO HOMOGENIZADO
           m_CTHomog = D_Taylor + D_Fluc;
           
       otherwise
           m_CTHomog = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,...
               true(nElem,1),true(nElem,1),e_VGMicro);
   end
   
   
   % CALCULO DE VARIABLES HOMOGENEIZADAS
   % %sigmaHomog = f_HomogArea({e_VarEst_new.sigma},ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VG);
   %defHomog = f_HomogArea({e_VarEst_new.eps},ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VG);
   
   %hvar_newMacro = struct('u',u,'c_GdlCond',{c_GdlCond},'Fint',Fint,'e_VarEst',e_VarEst_new,...
   %   'e_VarAux',e_VarAux,'m_LinCond',m_LinCond,'doff',doff,'dofl',dofl,'m_DefMacro',{c_DefMacro});
   
   switch e_VGMacro.MOD_TYPE
       case 1
           hvar_newMacro = struct('u',u,'EPS_ROMI',EPS_ROMI,'eta',eta,'c_GdlCond',{c_GdlCond},'Fint',Fint,'e_VarEst',e_VarEst_new,...
               'e_VarAux',e_VarAux,'m_LinCond',m_LinCond,'doff',doff,'dofl',dofl,...
               'c_DefMacro',{c_DefMacro});     %,'m_TensProy',m_TensProy
       case 2
           hvar_newMacro = struct('u',u,'EPS_ROMI',EPS_ROMI,'eta',eta,'c_GdlCond',{c_GdlCond},'Fint',Fint,'e_VarEst',[],...
               'e_VarEst_ROMII',e_VarEst_new,'e_VarAux',[],'e_VarAux_ROMII',e_VarAux,...
               'm_LinCond',m_LinCond,'doff',doff,'dofl',dofl,'c_DefMacro',{c_DefMacro});     % ORIGINAL
%            hvar_newMacro = struct('u',u,'c_GdlCond',{c_GdlCond},'Fint',Fint,'e_VarEst',e_VarEst_new,...
%                'e_VarAux',e_VarAux,'m_LinCond',m_LinCond,'doff',doff,'dofl',dofl,...
%                'c_DefMacro',{c_DefMacro});%              u: [72x1 double]
%      c_GdlCond: {2x4 cell}
%           Fint: [72x1 double]
%       e_VarEst: [2x1 struct]
%       e_VarAux: [2x1 struct]
%      m_LinCond: [24x48 double]
%           doff: [72x1 logical]
%           dofl: [72x1 logical]
%     c_DefMacro: {2x1 cell}
       otherwise
           hvar_newMacro = struct('u',u,'c_GdlCond',{c_GdlCond},'Fint',Fint,'e_VarEst',e_VarEst_new,...
               'e_VarAux',e_VarAux,'m_LinCond',m_LinCond,'doff',doff,'dofl',dofl,...
               'c_DefMacro',{c_DefMacro});     %,'m_TensProy',m_TensProy
   end
  
end