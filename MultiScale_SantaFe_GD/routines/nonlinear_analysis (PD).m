function [u,OTHEROUT] = nonlinear_analysis(in,xx,m_SetElem,f,funbc,e_DatSet,e_VG,isMICRO,varargin)

%***************************************,***************************************************
%*  RUTINA PARA RESOLVER EL PROBLEMA NO LINEAL MATERIAL MEDIANTE M.E.F.                   *
%*                                                                                        *                  
%*  ARGUMENTOS DE ENTRADA:                                                                *                  
%*  file    : nombre del archivo de datos para el problema en estudio                     *           
%*  in     : lista de nodos                                                              *
%*  xx      : lista de coordenadas                                                        *
%*  conec   : lista de conectividades                                                     *
%*  vfix    : vector de desplazamientos impuestos                                         *
%*  f       : vector de cargas externas aplicadas                                         *
%*  funbc   : funcion temporal para aplicar condiciones de borde                          *
%*  Eprop   : lista de propiedades de los elementos                                       *
%*                                                                                        *                  
%*  ARGUMENTOS DE SALIDA:                                                                 *                  
%*  u       : vector de desplazamientos                                                   *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

%varagin_INIT=varargin;
%[OTHEROUT,~,INDEX_INV_VAR_PLOT,USNAP,SIGMA_SNAP,ITERATION_COUNTER,BIFURCATION_FLAG,...
%    STATE_STEP,DISSIP_STEP,TIME_ONLINE,INV_VAR_PLOT] = INIT_DATA_TRAIN(e_DatSet,e_VG,varagin_INIT);
TR_VAR = false;

% VARIABLES GLOBALES
ndoft              = e_VG.ndoft;
nElem              = e_VG.nElem;
ndn                = e_VG.ndn;
np                 = e_VG.np;
postpro_impre_step = e_VG.postpro_impre_step;
CONTROL_STRAT      = e_VG.CONTROL_STRAT;
%Dtime              = e_VG.Dtime;
ndime              = e_VG.ndime;
ntens              = e_VG.ntens;
nSet               = e_VG.nSet;
Dtime              = f.Dtime;

% INICIALIZACION DE VARIABLES
Fint            = zeros(ndoft,1);             % Vector de fzas internas
Fext            = zeros(ndoft,1);             % Vector de fzas internas
Fext_inic       = zeros(ndoft,1);
u               = zeros(ndoft,1);             % Vector de desplazamientos totales
u_old           = zeros(ndoft,1);             % vector de desplazamientos totales del paso anterior
Du_step_old     = zeros(ndoft,1);             % Vector de incrementos de desplazamientos en el paso de tiempo previo
dofl_old        = zeros(ndoft,1);
u_old_Stage     = zeros(ndoft,1);
sigmaHomog      = zeros(ntens,1);
epsilon_Macro   = zeros(ntens,1);
%hvar_old        = zeros(sihvare,nElem);       % Vector de variables internas del paso previo
%eps_old         = zeros(sitvare,nElem);       % Vector de deformaciones previas
%sigma_old       = zeros(sitvare,nElem);       % Vector de tensiones previas
%Por simplicidad se considera una deformaci�n macro aplicada distinta por elemento, y no por PG.
%(no es necesario considerar una estructura para esta variable)
%DefMacro = zeros(ntens*npg,nElem);
%DefMacro = zeros(ntens,nElem);

% Only for ROM purposes
EPS_ROMI=[]; eta=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DefMacro = arrayfun(@(x)zeros(e_VG.ntens,x.e_DatElem.npg,x.nElem),e_DatSet,'UniformOutput',false);
if isMICRO.MICRO
    %epsilon_Macro0=[5e-2;0;0;0];
    %epsilon_Macro =[0;0;0;0];
    epsilon_Macro0 = isMICRO.epsilon_Macro0 ; % Increment of macro-strain
    epsilon_Macro  = isMICRO.epsilon_Macro ; % Initial imposed macro-strain
    
    % Only for ROM purposes
    if e_VG.MOD_TYPE~=0
        EPS_ROMI=zeros(e_VG.nModesEPS_TOT,1); eta=zeros(e_VG.ntens-1,1);
    %else
    %    EPS_ROMI=[]; eta=[];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%DefMacroT = zeros(ntens,nElem);

ELOCCalc = false(1,nElem);
%eps_fluct_old   = zeros(sitvare,nElem);       % Vector de deformaciones fluctuantes previas
if (~e_VG.esME && e_VG.MOD_TYPE==2)
    e_VarEst_old = f_eVarEstInic({'sigma','eps','eps_fluct','hvar','delta_eps'},e_VG.ROM_II.e_DatSet_ROM,1);
    e_VarAux = f_eVarAuxInic(e_VG.ROM_II.e_DatSet_ROM,1);
else
    %e_VarEst_new = f_eVarEstInic({'sigma','eps','eps_fluct','hvar','delta_eps'},e_DatSet,e_VG.nSet,0);
    e_VarEst_old = f_eVarEstInic({'sigma','eps','eps_fluct','hvar'},e_DatSet,nSet,e_VG.MOD_TYPE);
    e_VarAux = f_eVarAuxInic(e_DatSet,nSet);
end

%e_VarEst_new = e_VarEst_old;
c_GdlCond = f_cVarCondInic(e_DatSet,e_VG);

%Impresion de archivo de postprocesado de la malla y inicializaci�n del archivo de datos
matlab2gid_mesh(in,xx,e_DatSet,e_VG)
e_VG.istep = 0;
f_InicArchDat(in,m_SetElem,e_DatSet,e_VG)

%Ploteo en el paso 0 (se asume que tiene todos valores nulos en tiempo 0)
%Esto ultimo habra que ver las condiciones de borde y la funcion psi_value (por ejemplo que esta no
%sea nula en el tiempo 0)
%sitvare = e_VG.sitvare;
%nnod    = e_VG.nnod;
%matlab2gid_res(0,u,zeros(sitvare,nElem),zeros(sitvare,nElem),zeros(sitvare,nElem),DefMacro,...
%   zeros(sihvare,nElem),e_VG,zeros(nnod,ntens),zeros(ndime,nnod));
f_SaveDatGraf(u,c_GdlCond,Fint,e_VarEst_old,e_VarAux,e_DatSet,m_SetElem,sigmaHomog,epsilon_Macro,e_VG)

%DEF_SINPG = [];

%INTEGRACION TEMPORAL POR STAGE
for iStage = 1:e_VG.nStage
    
    % Variables correspondientes a cada Stage
    npStage    = np(iStage);
    e_VG.Dtime = Dtime(iStage);
    e_VG.iStage = iStage;
    
    % VALOR DE LA FUNCION Psi EN EL TIEMPO CERO
    psi_value_old      = get_psi_value_Stage(funbc.m_func{iStage},0,e_VG);
    %Fext               = f*psi_value_old;
    
    %Con istepSave se indica desde que paso se empieza a correr. Si lee el archivo mat, el valor nulo de
    %istepSave es pisado por el paso donde se guard� el workspace.
    istepSave = 0;
    [dirFile,nomFile] = fileparts(e_VG.fileCompleto);
    if 0
        load(fullfile(dirFile,'PasoSalvado'))
        fprintf('Se recuper� los valores de las variables del paso %d.\n',istepSave)
    end
    
    % INFORMACION ADICIONAL
    %Se lee un script de matlab. Esto permite cambiar algunas propiedades, matrices, etc., antes de entrar en el
    %calculo en forma rapida sin modificar demasiado el programa.
    %Ocurre un error al usar usar run, no se si pasa lo mismo con eval, ya que matlab no se da cuenta
    %que un script o funci�n fue modificada para precompilarla (usa la precompilaci�n anterior). Esto
    %hace que las modificaciones del script las ignora y usa por ejemplo valores de variable que son de
    %la version del script anterior.
    %Esto se arregla con clear all, pero eso puede traer muchos problemas, adem�s que se desaprovecha
    %las precompilaciones previas. Se borra solo la precompilaci�n de la funci�n (hay que usar el nombre
    %de la funci�n solo, sin camino).
    if exist([e_VG.fileCompleto,'.m'],'file')
        clear(nomFile)
        %run corre el script sin estar en el path o que en el directorio activo
        run([e_VG.fileCompleto,'.m'])
    end
    
    % FUNCION CONDICIONES DE BORDE
    [m_LinCond,vfix,m_InvCRR,doff,dofl,doffCondCte] = f_CondBord(e_VG.m_CondBord.m_SCondBord{iStage},xx,...
        ndoft,ndn,ndime,nSet,e_DatSet,e_VG.m_ConecFront,e_VG.RotacionGlobal);
    
    % Verificacion nuevos grados de libertad libres
    if iStage>1
        diff_dofl=dofl-dofl_old;
        positive=find(diff_dofl>0);
        Fext(positive)=Fint(positive);
        if (CONTROL_STRAT == 4)
            Fext_inic(positive)=Fint(positive);
            Fext_inic(dofl)=Fint(dofl);
        end
    end
    
    % INICIALIZACION DEL VECTOR DE CARGAS ACUMULADAS DEL STAGE
    Fext_Stage = zeros(ndoft,1);  
    
    % INTEGRACION TEMPORAL POR TIME-STEP EN CADA STAGE
    for istep = istepSave+1:npStage
        
        ticStep = tic;
        e_VG.istep = istep;
        time = Dtime(iStage)*istep;
        
        fprintf('STAGE: %-2d - STEP: %-3d\n',iStage,istep);
        
        % VALOR ACTUAL DE LA FUNCION TEMPORAL PARA ESCALAR CONDICIONES DE BORDE
        psi_value = get_psi_value_Stage(funbc.m_func{iStage},time,e_VG);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if isMICRO.MICRO
            epsilon_Macro = epsilon_Macro + epsilon_Macro0*(psi_value - psi_value_old);
            if e_VG.MOD_TYPE==2 % MONOSCALE + HROM
                DefMacro =arrayfun(@(x)repmat(epsilon_Macro,[1,x.e_DatElem.npg,x.nElem]),e_VG.ROM_II.e_DatSet_ROM,'UniformOutput',false);
            else % MONOSCALE + ROM & FE
                DefMacro =arrayfun(@(x)repmat(epsilon_Macro,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);
            end
            DefMacroT = epsilon_Macro;            
        else
            DefMacro = arrayfun(@(x)zeros(e_VG.ntens,x.e_DatElem.npg,x.nElem),e_DatSet,'UniformOutput',false);
            DefMacroT = zeros(ntens,1);            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % DESPLAZAMIENTO IMPUESTO
        vDeltaFixTemp  = vfix;
        Du_step_new    = zeros(ndoft,1);
        
        for iLoad=1:f.NTLoad(iStage)
            % FUERZA EXTERNA IMPUESTA EN CADA CASO DE CARGA
            if (CONTROL_STRAT == 1)
                doffStage                = e_VG.m_CondBord.m_GDL{iStage}{iLoad};
                Fext_Stage               = f.f_Stages{iStage}(:,iLoad)*(psi_value(iLoad) - psi_value_old(iLoad));
                vDeltaFixTemp(doffStage) = vfix(doffStage)*(psi_value(iLoad) - psi_value_old(iLoad));
            else
                if (istep == 1)
                    Fext_Stage           = Fext_Stage + f.f_Stages{iStage}(:,iLoad);
                end
            end
        end
 
        % FUERZA EXTERNA IMPUESTA TOTAL DEL STAGE (VALORES ACUMULADOS)
        if (CONTROL_STRAT == 1)
            Fext              = Fext + Fext_Stage;
            u(doff)           = u(doff) + m_InvCRR*vDeltaFixTemp(doff);
            Du_step_new(doff) = m_InvCRR*vDeltaFixTemp(doff);
        else
            Fext              = Fext_Stage;
            e_VG.vfix         =  vfix;
            e_VG.vfix_doff    =  m_InvCRR*vfix(doff);
            Du_step_new(doff) =  u_old_Stage(doff) + m_InvCRR*vfix(doff)*e_VG.lambda(iStage) - u(doff);
            u(doff)           =  u_old_Stage(doff) + m_InvCRR*vfix(doff)*e_VG.lambda(iStage);
        end
        
        ticIDNewt = tic; ticSPEEDUP = tic;
        % ESQUEMA DE NEWTON-RAPHSON
        [u,EPS_ROMI,eta,c_GdlCond,Fint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT,lambda,iter] = newton_raphson(xx,m_LinCond,...
            dofl,doff,u,EPS_ROMI,eta,u_old_Stage,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,Fext_inic,e_VarEst_old,e_VarAux,e_DatSet,DefMacro,e_VG);
        
        e_VG.TIME = toc(ticSPEEDUP) ;
        %Como dentro de las funciones que se llama a partir de newton_raphson, Fint representa el
        %incremento de las fuerzas internas, para impresi�n de los resultados se actualiza las fuerzas
        %internas.
        %Fint = Fint+DFint;
        % Fint_old =  Fint_new;
        fprintf('Tiempo del Newton: %f\n',toc(ticIDNewt));

        % IMPRESION DE RESULTADOS
        index_print = rem(istep,e_VG.IRES);
        if (index_print == 0)
            if mod(istep,postpro_impre_step)==0
                % DESPLAZAMIENTO TOTAL DE LA MICRO-CELDA
                uTotal = [DefMacroT(1),DefMacroT(4)/2;DefMacroT(4)/2,DefMacroT(2)]*xx(:,1:2)'...
                    +reshape(u,ndime,[]);
                matlab2gid_res(istep,iStage,in,u,u_old,c_GdlCond,e_DatSet,e_VarEst_new,e_VarEst_old,e_VarAux,DefMacro,uTotal,...
                    ELOCCalc,e_VG)
            end
        end
        
        % Operaciones Constitutivas despu�s de la convergencia del Newton
        %Se coloca al final de todo, despu�s de la impresi�n de los resultados para que se imprima los
        %valores con los datos que se utilizaron para obtenerlos. Por ejemplo, que las tensiones que se
        %grafica se corresponde con la normal indicada, y no la que se podr�a obtener del an�lisis dentro
        %de f_OperConst.
        [e_DatSet,e_VarEst_new,e_VarAux,sigmaHomog,e_VG] = ...
            f_OperPosConv(u,xx,e_VarEst_new,e_VarEst_old,e_VarAux,e_DatSet,c_CT,KT,m_LinCond,dofl,doff,isMICRO,e_VG);
        
        % Almacenamiento de datos para graficos X-Y
        f_SaveDatGraf(u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,e_DatSet,m_SetElem,sigmaHomog,epsilon_Macro,e_VG)        
        
        % ACTUALIZACION DE VARIABLES
        psi_value_old           = psi_value;
        Du_step_old             = Du_step_new;
        e_VarEst_old            = e_VarEst_new;
        u_old                   = u;
        if (CONTROL_STRAT == 4)
            e_VG.lambda(iStage) = lambda;
        end
               
        fprintf('FIN DE PASO: %-3d (tiempo: %f)\n',istep,toc(ticStep));
        disp('*******************************************************************')
        
    end
    
    if TR_VAR; break; end
    
    % revision de grados de libertad
    dofl_old = dofl;
    u_old_Stage = u;
end

OTHEROUT.tiempo_TOTAL_HF = toc(ticStep);

end
