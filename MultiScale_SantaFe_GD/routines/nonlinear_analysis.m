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

% VARIABLES GLOBALES
struhyp = e_VG.struhyp;
ndoft              = e_VG.ndoft;
nElem              = e_VG.nElem;
ndn = e_VG.ndn;
np                 = e_VG.np;
postpro_impre_step = e_VG.postpro_impre_step;
CONTROL_STRAT      = e_VG.CONTROL_STRAT;
Dtime              = e_VG.Dtime;
ndime              = e_VG.ndime;
ntens              = e_VG.ntens;
nSet = e_VG.nSet;

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
%Por simplicidad se considera una deformaciï¿½n macro aplicada distinta por elemento, y no por PG.
%(no es necesario considerar una estructura para esta variable)
%DefMacro = zeros(ntens*npg,nElem);
%DefMacro = zeros(ntens,nElem);
% %% Iniciacion Matriz de Snapshots (JLM)
% nset = e_DatSet.e_DatMat.e_VG.nSet ;
% ntens= e_DatSet.e_DatMat.e_VG.ntens ;
% nglT = 0 ;
% for iset=1:nset
%     nElem = e_DatSet.e_DatMat.e_DatSet(iset).nElem ;
%     npg = e_DatSet.e_DatMat.e_DatSet(iset).e_DatElem.npg ;
%     nglT = nglT + nElem*npg;%*ntens ;
% end
% SnapStrain = zeros(nglT*ntens,np) ;
% SnapEnergy = zeros(nglT,np) ;
% Snapflag       = zeros(1,np) ;
% %%
% Only for ROM purposes
EPS_ROMI=[]; eta=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DefMacro = arrayfun(@(x)zeros(e_VG.ntens,x.e_DatElem.npg,x.nElem),e_DatSet,'UniformOutput',false);
% if isMICRO.MICRO
%     epsilon_Macro0 = isMICRO.epsilon_Macro0 ; % Increment of macro-strain
%     epsilon_Macro  = isMICRO.epsilon_Macro ; % Initial imposed macro-strain
% end

if e_VG.isMICRO.MICRO
    switch struhyp %e_VG.conshyp
        case 20 %{100,110} % No linealidad geometrica Elastic Material neo-Hookean y J2, It is assumed next order Voigt tensor [Fxx;Fyy;Fzz;Fxy;Fyx]
            epsilon_Macro  = [1; 1; 1; 0; 0;];
            epsilon_Macro0 = e_VG.isMICRO.epsilon_Macro0 ; % + epsilon_Macro +  Increment of macro-strain
        otherwise % Pequeñas deformaciones
            epsilon_Macro  = zeros(ntens,1); % Initial imposed macro-strain
            epsilon_Macro0 = e_VG.isMICRO.epsilon_Macro0 ; % Increment of macro-strain
    end 
% Only for ROM purposes
    if e_VG.MOD_TYPE~=0
        EPS_ROMI=zeros(e_VG.nModesEPS_TOT,1); eta=zeros(e_VG.ntens-1,1);
    %else
    %    EPS_ROMI=[]; eta=[];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5

% switch struhyp
%    case {1,2,3,4,5}
%       DefMacroT = zeros(ntens,nElem);
%       %DefMacro = arrayfun(@(x)zeros(e_VG.ntens,x.e_DatElem.npg,x.nElem),e_DatSet,'UniformOutput',false);
%       %DefMacro = zeros(ntens,nElem);
%       DefMacro = arrayfun(@(x)zeros(ntens,x.nElem),e_DatSet,'UniformOutput',false);
%    case 20
%       %It is assumed next order Voigt tensor [Fxx;Fyy;Fzz;Fxy;Fyx].
%       DefMacroT = zeros(ntens,nElem);
%       %DefMacro = [ones(ndime,nElem);zeros(ntens-ndime,nElem)];
%       DefMacro = arrayfun(@(x)[ones(3,x.nElem);zeros(ntens-3,x.nElem)],e_DatSet,...
%          'UniformOutput',false);
%    otherwise
%       error(['Nonlinear analysis: Macro information variables definition: ',...
%          'Structural Hypothesis has not been implemented.'])
% end

% eps_fluct_old   = zeros(sitvare,nElem);       % Vector de deformaciones fluctuantes previas
% e_VarEst_old = f_eVarEstInic({'sigma','eps','eps_fluct','hvar'},e_DatSet,e_VG);
% %e_VarEst_new = e_VarEst_old;
% %Se envía xx solo para inicializar los matriz de deformación del punto central del elemento FBar_LD.
% e_VarAux = f_eVarAuxInic(xx,e_DatSet,e_VG);

ELOCCalc = false(1,nElem);
if e_VG.isMICRO.MICRO
    e_VarEst_old = f_eVarEstInic({'sigma','eps','eps_fluct','hvar','delta_eps'},e_VG.ROM_II.e_DatSet_ROM,1,e_VG.MOD_TYPE); %JLM
    %Se envía xx solo para inicializar los matriz de deformación del punto central del elemento FBar_LD.
    e_VarAux = f_eVarAuxInic(e_VG.ROM_II.e_DatSet_ROM,1,xx,e_VG);% VER!!! e_DatMatSet = e_DatSet(iSet).e_DatMat;
elseif(~e_VG.esME && e_VG.MOD_TYPE==2)   %JLM antes era ~e_VG.esME
%     e_VarEst_old =
%     f_eVarEstInic({'sigma','eps','eps_fluct','hvar','delta_eps'},e_VG.ROM_II.e_DatSet_ROM,1);     % Original
%     e_VarEst_old = f_eVarEstInic({'sigma','eps','eps_fluct','hvar','delta_eps'},e_VG.ROM_II.e_DatSet_ROM,1,e_VG.MOD_TYPE); %JLM para correr HPROM con isMICRO=1
    e_VarEst_old = f_eVarEstInic({'sigma','eps','eps_fluct','hvar','delta_eps'},e_DatSet,nSet,e_VG.MOD_TYPE); % e_DatSet.e_DatMat.ROM_II.e_DatSet_ROM,1,e_VG.MOD_TYPE); %JLM para correr HPROM en MS    
    %Se envía xx solo para inicializar los matriz de deformación del punto central del elemento FBar_LD.
%     e_VarAux = f_eVarAuxInic(e_VG.ROM_II.e_DatSet_ROM,1,xx,e_VG); %JLM para correr HPROM con isMICRO=1
    e_VarAux = f_eVarAuxInic(e_DatSet,nSet,xx,e_VG);%e_DatSet.e_DatMat.ROM_II.e_DatSet_ROM,1,xx,e_VG); %JLM para correr HPROM en MS  
% elseif e_VG.isMICRO.MICRO    
%     e_VarEst_old = f_eVarEstInic({'sigma','eps','eps_fluct','hvar'},e_DatSet,nSet,e_VG.MOD_TYPE);    %Se envía xx solo para inicializar los matriz de deformación del punto central del elemento FBar_LD.
% %     e_VarAux = f_eVarAuxInic(e_DatSet,nSet,xx);% VER!!! e_DatMatSet = e_DatSet(iSet).e_DatMat;
%     e_VarAux = f_eVarAuxInic(e_DatSet,nSet,xx,e_VG);% VER!!! e_DatMatSet = e_DatSet(iSet).e_DatMat;
else
    e_VarEst_old = f_eVarEstInic({'sigma','eps','eps_fluct','hvar'},e_DatSet,nSet,e_VG.MOD_TYPE);
    %Se envía xx solo para inicializar los matriz de deformación del punto central del elemento FBar_LD.
    e_VarAux = f_eVarAuxInic(e_DatSet,nSet,e_DatSet.e_DatMat.xx,e_VG);%,e_DatSet.e_DatMat); % VER!!! e_DatMatSet = e_DatSet(iSet).e_DatMat;
end 
% f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet,e_DatMatSet)
c_GdlCond = f_cVarCondInic(e_DatSet,e_VG);

%Impresiï¿½n de archivo de postprocesado de la malla y inicializaciï¿½n del archivo de datos
matlab2gid_mesh(in,xx,e_DatSet,e_VG)
e_VG.istep = 0;
f_InicArchDat(in,m_SetElem,e_DatSet,e_VG)

%Ploteo en el paso 0 (se asume que tiene todos valores nulos en tiempo 0)
%Esto ï¿½ltimo habrï¿½a que ver con las condiciones de borde y la funciï¿½n psi_value (por ejemplo que esta no
%sea nula en el tiempo 0)
%sitvare = e_VG.sitvare;
%nnod    = e_VG.nnod;
%matlab2gid_res(0,u,zeros(sitvare,nElem),zeros(sitvare,nElem),zeros(sitvare,nElem),DefMacro,...
%   zeros(sihvare,nElem),e_VG,zeros(nnod,ntens),zeros(ndime,nnod));
f_SaveDatGraf(u,c_GdlCond,Fint,e_VarEst_old,e_VarAux,e_DatSet,m_SetElem,sigmaHomog,epsilon_Macro,e_VG) %JLM

%INTEGRACION TEMPORAL POR STAGE
for iStage = 1:e_VG.nStage
    
    % Variables correspondientes a cada Stage
    npStage    = np(iStage);
    e_VG.Dtime = Dtime(iStage);
    e_VG.iStage = iStage;
    
    % VALOR DE LA FUNCION Psi EN EL TIEMPO CERO
%     psi_value_old      = get_psi_value(funbc,0,e_VG);
    psi_value_old      = get_psi_value_Stage(funbc.m_func{iStage},0,e_VG);
    %Fext               = f*psi_value_old;

    %Con istepSave se indica desde que paso se empieza a correr. Si lee el archivo mat, el valor nulo de
    %istepSave es pisado por el paso donde se guardï¿½ el workspace.
    istepSave = 0;
    [dirFile,nomFile] = fileparts(e_VG.fileCompleto);
    if 0
       %Hay algunas variables que tienen que recuperarse como vienen previamente al load del workspace, como el
       %path al archivo, ya que si se guarda en una computadora con directorios distintos a donde se hace el load
       %del workspace va tirar error posterior a la lectura. Además estas variables que se utilizan su valor
       %previo no tiene cambiar durante la corrida del programa, ya que justamente no se recuperan.
       %IGUAL FALLA PORQUE HAY QUE CAMBIAR EL E_VGMicro.fileCompleto, también falla en el handle a un función,
       %donde se guarda el path absoluto.
       %fileCompletoTemp = e_VG.fileCompleto;
       load(fullfile(dirFile,'PasoSalvado'))
       %e_VG.fileCompleto = fileCompletoTemp;
       fprintf('Se recuperó los valores de las variables del paso %d.\n',istepSave)
    end
    %Se indica cada cuÃ¡nto se salva los pasos
%     deltaPasoSave = 100;

% INFORMACION ADICIONAL
%Se lee un script de matlab. Esto permite cambiar algunas propiedades, matrices, etc., antes de entrar en el
%cï¿½lculo en forma rï¿½pida sin modificar demasiado el programa.
%Ocurre un error al usar usar run, no sï¿½ si pasa lo mismo con eval, ya que matlab no se da cuenta
%que un script o funciï¿½n fue modificada para precompilarla (usa la precompilaciï¿½n anterior). Esto
%hace que las modificaciones del script las ignora y usa por ejemplo valores de variable que son de
%la versiï¿½n del script anterior.
%Esto se arregla con clear all, pero eso puede traer muchos problemas, ademï¿½s que se desaprovecha
%las precompilaciones previas. Se borra solo la precompilaciï¿½n de la funciï¿½n (hay que usar el nombre
%de la funciï¿½n solo, sin camino).
    if exist([e_VG.fileCompleto,'.m'],'file')
       clear(nomFile)
       %run corre el script sin estar en el path o que en el directorio activo
       run([e_VG.fileCompleto,'.m'])
    end

% FUNCIÓN CONDICIONES DE BORDE
% [m_LinCond,vfix,m_InvCRR,doff,dofl,doffCondCte] = f_CondBord(e_VG,xx,e_DatSet,e_VG.m_ConecFront,e_VG.RotacionGlobal);
% Uso la funcion f_CondBord de PD hecha por Manuel
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
    
    
% stepLimInf = 270;
% facdPsiv = 1;

% Medición de tiempos en un cluster mediante parTicToc
    c_ParTT = cell(np,1);

% INTEGRACION TEMPORAL
    for istep = istepSave+1:npStage

       ticStep = tic;
       e_VG.istep = istep;
       time = Dtime*istep;
       fprintf('STEP: %-3d\n',istep);

       % VALOR ACTUAL DE LA FUNCION TEMPORAL PARA ESCALAR CONDICIONES DE BORDE
%        psi_value = get_psi_value(funbc,time,e_VG);
        psi_value = get_psi_value_Stage(funbc.m_func{iStage},time,e_VG);

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % epsilon_Macro = epsilon_Macro + epsilon_Macro0*(psi_value - psi_value_old);
       % DefMacro =arrayfun(@(x)repmat(epsilon_Macro,[1,x.e_DatElem.npg,x.nElem]),e_DatSet,'UniformOutput',false);
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
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
            switch struhyp %e_VG.conshyp
                case 20 %{100,110} % No linealidad geometrica Elastic Material neo-Hookean y J2
                    %It is assumed next order Voigt tensor [Fxx;Fyy;Fzz;Fxy;Fyx]
                    DefMacroT = zeros(ntens,nElem);
                    DefMacro = arrayfun(@(x)[ones(3,x.nElem);zeros(ntens-3,x.nElem)],e_DatSet,'UniformOutput',false);
%                     DefMacro = arrayfun(@(x)zeros(e_VG.ntens,x.e_DatElem.npg,x.nElem),e_DatSet,'UniformOutput',false);
%                     DefMacroT = zeros(ntens,1);
                otherwise % Pequeñas deformaciones
                    DefMacro = arrayfun(@(x)zeros(e_VG.ntens,x.e_DatElem.npg,x.nElem),e_DatSet,'UniformOutput',false);
                    DefMacroT = zeros(ntens,1);
            end            
        end        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
       % DESPLAZAMIENTO IMPUESTO
       vDeltaFixTemp = vfix;
       Du_step_new = zeros(ndoft,1);

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

       % FUERZA EXTERNA IMPUESTA
%        if (CONTROL_STRAT == 1)
%           Fext = f*psi_value;  f.f_Stages{1}
%           vDeltaFixTemp(doffCondCte) = vfix(doffCondCte)*(psi_value - psi_value_old);
%           u(doff)                    = u(doff) + m_InvCRR*vDeltaFixTemp(doff);
%           Du_step_new(doff)          = m_InvCRR*vDeltaFixTemp(doff);
%        else
%           Fext = f;
%           e_VG.vfix                  =  vfix;
%           e_VG.vfix_doff             =  m_InvCRR*vfix(doff);
%           Du_step_new(doff)          =  e_VG.vfix_doff*e_VG.lambda - u(doff) ;
%           u(doff)                    =  e_VG.vfix_doff*e_VG.lambda;
%        end

        % FUERZA EXTERNA IMPUESTA TOTAL DEL STAGE (VALORES ACUMULADOS)
        if (CONTROL_STRAT == 1)
%             Fext = f*psi_value;
%             vDeltaFixTemp(doffCondCte) = vfix(doffCondCte)*(psi_value - psi_value_old);
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
        
       ticIDNewt = tic;
       % ESQUEMA DE NEWTON-RAPHSON
%        [u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,~,lambda,o_Par] = newton_raphson(xx,m_LinCond,...
%           dofl,doff,u,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,e_VarEst_old,e_VarAux,e_DatSet,DefMacro,e_VG);

       [u,EPS_ROMI,eta,c_GdlCond,Fint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT,lambda,iter] = newton_raphson(xx,m_LinCond,...
            dofl,doff,u,EPS_ROMI,eta,u_old_Stage,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,Fext_inic,e_VarEst_old,e_VarAux,e_DatSet,DefMacro,e_VG);
  
      
       %Como dentro de las funciones que se llama a partir de newton_raphson, Fint representa el
       %incremento de las fuerzas internas, para impresiï¿½n de los resultados se actualiza las fuerzas
       %internas.
       %Fint = Fint+DFint;
       %Fint_old = Fint_new;
       fprintf('Tiempo del Newton: %f\n',toc(ticIDNewt));
       c_ParTT{istep} = [];%o_Par; % JLM

       % IMPRESIÓN DE RESULTADOS
       index_print = rem(istep,e_VG.IRES);
       if (index_print == 0)
          if mod(istep,postpro_impre_step)==0
             % DESPLAZAMIENTO TOTAL DE LA MICRO-CELDA
             uTotal = [DefMacroT(1),DefMacroT(4)/2;DefMacroT(4)/2,DefMacroT(2)]*xx(:,1:2)'...
                +reshape(u,ndn,[]);
%              matlab2gid_res(istep,iStage,in,u,u_old,c_GdlCond,e_DatSet,e_VarEst_new,e_VarEst_old,e_VarAux,DefMacro,uTotal,...
%                     ELOCCalc,e_VG)
            if istep==1 
                matlab2gid_PG_res(istep,iStage,in,u,u_old,c_GdlCond,e_DatSet,e_VarEst_new,e_VarEst_old,e_VarAux,DefMacro,uTotal,ELOCCalc,e_VG)
            end
          end
       end

       % Operaciones Constitutivas despuï¿½s de la convergencia del Newton
       %Se coloca al final de todo, despuï¿½s de la impresiï¿½n de los resultados para que se imprima los
       %valores con los datos que se utilizaron para obtenerlos. Por ejemplo, que las tensiones que se
       %grafica se corresponde con la normal indicada, y no la que se podrï¿½a obtener del anï¿½lisis dentro
       %de f_OperConst.
%        [e_DatSet,e_VarEst_new,e_VarAux,e_VG] = ...
%           f_OperPosConv(u,xx,e_VarEst_new,e_VarEst_old,e_VarAux,e_DatSet,c_CT,e_VG);

        [e_DatSet,e_VarEst_new,e_VarAux,sigmaHomog,e_VG] = ...
            f_OperPosConv(u,xx,e_VarEst_new,e_VarEst_old,e_VarAux,e_DatSet,c_CT,KT,m_LinCond,dofl,doff,e_VG);
        % Almacenamiento de datos para graficos X-Y
%         f_SaveDatGraf(u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,e_DatSet,m_SetElem,sigmaHomog,epsilon_Macro,e_VG)    ; % JLM    
       f_SaveDatGraf(u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,e_DatSet,m_SetElem,sigmaHomog,epsilon_Macro,e_VG) %JLM
                        
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
%    Strain_new = e_VarEst_new.hvar.e_VarEst(iset).eps_fluct ;
%    Strain_old = e_VarEst_old.hvar.e_VarEst(iset).eps_fluct ;
%    Stress_new = e_VarEst_new.hvar.e_VarEst(iset).sigma ;
%    Stress_old = e_VarEst_old.hvar.e_VarEst(iset).sigma ;
        % revision de grados de libertad
        dofl_old = dofl;
        u_old_Stage = u;

    end
OTHEROUT.tiempo_TOTAL_HF = toc(ticStep);
end


