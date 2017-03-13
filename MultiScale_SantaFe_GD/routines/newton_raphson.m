function [u,EPS_ROMI,eta,c_GdlCond,Fint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT,lambda,iter] = newton_raphson(xx,...
   m_LinCond,dofl,doff,u,EPS_ROMI,eta,u_old_Stage,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,Fext_inic,e_VarEst_old,e_VarAux,e_DatSet,...
   DefMacro,e_VG)

%******************************************************************************************
%*  ESQUEMA ITERATIVO NEWTON-RAPHSON                                                      *
%*                                                                                        *
%*  ARGUMENTOS DE ENTRADA:                                                                *
%*  inn         : lista de nodos                                                          *
%*  xx          : lista de coordenadas                                                    *
%*  conec       : lista de conectividades                                                 *
%*  u           : vector de desplazamientos                                               *                  
%*  Du_step_old : incremento de desplazamiento en el paso previo                          *
%*  Fint        : vector de fzas internas                                                 *
%*  Fext        : vector de fzas externa                                                  *
%*  hvar_old    : arreglo de variables historicas del paso previo                         *
%*  Eprop       : lista de propiedades de los elementos                                   *
%*                                                                                        *                  
%*  ARGUMENTOS DE SALIDA:                                                                 *                  
%*  u           : vector de desplazamientos                                               *                  
%*  Fint        : vector de fzas internas                                                 *
%*  sigma_new   : arreglo de tensiones del paso actual                                    *
%*  epsr_new    : arreglo de deformaciones del paso actual                                *
%*  hvar_new    : arreglo de variables historicas del paso actual                         *
%*  Du_step_new : incremento de desplazamiento en el paso actual                          *
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

% VARIABLES GLOBALES
ndoft = e_VG.ndoft;
tolnr = e_VG.tolnr;
iterMax = e_VG.iterMax;

if e_VG.CONTROL_STRAT ==4
    lambda = e_VG.lambda(e_VG.iStage);
    lambda_old = lambda;
    if e_VG.istep==1
        lambda_old=0;
    end
else
    lambda=1;
end

% INICIALIZACION DE VARIABLES
iter = 0;
%Se crea aca para evitar estar cre�ndola en cada iteraci�n del newton. Principalmente para el caso del modelo
%multiescala, donde hay que redifinir hvar como una estructura y no como una matriz como los otros casos. No
%interesa sus valores, solo que est� definida para indexar correctamente los argumentos de salida.
%Ver si no cambiar la estructura para que se pueda indicar que se necesita como variable de estado (sigma,
%eps, eps_fluct, etc.) por set. Esto evitar�a tener matrices que no se necesita en ciertos sets. Otra es tener
%variables que indiquen las dimensiones que se necesita de estas matrices, y en el caso que no se necesite que
%sea cero.
if e_VG.isMICRO.MICRO
    e_VarEst_new = f_eVarEstInic({'sigma','eps','eps_fluct','hvar','delta_eps'},e_VG.ROM_II.e_DatSet_ROM,1,e_VG.MOD_TYPE); %JLM
elseif (~e_VG.esME && e_VG.MOD_TYPE==2)
%     e_VarAux = f_eVarAuxInic(e_VG.ROM_II.e_DatSet_ROM,1,xx,e_VG); %JLM para correr HPROM con isMICRO=1
%     e_VarEst_new = f_eVarEstInic({'sigma','eps','eps_fluct','hvar','delta_eps'},e_VG.ROM_II.e_DatSet_ROM,1)
%     e_VarEst_new = f_eVarEstInic({'sigma','eps','eps_fluct','hvar','delta_eps'},e_VG.ROM_II.e_DatSet_ROM,1,e_VG.MOD_TYPE); %JLM para correr HPROM con isMICRO=1
    e_VarEst_new = f_eVarEstInic({'sigma','eps','eps_fluct','hvar','delta_eps'},e_DatSet,1,e_VG.MOD_TYPE) ;% e_DatSet.e_DatMat.ROM_II.e_DatSet_ROM,1,e_VG.MOD_TYPE); %JLM para correr HPROM en MS
else                                                                                       
    e_VarEst_new = f_eVarEstInic({'sigma','eps','eps_fluct','hvar','delta_eps'},e_DatSet,e_VG.nSet,e_VG.MOD_TYPE,0);
end
du_iter = zeros(ndoft,1);

% ESQUEMA ITERATIVO DE NEWTON
% Variables condensadas
c_GdlCond = f_cVar_initIncNR(e_DatSet,e_VG,c_GdlCond);

% MODIFICACI�N DEL VALOR INICIAL
%Parece que es muy lento, considerando que en multiescala se lo llama repetida veces (una vez por
%cada iteraci�n macro). Ver si se calcula esta condici�n inicial, no utilizar la iteraci�n cero y
%forzar que siempre se haga una iteraci�n m�s.
%ticIdVI = tic;
%u = f_ValorInicial(xx,m_LinCond,dofl,doff,u,c_GdlCond,Fext,e_VarEst_new,e_VarEst_old,...
%   e_VarAux,e_DatSet,DefMacro,e_VG);
%fprintf('Tiempo de c�lculo de valor inicial: %f\n',toc(ticIdVI));

% ITERACION "0", EVALUACION DE LA FUERZA INTERNA
e_VG.iter = iter;
%Para forzar la iteraci�n cero, se resuelva como un problema el�stico, independiente de
%como venga la variable e_VG.elast (luego en las siguientes iteraciones se recupera el valor 
%impuesto en esa variable). Esta resuelve el caso que durante la parte el�stica del problema, se
%itere dos veces ya que el valor inicial fuerza un comportamiento inel�stico ficticio.
%Notar que lo se consigue es realizar en la interaci�n 1 una avance el�stico, es decir 
%u(1) = u(0)+Du (u(i+1)=u(i)+Du), Du el�stico (cualquiera sea el estado actual del problema,
%que no corresponde al comportamiento de los modelos constitutivos elegidos) y que el residuo va
%corresponder a ese comportamiento el�stico asumido (es decir va ser incorrecto).
% elastOrig = e_VG.elast;
% e_VG.elast = 1;
%
%ticIDEnsam = tic;
[KT,Fint,c_GdlCond,e_VarEst_new,e_VarAux,c_CT] = f_Par_MatGlobales(xx,u,EPS_ROMI,eta,Du_step_new,c_GdlCond,e_DatSet,...
  e_VarEst_new,e_VarEst_old,e_VarAux,DefMacro,e_VG);
%fprintf('Tiempo de ensamblaje: %f\n',toc(ticIDEnsam));
% e_VG.elast = elastOrig;

% ITERACION "0", COMPUTO DEL RESIDUO
if (~e_VG.esME && (e_VG.MOD_TYPE==1 || e_VG.MOD_TYPE==2) && e_VG.isMICRO.MICRO==1) % JLM ~e_VG.esME && (e_VG.MOD_TYPE==1 || e_VG.MOD_TYPE==2)
    res = Fint;
else
    res = residuo(Fext,Fext_inic,Fint,m_LinCond,dofl,doff,lambda);
end
norm_res = norm([res;cell2mat((cellfun(@(x)reshape(x,[],1),c_GdlCond(:,2),'UniformOutPut',0)))],inf);

% IMPRESION de interaci�n 0
if ~e_VG.esME
   fprintf('ITERACION: %2d ; |res| = %-8.3e (VALOR INICIAL)\n',iter,norm_res);
else
   fprintf('Problema micro: Elemento Macro %d: Punto gauss macro %d:\n',e_VG.iElemNumMacro,e_VG.iPGMacro);
   fprintf('ITERACION: %2d ; |res| = %-8.3e (VALOR INICIAL)\n',iter,norm_res);
end
% INICIALIZACION DE VARIABLES PROPIAS DE ESQUEMA NEWTON-RAPHSON
norm_res = 1;
tol_rel  = tolnr;
max_iter = iterMax;

% LOOP DE NEWTON
while ((norm_res > tol_rel || iter < 1) && (iter <= max_iter))
   
   % CONTEO DE ITERACIONES
   iter = iter+1;
   e_VG.iter = iter;
   
   % INCREMENTO DE DESPLAZAMIENTO S/ESTRATEGIA DE CONTROL
   %ticIDResSist = tic;
   [du_iterl,dlambda] = incremental_disp(res,KT,Fext,m_LinCond,dofl,doff,Du_step_new,Du_step_old,e_VG);
   
   if ~e_VG.esME && (e_VG.MOD_TYPE==1) && e_VG.isMICRO.MICRO==1 % FIRST REDUCTION  JLM
       % ACTUALIZACION VARIABLES NR
       EPS_ROMI = EPS_ROMI + du_iterl(1:e_VG.nModesEPS_TOT,1) ;
       eta = eta + du_iterl(e_VG.nModesEPS_TOT+1:end,1) ;
       
   elseif ~e_VG.esME && (e_VG.MOD_TYPE==2) && e_VG.isMICRO.MICRO==1 % SECOND REDUCTION JLM
       % ACTUALIZACION VARIABLES NR
       EPS_ROMI = EPS_ROMI + du_iterl(1:e_VG.nModesEPS_TOT,1) ;
       eta = eta + du_iterl(e_VG.nModesEPS_TOT+1:e_VG.nModesEPS_TOT+e_VG.ntens-1,1) ;       
       %Q_ROMI = Q_ROMI + du_iterl(e_VG.nModesEPS_TOT+e_VG.ntens:end,1) ;
       
   else
       du_iter(dofl) = du_iterl;
       du_iter(doff) = m_LinCond*du_iterl;
       % INCREMENTO DE VARIABLES CONDENSADAS (gdl internos)
       c_GdlCond = f_cVar_IncrNR(du_iter,c_GdlCond,e_DatSet,e_VarEst_new,e_VarAux,e_VG);
       u = u + du_iter;
   end  
   
   %Du_step_new(dofl) = Du_step_new(dofl) + du_iter;
   %Du_step_new(doff) = Du_step_new(doff) + m_LinCond*du_iter;
   Du_step_new = Du_step_new+du_iter;
   
   lambda = lambda + dlambda;
   
   if e_VG.CONTROL_STRAT ==4
       u(doff) = u_old_Stage(doff) + e_VG.vfix_doff*lambda;
       Du_step_new(doff) =  e_VG.vfix_doff*(lambda-lambda_old);
   end
   
   % VARIABLES CONDESADAS TOTALES
   % c_GdlCond(:,1) = cellfun(@plus,c_GdlCond(:,1),c_dGdlCond,'UniformOutput',false);
   
   % MATRIZ DE RIGIDEZ y COMPUTO DE LA FUERZA INTERNA
   %ticIDEnsam = tic;
   [KT,Fint,c_GdlCond,e_VarEst_new,e_VarAux,c_CT] = f_Par_MatGlobales(xx,u,EPS_ROMI,eta,Du_step_new,c_GdlCond,e_DatSet,...
      e_VarEst_new,e_VarEst_old,e_VarAux,DefMacro,e_VG);
   %fprintf('Tiempo de ensamblaje: %f\n',toc(ticIDEnsam));
   
   % EVALUACION DEL RESIDUO
   %res = residuo(Fext,Fext_inic,Fint,m_LinCond,dofl,doff,lambda);
    if ~e_VG.esME && (e_VG.MOD_TYPE==1 || e_VG.MOD_TYPE==2) && e_VG.isMICRO.MICRO==1 % JLM
      res = Fint;
   else
       res = residuo(Fext,Fext_inic,Fint,m_LinCond,dofl,doff,lambda);
   end
   
   % NORMA DEL RESIDUO (se considera los residuos de los gdl sin condensar)
   norm_res = norm(res)  +  f_cVar_normRes(e_DatSet,e_VG, c_GdlCond);

 %  norm_res = norm([res;cell2mat((cellfun(@(x)reshape(x,[],1),c_GdlCond(:,2),...
 %     'UniformOutPut',0)))],inf);
   
   % IMPRESION POR PANTALLA DEL PROCESO DE CONVERGENCIA
    if ~e_VG.esME
      fprintf('ITERACION: %2d ; |res| = %-8.3e\n',iter,norm_res);
    end
end

% Weak problem (recovering displacements from strains)
if (~e_VG.esME && (e_VG.MOD_TYPE==1 || e_VG.MOD_TYPE==2) && e_VG.isMICRO.MICRO==1) && e_VG.INV_PROB % JLM
    if e_VG.MOD_TYPE==1
        MAT_WEAK = e_VG.MAT_WEAK;
        ValRM = e_VG.ValRM;
    else
        MAT_WEAK = e_VG.ROM_II.MAT_WEAK;
        ValRM = e_VG.ROM_II.ValRM;
    end
    u(ValRM)=MAT_WEAK*EPS_ROMI;
end

% IMPRESION POR PANTALLA DEL PROCESO DE CONVERGENCIA
if (iter > max_iter)
   %fprintf('ESQUEMA DE NEWTON NO CONVERGE\n');
   %pause
   %En el debug se puede hacer que se active un breakpoint si usa un warning como mensaje.
   warning('off','backtrace')
   if e_VG.esME
      %En el caso que se paralize en varias m�quinas esto va fallar, ya que no va encontrar el directorio
      %donde est� el archivo
%       fileCompleto = [e_VG.fileCompleto,'_EM',int2str(e_VG.iElemNumMacro),'PGM',int2str(e_VG.iPGMacro)];
%       nomArch = [fileCompleto,'.pasNoConv'];
%       fId = fopen(nomArch,'at');
%       fprintf(fId,'Paso %d - IterMacro %d - Iter %d: |res| = %e\n',e_VG.istep,e_VG.iterMacro,...
%          iter,norm_res);
%       fclose(fId);
      if e_VG.nLab==0||strcmp(e_VG.tipoMPool,'local')
         nomArch = [e_VG.fileCompleto,'.pasNoConv'];
         fId = fopen(nomArch,'at');
         fprintf(fId,'Paso %d - IterMacro %d - ElemMacro %d - PGMacro %d - Iter %d: |res| = %e\n',e_VG.istep,...
            e_VG.iterMacro,e_VG.iElemNumMacro,e_VG.iPGMacro,iter,norm_res);
         fclose(fId);
      end
      %
      warning(['Problema micro: Elemento Macro %d: Punto de gauss macro %d: Esquema de newton no converge ',...
         '(|res|=%e).'],e_VG.iElemNumMacro,e_VG.iPGMacro,norm_res);
   else
      nomArch = [e_VG.fileCompleto,'.pasNoConv'];
      fId = fopen(nomArch,'at');
      fprintf(fId,'Paso %d - Iter %d: |res| = %e\n',e_VG.istep,iter,norm_res);
      fclose(fId);
      warning('ESQUEMA DE NEWTON NO CONVERGE'); %#ok<WNTAG>
   end
   warning('on','backtrace')
else
   if ~e_VG.esME
      fprintf('CONVERGE\n');
   else
      fprintf('Problema micro: Converge\n')
   end
   %fprintf('----------------------------------\n');
end