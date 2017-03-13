function [T_TOTAL,OTHEROUT] = analysis(file,path_file,nLabImp,isMICRO,MOD_TYPE,varargin)

%******************************************************************************************
%*  PROGRAMA DE ELEMENTOS FINITOS PARA LA RESOLUCION DE PROBLEMAS EN MECANICA DE SOLIDOS  *
%*                                                                                        *                  
%*  MATERIAL BASICO DEL CURSO DE POSTGRADO:                                               *                  
%*  MODELOS CONSTITUTIVOS PARA MATERIALES DISIPATIVOS. APLICACION A MECANICA DE SOLIDOS   *                  
%*                                                                                        *                  
%*  DOCTORADO EN INGENIERIA - MENCION MECANICA COMPUTACIONAL. FICH-UNL                    *                  
%*                                                                                        *                  
%*  MODELOS CONSTITUTIVOS IMPLEMENTADOS                                                   *
%*    1) Elasticidad lineal                                                               *
%*                                                                                        *                  
%*  MODELOS CONSTITUTIVOS A IMPLEMENTAR EN EL DESARROLLO DEL CURSO                        *
%*    2) Elasto-plasticidad (teoria J2) con endurecimiento isotropo                       *
%*    3) Visco-plasticidad (teoria J2) con endurecimiento isotropo                        *
%*    4) Da�o isotropo                                                                    *
%*    5) Da�o isotropo solo en traccion                                                   *
%*    6) Visco-danio isotropo                                                              *
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

%Para separar cuando est� corriendo el problema macro del problema micro, se pinta de colores
%distinto el lado izquierdo de la l�nea de comando seg�n en que parte est�. Utilizar para debug.
% cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
% listeners = cmdWinDoc.getDocumentListeners;
% jFxCommandArea = listeners(3);
% set(jFxCommandArea,'Background','yellow');

% ****************************
% * CONFIGURACION DE ENTORNO *
% ****************************
%Para cerrar el MatLab Pool cuando termina de correr.
cerMPool = 0;

if 0
%Definicion de los argumentos de entrada
switch nargin 
   case 0
      %Esto es el caso que no se pasa argumento, donde se se abre el di�logo para abrir el archivo
      %y se utiliza el n�mero de laboratorios por defecto.
   case 1
      if ischar(varargin{1})
         %Si se pasa como argumento el nombre del archivo con el directorio completo o solo el 
         %nombre si el mismo est� en el path del matlab.
         %analysis(file&Dir)
         [path_file,file,ext] = fileparts(varargin{1});
         file = [file,ext];
      else
         %Si el argumento es un n�mero, lo toma como el n�mero de lab a activar en matlabpool.
         %Ver si no verificar el argumento que sea num�rico en este caso.
         %analysis(nLabImp)
         nLabImp = varargin{1};         
      end
   case 2
      %Se asume que el primer argumento es un char, indicando el nombre o directorio completo, o
      %solo el directorio.
      if ischar(varargin{2})
         %En el caso que el 2do argumento sea el nombre del archivo.
         %analysis(file,path_file)         
         path_file = varargin{1};
         file = varargin{2};
      else
         %Si el segundo argumento es el n�mero de laborarios (y el primero es el nombre y directorio
         %del archivo)
         %analysis(file&Dir,nLabImp)
         [path_file,file,ext] = fileparts(varargin{1});
         file = [file,ext];
         nLabImp = varargin{2};   
      end      
   case 3
      %En el caso que se pase los tres argumentos
      %analysis(file,path_file,nLabImp)
      path_file = varargin{1};
      file = varargin{2};
      nLabImp = varargin{3};
   case 4
      %the case of 4 arguments including the microscale flag
      %analysis(file,path_file,nLabImp)
      path_file = varargin{1};
      file = varargin{2};
      nLabImp = varargin{3};
      isMICRO = varargin{4};
   otherwise
      error('Analyis: Numero de argumentos de llamada incorrectos.')
end
end

%INICIALIZACION DE LAS VARIABLES
INIT_VARS;

%En el caso que no se defina el directorio y el nombre del archivo, se busca por medio de un
%dialogo.
if ~exist('file','var')&&~exist('path_file','var')
    [file,path_file] = uigetfile('*.mfl','Definir archivo de calculo');
end

%Si se quiere imponer una cierta cantidad de laboratorios a crear, se pone un valor distinto mayor
%que cero. Con cero se utiliza el valor por defecto del MatLab o si hay un matlabpool abierto, lo
%utiliza (aunque los labs no sean los por defecto). Un n�mero negativo hace que no se active el
%matlabpool (si est� abierto, lo cierra) (esto sirve para hacer debug dentro de las funciones que
%llama el parfor, pero no en el mismo parfor, para hacer ello hay que cambiar por el for).
if ~exist('nLabImp','var')
   nLabImp = 0;
end

%clear; 
%close all
%Es importante la siguiente linea para que cuando se corta la ejecuci�n del programa antes de cerrar
%todos los archivos.
fclose all;
clc
%p = pwd;
p = fileparts(mfilename('fullpath'));
addpath(genpath(p));

% ********************
% * LECTURA DE DATOS *
% ********************
if ischar(file)&&ischar(path_file)
   
   % Lectura de datos del problema
   [in,xx,m_SetElem,f,funbc,e_DatSet,e_VG] = read_data(file,path_file);
          
   ADD_CONFIGS;
       
   % *************************************
   % * RESOLUCION DEL PROBLEMA NO LINEAL *
   % *************************************
%    nLab = matlabpool('size');
%    %Se guarda el objeto scheduler para la configuraci�n activa para ver si la corrida es de tipo local,
%    %as� no se realiza la transferencia de archivos.
%    objSch = findResource();
%    if (nLabImp==0||nLabImp==nLab)&&nLab>0
%       disp(['Se utiliza la configuraci�n activa: ',get(findResource(),'configuration'),...
%          ' de ',num2str(nLab),' laboratorios.']);
%    else
%       %Esto tambi�n cierra el MatLab si nLabImp==-1
%       if nLab>0&&nLab~=nLabImp
%          matlabpool close
%       end
% 
%       %Recordar que para pocos elementos no conviene paralelizar.
%       if nLabImp==0
%          if strcmp(objSch.type,'local')
%             matlabpool open
%          else
%             matlabpool('open','FileDependencies',{p})
%          end
%       elseif nLabImp>0
%          if strcmp(objSch.type,'local')
%             matlabpool('open',nLabImp)
%          else
%             matlabpool('open',nLabImp,'FileDependencies',{p})
%          end
%       end
%    end
%    %
%    nLab = matlabpool('size');
%    if nLab>0&&~strcmp(objSch.type,'local')
%       %Para que funcione en un cluster hay que indicar que archivos se debe copiar (se indica el
%       %directorio del programa)
%       %matlabpool('addfiledependencies',{p})
%       %El directorio por defecto no es el directorio donde copia los archivos, sino que otro dado
%       %por pctRunOnAll getFileDependencyDir (directorio donde descomprime los archivos de la
%       %dependencia en cada nodo), por lo que en cada nodo se debe modificar para que el directorio
%       %activo (pwd) sea el del programa y luego generar los caminos.
%       stringCmd = 'pDep=getFileDependencyDir;if ~isempty(pDep),cd(getFileDependencyDir),end;';
%       pctRunOnAll(stringCmd)
%       %pctRunOnAll pwd
%       pctRunOnAll('addpath(genpath(pwd))')
%       %pctRunOnAll path
%       %Para asegurar que se tenga las �ltimas copias del c�digo (si no se hace esto cuando se
%       %modifica un archivo y el pool se deja abierto, se ignora los cambios, verificar esto)
%       matlabpool updatefiledependencies
%    end

   %Se almacena los n�meros de laborarios y el tipo de conexi�n para usar imprimir en el newton en un archivo
   %los pasos que no convergen (en el caso que es un cluster tira un error porque no puede crear el archivo en
   %en los nodos, en el camino indicado).
%    e_VG.nLab = nLab;
%    e_VG.tipoMPool = objSch.type;
   
   e_VG.nLab = 1;
   e_VG.tipoMPool = 'local'; 
   e_VG.isMICRO = isMICRO ;
   e_VG.INV_PROB = INV_PROB ;
   
   ticIDNLA = tic;
   [u,OTHEROUT] = nonlinear_analysis(in,xx,m_SetElem,f,funbc,e_DatSet,e_VG,isMICRO,...
       'SEL_SNPSH_DISP',SEL_SNPSH_DISP,'INDEX_INV_VAR_PLOT',INDEX_INV_VAR_PLOT,...
       'INDEX_VAR_STATE_STEP',INDEX_VAR_STATE_STEP);
   
    %if OTHEROUT.CONVERGE == 0
    %    warning('Not-conveged step...')
    %    return
    %end
    T_TOTAL = toc(ticIDNLA);
    
    fprintf('**************************************************************\n')
    fprintf('Tiempo de calculo de nonlinear analysis: %f seg.\n',T_TOTAL)
    fprintf('**************************************************************\n')
   
   if cerMPool
      matlabpool close
   end
   
end

end
