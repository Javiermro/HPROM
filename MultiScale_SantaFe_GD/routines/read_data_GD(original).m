function [in,xx,m_SetElem,f,funbc,e_DatSet,e_VG] = read_data(file,path_file)

%******************************************************************************************
%*  FUNCION DE LECTURA DE DATOS Y GENERACION DE VARIABLES PRINCIPALES                     *
%*                                                                                        *
%*  ARGUMENTOS DE ENTRADA:                                                                *
%*  file      : nombre del archivo de datos                                               *
%*  path_file : path completo del archivo de datos                                        *
%*                                                                                        *
%*  ARGUMENTOS DE SALIDA:                                                                 *
%*  in        : lista de nodos                                                            *
%*  inn       : lista de nodos                                                            *
%*  xx        : lista de coordenadas                                                      *
%*  iel       : lista de elementos                                                        *
%*  ieln      : lista de elementos                                                        *
%*  conec     : lista de conectividades                                                   *
%*  vfix      : vector de desplazamientos impuestos                                       *
%*  f         : vector de cargas externas aplicadas                                       *
%*  funbc     : funcion temporal para aplicar condiciones de borde                        *
%*  Eprop     : lista de propiedades de los elementos                                     *
%*  flags     : indicador de error de lectura de datos                                    *
%*                                                                                        *
%*  A.E. Huespe, P.J.Sanchez                                                              *
%*  CIMEC-INTEC-UNL-CONICET                                                               *
%******************************************************************************************

% Variable que almacena si el problema se puede considerar sim�trico
%(el valor por defecto es que sea un problema sim�trico, y seg�n el modelo constitutivo elemento o
%directamente hip�tesis estructural, si algunos de ellos hacen el problema no sim�trico, se la modifica)
esKTSim = 1;

%*********************************
%* APERTURA DEL ARCHIVO DE DATOS *
%*********************************
file_open = fullfile(path_file,file);
fid = fopen(file_open,'rt');
[~,file,ext] = fileparts(file);

%*********************************
%* LECTURA DEL TIPO DE PROBLEMA  *
%*********************************
%Por ahora no se utiliza para nada.
seccion = f_ProxString(fid);
f_VerifNom(seccion,'CSDA','Tipo de problema: Se debe ingresar "CSDA".')
tipoProbl = [f_ProxString(fid),' ',f_ProxString(fid)];

%***************************************
%* LECTURA DE C�MO EJECUTAR EL PROBLEMA  *
%***************************************
%Por ahora no se utiliza para nada.
seccion = f_ProxString(fid);
f_VerifNom(seccion,'START','Definici�n de tipo de ejecuci�n: Se debe ingresar "START".')
tipoCorr = f_ProxString(fid);

%************************************************
%* LECTURA DE VARIABLES INICIALES DEL PROBLEMA  *
%************************************************
seccion = f_ProxString(fid);
f_VerifNom(seccion,'CONTROL_DATA',...
   'Variables iniciales: No est� definido el inicio con "CONTROL_DATA".')
seccion = f_ProxString(fid);
f_VerifNom(seccion,'GEOMETRY','Variables iniciales: No est� definido "GEOMETRY".')
struhyp = textscan(fid,'%f',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
struhyp = struhyp{1};
seccion = f_ProxString(fid);
f_VerifNom(seccion,'DIMENSIONS','Variables iniciales: No est� definido "DIMENSIONS".')
c_ValVarInic = f_ExtrVar(fid);
seccion = f_ProxString(fid);
if strcmp(seccion,'ADD_PARAMETERS')
   c_ValVarInic2 = f_ExtrVar(fid);
   seccion = f_ProxString(fid);
end
f_VerifNom(seccion,'END_CONTROL_DATA',...
   'Variables iniciales: Se debe terminar con "END_CONTROL_DATA".')

%*******************************************************************************
%* HIPOTESIS ESTRUCTURAL:                                                      *
%*   1 = Estado plano de deformacion                                           *
%*   2 = Estado plano de tension                                               *
%*   3 = Estado tridimensional                                                 *
%*   4 = Axisimetria                                                           *
%*   5 = Barras en 2D                                                          *
%*******************************************************************************
switch struhyp
   %% SMALL DEFORMATIONS (Infinitesimal Strain)
   case {1,2}
      ndime = 2;
      ndn = 2;
      %ntens por ahora se va considerar como una variable global, dependiente de la hipot�sis
      %estructural y no del tipo de elemento.
      ntens = 4;
   case 3
      ndime = 3;
      ndn = 3;
   case 4
      ndime = 2;
      error('Lectura de datos: Hip�tesis estructural no completada')
   case 5
      ndime = 2;
      ndn = 2;
      ntens = 1;
   %% LARGE DEFORMATIONS (Finite Strain)
   case 20   %Plane Deformations State (2D).
      %The adopted measures are, for the stresses, the First Piola�Kirchhoff tensor (P) and for the strains,
      %the Deformation gradient tensor (F).
      ndime = 2;
      ndn = 2;
      %Because the non symmetry tensors. Voigt Notation Tensor Order is [xx,yy,zz,xy,yx].
      ntens = 5;
      %Como este problema se resuelve el equilibrio con el tensores no sim�tricos, se impone directamente
      %la falta de simetr�a del problema
      esKTSim = 0;
   otherwise
      fclose(fid);
      error('Lectura de datos: Hip�tesis estructural no implementada')
end

%*******************************************************************************
%* NUMERO DE NODOS                                                             *
%*******************************************************************************
nnod = f_ValVar('NPOIN',true,c_ValVarInic,'Variables iniciales');

%*******************************************************************************
%* NUMERO DE ELEMENTOS                                                         *
%*******************************************************************************
nElem = f_ValVar('NELEM',true,c_ValVarInic,'Variables iniciales');

%*******************************************************************************
%* ELEMENTOS DE FRONTERA                                                       *
%*******************************************************************************
nElemFront = f_ValVar('NELFR',false,c_ValVarInic,'Variables iniciales');
%Se considera opcional la cantidad de elementos de frontera, si no se ingresa se asume que es nula.
if isnan(nElemFront)
   nElemFront = 0;
end

%*******************************************************************************
%* NUMERO DE NODOS POR ELEMENTO                                                *
%*******************************************************************************
npe = f_ValVar('NDPEL',true,c_ValVarInic,'Variables iniciales');

%*******************************************************************************
%* NUMERO DE PUNTOS DE GAUSS                                                   *
%*******************************************************************************
%Este valor de cantidad de puntos de gauss se asume como el valor por defecto, es decir si no se
%lo especifica en los SETS, se adopta este valor (por eso ya no se la exige como variable
%obligatoria).
npg = f_ValVar('NGAUS',false,c_ValVarInic,'Variables iniciales');

%*******************************************************************************
%* NUMERO DE SET DISTINTOS                                                     *
%*******************************************************************************
nSet = f_ValVar('NSETS',true,c_ValVarInic,'Variables iniciales');

%*******************************************************************************
%* ADDITIONAL GEOMETRICAL PARAMETERS *
%*******************************************************************************
if exist('c_ValVarInic2','var')&&~isempty(c_ValVarInic2)
   nRefSmoothing  = f_ValVar('N_REF_SMOOTHING',false,c_ValVarInic2,'Variables addicionales');
   omega_micro    = f_ValVar('OMEGA_MICRO',false,c_ValVarInic2,'Variables addicionales');
   n_selec_mode   = f_ValVar('N_SELECT_MODE',false,c_ValVarInic2,'Variables addicionales');
   angSmoothingImp = f_ValVar('ANGSMOOTH',false,c_ValVarInic2,'Variables addicionales');
end


seccion = f_ProxString(fid);
f_VerifNom(seccion,'GENERAL_DATA',...
   'Datos generales: No est� definido el inicio con "GENERAL_DATA".')
seccion = f_ProxString(fid);
f_VerifNom(seccion,'GEOMETRY',...
   'Datos de geometr�a: No est� definido el inicio con "GEOMETRY".')

%*******************************************************************************
%* LISTA DE CONECTIVIDADES                                                     *
%*******************************************************************************
%Se debe ingresar en el archivo de datos, para las conectividades, la misma cantidad de columnas
%en todos los elementos. Estas son iguales a npe+2, y npe debe corresponder al elemento con mayor
%cantidad de nodos. En el archivo de datos, en los elementos con menos nodos se debe llenar las
%�ltimas columnas sobrantes con ceros (es decir, los nodos v�lidos deben estar primero, ya se
%descartan las �ltimas columnas sobrantes).
format = ['%f %f',repmat(' %f',1,npe)];
conec = textscan(fid,format,nElem,'CollectOutput',1,'CommentStyle','$');
conec = conec{1};
%En la lista de elementos se guarda la numeraci�n de nodos (no tiene que ser correlativa y
%completa), pero tiene que ser �nica.
m_NumElem = conec(:,1);
if length(m_NumElem)~=length(unique(m_NumElem))
   error(['Lectura de datos: Lectura de conectividades: La numeraci�n de los elementos ',...
      'debe ser �nica.'])
end
m_SetElem = conec(:,2)';
conec = conec(:,3:npe+2);

%*******************************************************************************
%* LISTA DE COORDENADAS Y NODOS                                                *
%*******************************************************************************
%Ver si no cambiar y en que 2D no sea necesario agregar la columna con ceros.
format = '%f %f %f %f';
xx = textscan(fid,format,nnod,'CollectOutput',1,'CommentStyle','$');
xx = xx{1};
in = xx(:,1);
xx = xx(:,2:4);
%inn = in;
if length(in)~=length(unique(in))
   error(['Lectura de datos: Lectura de conectividades: La numeraci�n de los nodos ',...
      'debe ser �nica.'])
end
%Se cambia las conectividades seg�n la numeraci�n interna del programa de los nodos (esta es
%seg�n el orden de las coordenadas ingresadas en el archivo de datos, y es completa), y no la
%adoptada en el archivo de datos y guardada en in (que puede ser una lista no ordenada y no
%completa, pero s� �nica).
%Se coloca conec>0 ya que se mantiene la convenci�n de que se usa malla distintas y uno de ellos
%tiene menos nodos, se rellena con ceros.
%Esta forma de hacer la operaci�n ocupa mucha memoria para el caso de mallas grandes, y la velocidad con
%que se resuelve la operaci�n vectorial se pierde si tiene que acceder al disco.
%[conec(conec~=0),~] = find(bsxfun(@eq,in,conec(conec~=0)'));
%La siguiente es menos exigente en memoria, pero m�s cpu intensivo.
%Como en mallas grandes esta operaci�n de cambio de numeraci�n de la conectividad puede llevar un tiempo se
%hace una verificaci�n para ver si es necesaria realizarla (en el caso que se ingrese un numeraci�n de nodo
%desordenada y/o discontinua).
%Esta operaci�n en una malla de 193101 nodos ordenada redujo los tiempos de 228.088387 seg. a 0.000741 seg
%en el MatLab 2015a. Esto en el caso que la numeraci�n de los nodos sea igual al �ndice de su posici�n.
ticIDCambNumRed = tic;
esComplNumNod = any(in~=(1:nnod)');
if esComplNumNod
   ticIDCambNum = tic;
   for iNod = 1:nnod
      conec(conec==in(iNod)) = iNod;
   end
   fprintf('Tiempo de cambio de numeraci�n de nodos de la conectividad: %f\n',toc(ticIDCambNum));
end
fprintf('Tiempo total del cambio de numeraci�n reducido: %f\n',toc(ticIDCambNumRed));

%*******************************************************************************
%* LISTA DE NODOS DE LA FRONTERA (Conectividades)                              *
%*******************************************************************************
%Se asume que los elementos de frontera tiene dos nodos (esto se tendr�a que leer del archivo de
%datos).
nNodElemFront = 2;
format = ['%f',repmat(' %f',1,nNodElemFront)];
m_ConecFront = textscan(fid,format,nElemFront,'CollectOutput',1,'CommentStyle','$');
%Se elimina la numeraci�n de los elementos de frontera (esta deber�a mantenerse si se imprimiera)
m_ConecFront = m_ConecFront{1}(:,2:nNodElemFront+1);
%Se cambia la enumercaci�n de los nodos a la utilizada internamente dentro del programa.
%Ver aclaraci�n para conec.
%[m_ConecFront(:),~] = find(bsxfun(@eq,in,m_ConecFront(:)'));
if esComplNumNod
for iNod = 1:nnod
   m_ConecFront(m_ConecFront==in(iNod)) = iNod;
end
end
%
seccion = f_ProxString(fid);
f_VerifNom(seccion,'END_GEOMETRY','Geometr�a de la Malla: Debe terminar con "END_GEOMETRY".')

%*******************************************************************************
%* PROPIEDADES DEL SET                                                         *
%*******************************************************************************
seccion = f_ProxString(fid);
exist_CrackPath=false;
f_VerifNom(seccion,'SETS','Propiedades de materiales: Se debe definir el inicio de "SETS".')
e_DatSet(nSet,1) = struct('e_DatElem',[],'e_DatMat',[],'conec',[],'m_IndElemSet',[],...
   'm_NumElem',[],'nElem',[],'sihvare',[],'siavare',[],'sitvare',[],...
   'm_BT',[],'m_DetJT',[],'m_VolElem',[]);
for iSet = 1:nSet
   
   %Sirve para que cada Set puedan
   %tener las estructuras e_DatElem y e_DatMat, definidas con distintos fields.
   clear e_DatElem e_DatMat
   %
   seccion = f_ProxString(fid);
   f_VerifNom(seccion,'SET','Propiedades Materiales: Cada material debe empezar con "SET".')
   set = textscan(fid,'%f',1,'Delimiter',' =','MultipleDelimsAsOne',1,'CommentStyle','$');
   set = set{1};
   %Verifica que los n�meros asignados al set son menor a la cantidad de set (esta forma permite
   %programar de forma m�s f�cil). Se utiliza para indexar set, en lugar de iSet, por si la lista
   %de sets en forma desordenada.
   if isempty(set)||set<=0||set>nSet
      %ver si hacer que se pueda usar n�mero cualquiera para identificar el material
      error(['Lectura de datos: Propiedades Materiales: Se debe indicar el n�mero del "SET" ',...
         'y el n�mero de set debe ser mayor que ser cero y menor o igual a la cantidad de sets.'])
   end
   
   %Lectura de los datos del elemento del set
   seccion = f_ProxString(fid);
   f_VerifNom(seccion,'ELEMENT_DATA',...
      'Propiedades Materiales: Falta los datos de elemento, "ELEMENT_DATA".')
   c_ValVarElem = f_ExtrVar(fid);
   
   %Tipo y modelo constitutivo del elemento
   eltype = f_ValVar('TYPE',true,c_ValVarElem,'Propiedades Materiales');
   conshyp = f_ValVar('MODEL',true,c_ValVarElem,'Propiedades Materiales');
   
   %Puntos de gauss del SET
   e_DatElem.npg = f_ValVar('NGAUS',false,c_ValVarElem,'Propiedades Materiales');
   %Se verifica que se indic� la cantidad de puntos de gauss en SET, y si no se usa los indicados
   %por defecto.
   if isnan(e_DatElem.npg)
      e_DatElem.npg = npg;
   end
   if isnan(e_DatElem.npg)||e_DatElem.npg<=0
      fclose(fid);
      error(['Lectura de datos: Lectura de los SETs: El n�mero de PG debe ser indicado en ',...
         'DIMENSIONS o en el SET, y debe ser mayor que cero'])
   end
   
   %Determinaci�n de cantidad de elementos correspondiente a este set
   %m_IndElemSet = m_SetElem==set;
   m_IndElemSet = find(m_SetElem==set);
   %e_DatSet(set).nElem = sum(m_IndElemSet);
   nElemSet = length(m_IndElemSet);
   m_NumElemSet = m_NumElem(m_IndElemSet)';
   e_DatSet(set).nElem = nElemSet;
   
   %Variables exigidas seg�n el tipo de elemento
   switch eltype
      case 2    % Triangulo de tres nodos
         %Espesor de elemento.
         e_DatElem.thickness = f_ValVar('THICKNESS',true,c_ValVarElem,'Propiedades Materiales');
         %e_DatElem.ntens = 4;
         nVarAuxElem = 0;
         nVarHistElem = 0;
         e_DatElem.npe = 3;
         [e_DatElem.xg,e_DatElem.wg] = set_var(eltype,e_DatElem.npg,fid);
         e_DatElem.wg = e_DatElem.wg*e_DatElem.thickness;
      case {4,8}  % Cuadrangulos de 4 nodos en desplazamientos
         %Espesor de elemento.
         e_DatElem.thickness = f_ValVar('THICKNESS',true,c_ValVarElem,'Propiedades Materiales');
         %e_DatElem.ntens = 4;
         nVarAuxElem = 0;
         nVarHistElem = 0;
         e_DatElem.npe = 4;
         [e_DatElem.xg,e_DatElem.wg] = set_var(eltype,e_DatElem.npg,fid);
         %Lo que se ingresa son los puntos de gauss por lado, y se necesita los totales.
         e_DatElem.npg = e_DatElem.npg*e_DatElem.npg;
         e_DatElem.wg = e_DatElem.wg*e_DatElem.thickness;
         %Esto sirve para cuando se usa con elementos micros de un problema multiescala donde el elemento
         %macro son elementos Q1-SD (23). Se debe indicar la el �ngulo de smoothing. Soluci�n es usar siempre
         %los elementos 31 � 32, pero de tipo el�stico, ya que lo tiene implementado.
         %p_IntDissip       = 1;
         %p_IntEnergy       = p_IntDissip + 1;
         %nVarHistElem = p_IntEnergy + e_DatElem.npg -1;
         %pointersVHE.i_IntDissip  = p_IntDissip    : p_IntEnergy-1 ;
         %pointersVHE.i_pIntEnergy  = p_IntEnergy    : nVarHistElem ;
         %e_DatElem.pointersVHE   = pointersVHE;

         %p_nSmoothing    = 1   ;
         %nVarAuxElem     = p_nSmoothing + 8-1 ; 
         %pointersVAE.p_nSmoothing   =  p_nSmoothing   : nVarAuxElem       ;
         %e_DatElem.pointersVAE      = pointersVAE;

      case 5     % Elemento de barra de 2 nodos en desplazamientos (2D)
         %�rea de la secci�n transversal del elemento.
         e_DatElem.area = f_ValVar('AREA',true,c_ValVarElem,'Propiedades Materiales');
         %e_DatElem.ntens = 1;
         nVarAuxElem = 0;
         nVarHistElem = 0;
         e_DatElem.npe = 2;
         [e_DatElem.xg,e_DatElem.wg] = set_var(eltype,e_DatElem.npg,fid);
         e_DatElem.wg = e_DatElem.wg*e_DatElem.area;
      case 7     % Elemento hexaedrico de 8 nodos (3D)
         nVarAuxElem = 0;
         nVarHistElem = 0;
         %No tiene una propiedad geom�trica adicional.
         %e_DatElem.ntens = 6;
         e_DatElem.npe = 8;
         [e_DatElem.xg,e_DatElem.wg] = set_var(eltype,e_DatElem.npg,fid);
         e_DatElem.npg = e_DatElem.npg*e_DatElem.npg*e_DatElem.npg;
      case 10    % Tri�ngulo de 3 nodos con discontinuidades fuertes (SDA)
         %Ancho de la zona de localizaci�n
         %(no se exige ingresar en el modelo multiescala 51 porque no se usa en este caso)
         e_DatElem.anchoLoc = f_ValVar('ANLOC',false,c_ValVarElem,'Propiedades Materiales');
         if conshyp~=51&&isnan(e_DatElem.anchoLoc)
            error(['Lectura de datos: Propiedades del EF: Se debe ingresar el ancho de ',...
               'localizaci�n ANLOC.'])
         end
         %Espesor de elemento.
         e_DatElem.thickness = f_ValVar('THICKNESS',true,c_ValVarElem,'Propiedades Materiales');
         %Si se quiere que se resuelva el equilibrio de tracciones con la metodolog�a sim�trica
         %(1) o no sim�trica (0, por defecto).
         %e_DatElem.simetrico = f_ValDefecto(f_ValVar('SIM',false,c_ValVarElem,...
         %   'Propiedades Materiales'),0,set,'SIM');
         %Se indica que este elemento hace que la matriz de rigidez global que no sea sim�trica.
         esKTSim = 0;
         %Numero de variables auxiliares elementales (componentes)
         %Se esta pasando una variable que indica la condicion de bifurcacion (1 componente), la normal a
         %la discontinuidad n en notacion de voigt (8 componentes), el tensor gradiente de phi de los nodos
         %solitarios en notacion de voigt (8 componentes), el paso que se produjo la bifurcacion (1
         %componente), se guarda la traccion (como vector), aunque no se utiliza en el calculo directamente
         %(se utiliza delta de sigma), otro valor (1 componente) indicando que no se permite que ese
         %elemento active la discontinuidad fuerte y que despu�s que bifurc� sea el�stico, y los dos
         %�ngulos obtenidos de la bifurcaci�n (2 componentes), colocando primero el de la normal a la
         %fisura.
         nVarAuxElem = 23;
         e_DatElem.npe = 3;
         %Est� programado para que el elemento tenga dos puntos de gauss en la misma posici�n,
         %por lo que se ignora cualquier valor impuesto en el archivo de datos.
         warning(['Lectura de datos: Propiedades del EF: Set %d: Se ignora los PG indicados, ',...
            'ya que este elemento est� realizado para 2 PG en la misma posici�n.'],set)
         [xgUnico,wgUnico] = set_var(eltype,1,fid);
         e_DatElem.xg = [xgUnico;xgUnico];
         e_DatElem.wg = [wgUnico;wgUnico];
         e_DatElem.npg = 2;
         e_DatElem.wg = e_DatElem.wg*e_DatElem.thickness;
         % Se utiliza para almacenar el salto y tracciones.
         nVarHistElem = 2*ndime;
      case 20       %Cuadr�ngulo de 4 nodos mixed con strain injection
         e_DatElem.npe = 4;
         %Extracci�n de variables
         m_ValVar = f_ValVar({'THICKNESS','ESTNOBIF','ESTBIF'},true,c_ValVarElem,'Propiedades Materiales');
         %Espesor de elemento.
         e_DatElem.thickness = m_ValVar(1);
         %Factores de estabilizaci�n para la zona que no tiene activado la bifurcaci�n y la que no tiene
         %activado la bifurcaci�n
         e_DatElem.estabNoBif = m_ValVar(2);
         e_DatElem.estabBif = m_ValVar(3);
         %Est� programado para que el elemento tenga cinco puntos de gauss, 4 ubicados seg�n la cuadratura
         %de gauss y otro en el centro por lo que se ignora cualquier valor impuesto en el archivo de
         %datos.
         warning(['Lectura de datos: Propiedades del EF: Set %d: Se ignora los PG indicados, ',...
            'ya que este elemento est� realizado para 5 PG.'],set)
         %Lo que se ingresa son los puntos de gauss por lado, y se necesita los totales.
         [xgEst,wgEst] = set_var(eltype,2,fid);
         xgCen = set_var(eltype,1,fid);
         e_DatElem.xg = [xgEst;xgCen];
         %Se aplica peso nulo al punto de gauss 5, ya que se utiliza para llevar la historia de tensiones,
         %deformaciones, variables internas del modelo constitutivo en ese punto.
         e_DatElem.wg = [wgEst;0];
         e_DatElem.npg = 5;
         e_DatElem.wg = e_DatElem.wg*e_DatElem.thickness;
         %N�mero de variables auxiliares del elemento.
         %1er Componente: Para almacenar la condici�n de bifurcaci�n.
         nVarAuxElem = 1;
         %N�mero de variables hist�ricas del elemento.
         %Se utiliza para almacenar las tensiones estabilizadas (un tensor de tensi�n por punto de gauss).
         nVarHistElem = ntens*e_DatElem.npg;
      case {21,22,23}       %Cuadr�ngulo de 4 nodos mixed con strain injection & SD MANUEL
         %
         esKTSim = 0;         
         exist_CrackPath = true;
         %
         e_DatElem.npe = 4;
         %Extracci�n de variables
         m_ExigProp = true(5,1);
         if eltype==23
            %En el elemento 23 actualmente no se utiliza LE.
            m_ExigProp(3) = false;
         end
         %La energ�a de fractura de referencia para la determinaci�n del instante de inyecci�n (GFVREF) no se
         %exige.
         m_ExigProp(5) = false;
         %
         m_ValVar = f_ValVar({'THICKNESS','KINF','LE','FACTINY','GFVREF'},m_ExigProp,c_ValVarElem,...
            'Propiedades Materiales');
         %Espesor de elemento.
         e_DatElem.thickness = m_ValVar(1);
         %
         e_DatElem.kinf = m_ValVar(2);
         %En el elemento 23 actualmente no se utiliza LE.
         if eltype==21||eltype==22
            e_DatElem.le = m_ValVar(3);
         end
         %Factor de inyecci�n
         e_DatElem.facIny = m_ValVar(4);
         %Energ�a de fractura de referencia
         e_DatElem.gfvRef = f_ValDefecto(m_ValVar(5),NaN,set,'GFVREF');
         %Est� programado para que el elemento tenga seis puntos de gauss, 4 ubicados seg�n la cuadratura
         %de gauss y otro dos en el centro por lo que se ignora cualquier valor impuesto en el archivo de
         %datos.
         warning(['Lectura de datos: Propiedades del EF: Set %d: Se ignora los PG indicados, ',...
            'ya que este elemento est� realizado para 6 PG.'],set)
         %Lo que se ingresa son los puntos de gauss por lado, y se necesita los totales.
         [xgEst,wgEst] = set_var(eltype,2,fid);
         xgCen = set_var(eltype,1,fid);
         e_DatElem.xg = [xgEst;xgCen];
         %Se aplica peso nulo al punto de gauss 5 y 6, ya que se utiliza para llevar la historia de tensiones,
         %deformaciones, variables internas del modelo constitutivo en ese punto.
         e_DatElem.wg = [wgEst;0;0];
         e_DatElem.npg = 6;
         e_DatElem.wg = e_DatElem.wg*e_DatElem.thickness;
         %N�mero de variables auxiliares del elemento.
         %Indice inicial de los variables auxiliares.
         %1er Componente: Para almacenar la condici�n de bifurcaci�n.
         p_condBif       = 1;
         p_elem_BifType  = p_condBif      + 1  ;
         p_leq_elem      = p_elem_BifType + 1 ;
         p_kSD           = p_leq_elem     + 1 ;
         p_phi_grad      = p_kSD          + 1 ;
         %Se almacena los dos �ngulos de bifurcaci�n obtenidos.
         p_n_tens        = p_phi_grad     + 8   ;
         %p_normal_bif    = p_n_tens       + 8   ;
         pAngBif = p_n_tens+8;
         pCeBif  = pAngBif+2;
         %p_fii           = p_normal_bif   + 4   ;
         %p_fii = pAngBif+2;
         p_fii = pCeBif+ntens^2;
         p_injFactor     = p_fii          + 4   ;
         p_ref_vector    = p_injFactor    + 1   ;
         p_nSmoothing    = p_ref_vector   + 2   ;         
         %nVarAuxElem     = p_nSmoothing   + 7   ;
         %Se almacena la tracci�n del PGS.
         p_Traccion = p_nSmoothing+8;
         nVarAuxElem = p_Traccion+1;         
         
         pointersVAE.p_condBif      =  p_condBif      : p_elem_BifType-1  ;
         pointersVAE.p_elem_Biftype =  p_elem_BifType : p_leq_elem-1      ;
         pointersVAE.p_leq_elem     =  p_leq_elem     : p_kSD-1           ;
         pointersVAE.p_kSD          =  p_kSD          : p_phi_grad-1      ;
         pointersVAE.p_phi_grad     =  p_phi_grad     : p_n_tens-1        ;
         %pointersVAE.p_n_tens       =  p_n_tens       : p_normal_bif-1    ;
         pointersVAE.p_n_tens       =  p_n_tens       : pAngBif-1;
         %pointersVAE.p_normal_bif   =  p_normal_bif   : p_fii-1           ;
         %pointersVAE.pAngBif        =  pAngBif        : p_fii-1;
         pointersVAE.pAngBif        =  pAngBif        : pCeBif-1;
         pointersVAE.pCeBif         =  pCeBif         : p_fii-1;
         pointersVAE.p_fii          =  p_fii          : p_injFactor-1     ;
        % pointersVAE.p_fii          =  p_fii          : p_injFactor-1     ;
         pointersVAE.p_injFactor    =  p_injFactor    : p_ref_vector-1    ;
         pointersVAE.p_ref_vector   =  p_ref_vector   : p_nSmoothing-1    ;
         %pointersVAE.p_nSmoothing = p_nSmoothing:nVarAuxElem;
         pointersVAE.p_nSmoothing = p_nSmoothing:p_Traccion-1;
         pointersVAE.p_Traccion = p_Traccion:nVarAuxElem;
         
         %N�mero de variables hist�ricas del elemento.
         %Se utiliza para almacenar las tensiones estabilizadas (un tensor de tensi�n por punto de gauss).
         p_indSTmacro    = 1   ;
         p_indActSTmacro = p_indSTmacro + 1   ;
         p_vectVHElem    = p_indActSTmacro + 1 ;
         %p_stressTilde   = p_vectVHElem  + 8 ;
         %Ver que el �ltimo valor guardado hist�rico no se est� usando (9), as� que se puede sacar. Hay que
         %cambiar las funci�n del elemento donde se almacena.
         p_stressTilde = p_vectVHElem+9;
         nVarHistElem = p_stressTilde+ntens*e_DatElem.npg-1;
         %
         pointersVHE.i_indST          = p_indSTmacro    : p_indActSTmacro -1 ;
         pointersVHE.p_indActSTmacro  = p_indActSTmacro : p_vectVHElem -1 ;
         pointersVHE.i_vectVHElem     = p_vectVHElem    : p_stressTilde -1 ;
         pointersVHE.i_stressTilde    = p_stressTilde   : nVarHistElem ;
         %
         %%
         e_DatElem.pointersVAE   = pointersVAE;
         e_DatElem.pointersVHE   = pointersVHE;
         
         % Variables asociadas al numero de cortes que hace el path en los elementos finitos
         e_DatElem.EF_sidesCutCPF = zeros(e_DatSet(set).nElem,1);
         e_DatElem.NumsidesCutCPF = zeros(e_DatSet(set).nElem,1);
         
      case 31  % Cuadrangulos de 4 nodos en desplazamientos con c�lculo de normal incorporada.
         %Espesor de elemento.
         e_DatElem.thickness = f_ValVar('THICKNESS',true,c_ValVarElem,'Propiedades Materiales');
         %e_DatElem.ntens = 4;
         nVarAuxElem = 0 ;
         e_DatElem.npe = 4;
         [e_DatElem.xg,e_DatElem.wg] = set_var(eltype,e_DatElem.npg,fid);
         %Lo que se ingresa son los puntos de gauss por lado, y se necesita los totales.
         e_DatElem.npg = e_DatElem.npg*e_DatElem.npg;
         e_DatElem.wg = e_DatElem.wg*e_DatElem.thickness;
         %
         p_IntDissip       = 1;
         p_IntEnergy       = p_IntDissip + 1;
         nVarHistElem = p_IntEnergy + e_DatElem.npg -1;
         pointersVHE.i_IntDissip  = p_IntDissip    : p_IntEnergy-1 ;
         pointersVHE.i_pIntEnergy  = p_IntEnergy    : nVarHistElem ;
         e_DatElem.pointersVHE   = pointersVHE;
      case 32  % Tri�ngulo de 3 nodos en desplazamientos con c�lculo de normal incorporada.
         %Espesor de elemento.
         e_DatElem.thickness = f_ValVar('THICKNESS',true,c_ValVarElem,'Propiedades Materiales');
         %e_DatElem.ntens = 4;
         %Se incorpora una variable auxiliar para almacenar el sentido actual de la normal al elemento.
         nVarAuxElem = 1;
         e_DatElem.npe = 3;
         [e_DatElem.xg,e_DatElem.wg] = set_var(eltype,e_DatElem.npg,fid);
         %Lo que se ingresa son los puntos de gauss por lado, y se necesita los totales.
         e_DatElem.npg = e_DatElem.npg;
         e_DatElem.wg = e_DatElem.wg*e_DatElem.thickness;
         %
         p_IntDissip       = 1;
         p_IntEnergy       = p_IntDissip + 1;
         nVarHistElem = p_IntEnergy + e_DatElem.npg -1;
         pointersVHE.i_IntDissip  = p_IntDissip    : p_IntEnergy-1 ;
         pointersVHE.i_pIntEnergy  = p_IntEnergy    : nVarHistElem ;
         e_DatElem.pointersVHE   = pointersVHE;
      %
      %%LARGE DEFORMATIONS     
      case 108     %Cuadr�ngulo bilineal FBar
         %Se est� incorporado en el tipo de problema, pero por si acaso se define de nuevo.
         %Con este elemento la matriz resulta no sim�trica.
         esKTSim = 0;
         %Espesor de elemento.
         e_DatElem.thickness = f_ValVar('THICKNESS',true,c_ValVarElem,'Propiedades Materiales');
         %
         nVarAuxElem = 0;
         nVarHistElem = 0;
         %
         e_DatElem.npe = 4;
         [e_DatElem.xg,e_DatElem.wg] = set_var(eltype,e_DatElem.npg,fid);
         %Lo que se ingresa son los puntos de gauss por lado, y se necesita los totales.
         e_DatElem.npg = e_DatElem.npg*e_DatElem.npg;
         e_DatElem.wg = e_DatElem.wg*e_DatElem.thickness;
         %Se almacena la matriz de deformaci�n del punto central para c�lculo del gradiente de deformaci�n del
         %punto central en las variables auxiliares.
         %Ver si no conviene crear un 5to PG.
         dofpeFBar = ndn*e_DatElem.npe;
         nVarAuxElem = nVarAuxElem+ntens*dofpeFBar;
      otherwise
         error('Propiedades del EF: Elemento finito no definido.')
   end
   %Cantidad de grados de libertad por elemento
   e_DatElem.dofpe = ndn*e_DatElem.npe;
   %Tipo de elemento
   e_DatElem.eltype = eltype;
   %Cantidad de variables auxiliares del elemento.
   e_DatElem.nVarAuxElem = nVarAuxElem;
   %Cantidad de variables hist�ricas del elemento.
   e_DatElem.nVarHistElem = nVarHistElem;
   
   %Lectura de los datos de los materiales del set
   seccion = f_ProxString(fid);
   f_VerifNom(seccion,'MATERIAL_DATA',...
      'Propiedades Materiales: Falta los datos del material, "MATERIAL_DATA".')
   c_ValVarMat = f_ExtrVar(fid);
   
   % Variables exigidas seg�n el tipo de material
   %(tambi�n se podr�a configurar variables que si no son ingresadas tiene un valor por defecto).
   switch conshyp
      case 1      % ELASTICIDAD LINEAL
         nhvart = 0;      % numero de variables historicas tensoriales por punto de Gauss
         nhvare = 0;      % numero de variables historicas escalares por punto de Gauss
         nauxvar = 0;     % numero de variables auxiliares escalares por punto de Gauss
         %
         m_ValVar = f_ValVar({'YOUNG';'POISS'},true,c_ValVarMat,'Propiedades Materiales');
         e_DatMat = struct('young',m_ValVar(1),'poiss',m_ValVar(2));
      case 2      % ELASTO - PLASTICIDAD J2: HARDENING-SOFTENING ISOTROPICO
         nhvart = 2;
         nhvare = 2;
         nauxvar = 1;
         %
         %La funci�n de m�dulo de endurecimiento isotr�pico es (p�g. 91 del Computational Elasticity - Simo):
         %K(alpha) = FTULT+TIT*HBA*alfa+(KIN-KCE)*(1-exp(-DEL*alfa))
         m_ValVar = f_ValVar({'YOUNG';'POISS';'FTULT';'HBA';'TIT';'KIN';'KCE';'DEL'},...
            true,c_ValVarMat,'Propiedades Materiales');
         e_DatMat = struct('young',m_ValVar(1),'poiss',m_ValVar(2),'ftult',m_ValVar(3),...
            'hba',m_ValVar(4),'tit',m_ValVar(5),'kin',m_ValVar(6),'kce',m_ValVar(7),...
            'del',m_ValVar(8));
      case 3      % VISCO - PLASTICIDAD J2: HARDENING ISOTROPICO
         error('Para este modelo falta ver cu�les son las variables que hay que definir')
         %nhvart = 1;
         %nhvare = 1;
         %nauxvar = 1;
      case 4     % DA�O IS�TROPO
         nhvart = 0;
         nhvare = 2;
         nauxvar = 2;
         %
         m_ValVar = f_ValVar({'YOUNG';'POISS';'FTULT';'HBA';'TIT'},[true,true,true,true,false],...
            c_ValVarMat,'Propiedades Materiales');
         e_DatMat = struct('young',m_ValVar(1),'poiss',m_ValVar(2),'ftult',m_ValVar(3),...
            'hba',m_ValVar(4),'tit',[]);
         %TIT toma por defecto el valor 0 (ya que se inicializa Eprop con ceros)
         e_DatMat.tit = f_ValDefecto(m_ValVar(5),0,set,'TIT');
         %Se agrega el valores iniciales de las variables internas de evoluci�n para no tener calcularlas
         %cada vez que se entra a la funci�n de da�o (ver que conviene en el caso del cluster).
         e_DatMat.r_0 = e_DatMat.ftult/sqrt(e_DatMat.young);
      case 5     % DA�O IS�TROPO SOLO TRACCI�N
         nhvart = 0;
         nhvare = 2;
         nauxvar = 2;
         %
         m_ValVar = f_ValVar({'YOUNG';'POISS';'FTULT';'HBA';'TIT'},[true,true,true,true,false],...
            c_ValVarMat,'Propiedades Materiales');
         e_DatMat = struct('young',m_ValVar(1),'poiss',m_ValVar(2),'ftult',m_ValVar(3),...
            'hba',m_ValVar(4),'tit',[]);
         %TIT toma por defecto el valor 0 (ya que se inicializa Eprop con ceros)
         e_DatMat.tit = f_ValDefecto(m_ValVar(5),0,set,'TIT');
         %Se agrega el valores iniciales de las variables internas de evoluci�n para no tener calcularlas
         %cada vez que se entra a la funci�n de da�o (ver que conviene en el caso del cluster).
         e_DatMat.r_0 = e_DatMat.ftult/sqrt(e_DatMat.young);
         %Se indica que este elemento hace que la matriz de rigidez global que no sea sim�trica.
         esKTSim = 0;
      case 6     % VISCO - DA�O
         error('Para este modelo falta ver cu�les son las variables que hay que definir')
         %nhvart = [];
         %nhvare = [];
         %nauxvar = [];
      case 7     % DA�O - PLASTICIDAD DP
         error('Para este modelo falta ver cu�les son las variables que hay que definir')
         %nhvart = 1;
         %nhvare = 5;
         %nauxvar = 6;
      case 8     % DA�O - FUERZAS CENTRADAS
         nhvart = 0;
         nhvare = 13;
         nauxvar = 10;
         %
         m_ValVar = f_ValVar({'YOUNG';'POISS';'FTULT';'GFV'},true,c_ValVarMat,...
            'Propiedades Materiales');
         e_DatMat = struct('young',m_ValVar(1),'poiss',m_ValVar(2),'ftult',m_ValVar(3),...
            'gfv',m_ValVar(4));
         e_DatMat.tit = f_ValDefecto(m_ValVar(5),0,set,'TIT');
         %Se indica que este elemento hace que la matriz de rigidez global que no sea sim�trica.
         esKTSim = 0;
      case 9     % ELASTO - PLASTICIDAD J2: HARDENING Y SOFTENING ISOTROPICO BILINEAL
         error('Para este modelo falta ver cu�les son las variables que hay que definir')
         %nhvart = 2;
         %nhvare = 2;
         %nauxvar = 1;
         %nauxvar = nauxvar+ntens^2; % nauxvar+e_DatElem.ntens^2;
      case 10    % DA�O IS�TROPO REGULARIZADO
         nhvart = 0;
         nhvare = 2;
         nauxvar = 2;
         %
         m_ValVar = f_ValVar({'YOUNG';'POISS';'FTULT';'GFV';'TIT'},[true,true,true,true,false],...
            c_ValVarMat,'Propiedades Materiales');
         e_DatMat = struct('young',m_ValVar(1),'poiss',m_ValVar(2),'ftult',m_ValVar(3),...
            'gfv',m_ValVar(4),'tit',[]);
         %TIT toma por defecto el valor 0 (ya que se inicializa Eprop con ceros)
         e_DatMat.tit = f_ValDefecto(m_ValVar(5),0,set,'TIT');
         %Se lee la distancia regularizaci�n impuesta, si es una n�mero se impone la misma a
         %todos los elementos, si es string, se interpreta como una archivo de texto que indica
         %por fila el n�mero de elemento y la distancia de regularizaci�n. Si no se indica nada,
         %significa que se calcula una regularizaci�n seg�n una geometr�a del elemento est�ndar
         %(un tri�ngulo is�sceles recto o un cuadrado, seg�n el elemento, por ejemplo).
         %Se lleva una matriz de elementos del set.
         hReg = f_ValVar('HREG',false,c_ValVarMat,'Propiedades Materiales');
         hReg = f_ValDefecto(hReg(1),NaN,set,'HREG');
         nomArchHReg = f_ValStrVar('ARCHHREG',false,c_ValVarMat,'Propiedades Materiales');
         nomArchHReg = f_ValDefecto(nomArchHReg{1},NaN,set,'ARCHHREG');
         e_DatMat.m_hReg = f_LongReg(hReg,nomArchHReg,path_file,m_NumElemSet,nElemSet);
         %Se agrega el valores iniciales de las variables internas de evoluci�n para no tener calcularlas
         %cada vez que se entra a la funci�n de da�o (ver que conviene en el caso del cluster).
         e_DatMat.r_0 = e_DatMat.ftult/sqrt(e_DatMat.young);
         %Se indica que este elemento hace que la matriz de rigidez global que no sea sim�trica.
         esKTSim = 0;
      case 11   % DA�O IS�TROPO S�LO TRACCI�N REGULARIZADO
         %N�mero de variables hist�ricas tensoriales.
         nhvart = 0;
         %N�mero de variables hist�ricas escalares.
         % 4 variables indice de inyeccion y 3 del modelo implex y 3 de nueva version implex tau_n y tau_n-1
         %nhvare  = 7;
         %nhvare  = 8;
         nhvare  = 9;
         %Numero de variables auxiliares.
         %Para guardar el factor de carga (fLoad).
         nauxvar = 1;
         %
         m_ValVar = f_ValVar({'YOUNG';'POISS';'FTULT';'GFV';'TIT';'HREG';'IMPLEX'},...
            [true,true,true,true,false,false,false],c_ValVarMat,'Propiedades Materiales');
         e_DatMat = struct('young',m_ValVar(1),'poiss',m_ValVar(2),'ftult',m_ValVar(3),...
            'gfv',m_ValVar(4),'tit',[]);
         %TIT toma por defecto el valor 0 (ya que se inicializa Eprop con ceros)
         e_DatMat.tit = f_ValDefecto(m_ValVar(5),0,set,'TIT');
         %Se lee la distancia regularizaci�n impuesta, si es un n�mero se impone la misma a todos los
         %elementos, si es string, se interpreta como una archivo de texto que indica por fila el n�mero de
         %elemento y la distancia de regularizaci�n. Si no se indica nada, significa que se calcula una
         %regularizaci�n seg�n una geometr�a del elemento est�ndar (un tri�ngulo is�sceles recto o un
         %cuadrado, seg�n el elemento, por ejemplo). Se lleva una matriz de elementos del set.
         hReg = f_ValDefecto(m_ValVar(6),NaN,set,'HREG');
         nomArchHReg = f_ValStrVar('ARCHHREG',false,c_ValVarMat,'Propiedades Materiales');
         nomArchHReg = f_ValDefecto(nomArchHReg{1},NaN,set,'ARCHHREG');
         e_DatMat.m_hReg = f_LongReg(hReg,nomArchHReg,path_file,m_NumElemSet,nElemSet);
         %Se define si se resuelve en forma Impl�cita-Expl�cita (Implex).
         %Por defecto es impl�cito (0 � false)
         e_DatMat.esImplex = f_ValDefecto(m_ValVar(7),false,set,'IMPLEX');            
         %Se agrega el valores iniciales de las variables internas de evoluci�n para no tener calcularlas
         %cada vez que se entra a la funci�n de da�o (ver que conviene en el caso del cluster).
         e_DatMat.r_0 = e_DatMat.ftult/sqrt(e_DatMat.young);
         %En el caso del exponencial para una correcta regularizaci�n cuando cambie el espesor de
         %regularizaci�n.
         %Notar que se asume que estas dos variables hist�ricas adicionales est�n al final de la matriz, es
         %decir en la posici�n end y end-1.
         if e_DatMat.tit==1  %Exponencial
            nhvare = nhvare+2;
         end            
         %Se indica que este elemento hace que la matriz de rigidez global que no sea sim�trica.
         esKTSim = 0;
      case {12,13}  % tipo:12=DANIO SOLO TRACCION #### tipo:13=DANIO rankine
         %N�mero de variables hist�ricas tensoriales.
         nhvart = 0;
         %N�mero de variables hist�ricas escalares.
         %Se agrega una para almacenar los dos �ngulos de bifurcaci�n y si bifurc� o no.
         nhvare = 5;
         %N�mero de variables auxiliares.
         %nauxvar = 2;
         %Se agrega una para indicarle si el punto de gauss tiene que compartarse el�sticamente.
         nauxvar = 3;
         %
         m_ValVar = f_ValVar({'YOUNG';'POISS';'FTULT';'GFV';'TIT';'HREG';'IMPLEX'},...
            [true,true,true,true,false,false,false],c_ValVarMat,'Propiedades Materiales');
         e_DatMat = struct('young',m_ValVar(1),'poiss',m_ValVar(2),'ftult',m_ValVar(3),...
            'gfv',m_ValVar(4),'tit',[]);
         %TIT toma por defecto el valor 0 (ya que se inicializa Eprop con ceros)
         e_DatMat.tit = f_ValDefecto(m_ValVar(5),0,set,'TIT');
         %Se lee la distancia regularizaci�n impuesta, si es un n�mero se impone la misma a todos los
         %elementos, si es string, se interpreta como una archivo de texto que indica por fila el n�mero de
         %elemento y la distancia de regularizaci�n. Si no se indica nada, significa que se calcula una
         %regularizaci�n seg�n una geometr�a del elemento est�ndar (un tri�ngulo is�sceles recto o un
         %cuadrado, seg�n el elemento, por ejemplo). Se lleva una matriz de elementos del set.
         hReg = f_ValDefecto(m_ValVar(6),NaN,set,'HREG');
         nomArchHReg = f_ValStrVar('ARCHHREG',false,c_ValVarMat,'Propiedades Materiales');
         nomArchHReg = f_ValDefecto(nomArchHReg{1},NaN,set,'ARCHHREG');
         e_DatMat.m_hReg = f_LongReg(hReg,nomArchHReg,path_file,m_NumElemSet,nElemSet);
         %Se define si se resuelve en forma Impl�cita-Expl�cita (Implex).
         %Por defecto es impl�cito (0 � false)
         e_DatMat.esImplex = f_ValDefecto(m_ValVar(7),false,set,'IMPLEX');
         if e_DatMat.esImplex
            %En el implex se necesita almacenar el incremento de r como variable hist�rica.
            nhvare = nhvare+1;
            %Si es implex se almacena tambi�n el factor del incremento de r.
            nauxvar = nauxvar+1;
         end
         %Se agrega el valores iniciales de las variables internas de evoluci�n para no tener calcularlas
         %cada vez que se entra a la funci�n de da�o (ver que conviene en el caso del cluster).
         e_DatMat.r_0 = e_DatMat.ftult/sqrt(e_DatMat.young);
         %Se indica que este elemento hace que la matriz de rigidez global que no sea sim�trica.
         esKTSim = 0;
      case {50,55}    %50: MULTIESCALA MODELO CL�SICO
                      %55: MULTIESCALA MODELO CL�SICO CON AN�LISIS DE BIFURCACI�N
         nhvart = 0;
         nhvare = 1;
         nauxvar = 0;
         c_ValVar = f_ValStrVar({'ARCHDAT','IMPRCU'},[true,false],c_ValVarMat,...
            'Propiedades Materiales');
         % Archivo de elementos y PG que tiene se imprimir el postproceso de la celda unitaria
         imprResCU = f_ValDefecto(c_ValVar{2},NaN,set,'IMPRCU');
         m_ElemPGImpr = f_LectElemPGImpr(imprResCU,path_file,m_NumElemSet);
         %Propiedades del modelo multiescala cl�sico
         e_DatMat = f_VarME(c_ValVar{1},path_file,m_ElemPGImpr);
      case 51    % MULTIESCALA MODELO COHESIVO OBJETIVO
         nhvart = 0;
         nhvare = 1;
         nauxvar = 1;
         c_ValVar = f_ValStrVar({'ARCHDAT','IMPRCU'},[true,false],c_ValVarMat,'Propiedades Materiales');
         % Archivo de elementos y PG que tiene se imprimir el postproceso de la celda unitaria
         imprResCU = f_ValDefecto(c_ValVar{2},NaN,set,'IMPRCU');
         m_ElemPGImpr = f_LectElemPGImpr(imprResCU,path_file,m_NumElemSet);
         % IMPLEX: Esquema de integraci�n temporal impl�cita-expl�cita
         %Si en la microestructura se consider� modelos constitutivos integrada en forma IMPLEX se debe
         %activar esta opci�n para realizar el an�lisis de bifurcaci�n con el tensor constitutivo impl�cito
         %y realizar el reemplazo de las tensiones olds implex con las olds impl�citas. Por ello para que
         %la estructura macro se resuelva con Implex, en forma adicional, se debe definir los modelos
         %constitutivos que as� lo sean (se podr�a realizar una mezcla entre modelos impl�citos e implex en
         %la microestructura, ya que en ese caso si un set micro no tiene activado el implex se toma el
         %t�rmino ct cl�sico).
         %Por defecto es impl�cito (0 � false)
         esImplex = f_ValDefecto(f_ValVar('IMPLEX',false,c_ValVarMat,'Propiedades Materiales'),...
            false,set,'IMPLEX');
         %Propiedades del modelo multiescala cohesivo
         e_DatMat = f_VarMECohesivo(c_ValVar{1},path_file,m_ElemPGImpr,esImplex);
      case 52    % MULTIESCALA, CELDA UNITARIA EL�STICA
         nhvart = 0;
         nhvare = 0;
         nauxvar = 0;
         %A este modelo constitutivo se lo trata como uno el�stico, pero en lugar de usar el tensor
         %constitutivo el�stico est�ndar se utiliza el determinado de la resoluci�n de la celda unitaria.
         c_ValVar = f_ValStrVar('ARCHDAT',true,c_ValVarMat,'Propiedades Materiales');
         %Se lee los datos de la celda unitaria (se utiliza la funci�n del modelo constitutivo 50).
         %Notar que esta estructura se borra, y se deja solo que tenga como field el tensor
         %tangente constitutivo ce. Se podr�a guardar todos los datos de la microescala, pero no ocupar
         %espacio innecesario se los elimina.
         %Ver si no integrarlo con el modelo el�stico agregando la propiedad ARCHDAT.
         e_DatMat = f_VarME(c_ValVar{1},path_file,[]);
      case 53    % MULTIESCALA, Barcelona
         nhvart = 0;
         nhvare = 1;
         nauxvar = 0;
         %A este modelo constitutivo se lo trata como uno el�stico, pero en lugar de usar el tensor
         %constitutivo el�stico est�ndar se utiliza el determinado de la resoluci�n de la celda unitaria.
         c_ValVar = f_ValStrVar({'ARCHDAT','IMPRCU','IMPLEX'},[true,false,false],c_ValVarMat,...
            'Propiedades Materiales');
         % Archivo de elementos y PG que tiene se imprimir el postproceso de la celda unitaria
         imprResCU = f_ValDefecto(c_ValVar{2},NaN,set,'IMPRCU');
         m_ElemPGImpr = f_LectElemPGImpr(imprResCU,path_file,m_NumElemSet);
         %Por defecto es impl�cito (0 � false)
         esImplex = f_ValDefecto(c_ValVar{3},false,set,'IMPLEX');
         %Se lee los datos de la celda unitaria (se utiliza la funci�n del modelo constitutivo 50).
         %Notar que esta estructura se borra, y se deja solo que tenga como field el tensor
         %tangente constitutivo ce. Se podr�a guardar todos los datos de la microescala, pero no ocupar
         %espacio innecesario se los elimina.
         %Ver si no integrarlo con el modelo el�stico agregando la propiedad ARCHDAT.
         e_DatMat = f_VarMEBcna(c_ValVar{1},path_file,m_ElemPGImpr,esImplex);
      case 54    % MULTIESCALA, SANTA FE
         nhvart = 0;
         nhvare = 1;
         nauxvar = 1;
         c_ValVar = f_ValStrVar({'ARCHDAT','IMPRCU'},[true,false],c_ValVarMat,'Propiedades Materiales');
         imprResCU = f_ValDefecto(c_ValVar{2},NaN,set,'IMPRCU');
         m_ElemPGImpr = f_LectElemPGImpr(imprResCU,path_file,m_NumElemSet);
         esImplex = f_ValDefecto(f_ValVar('IMPLEX',false,c_ValVarMat,'Propiedades Materiales'),...
            false,set,'IMPLEX');
         e_DatMat = f_VarMESFe(c_ValVar{1},path_file,m_ElemPGImpr,esImplex);
      %
      %%LARGE DEFORMATIONS   
      case 100 %Elastic Material neo-Hookean
         nhvart = 0;
         nhvare = 0;
         nauxvar = 0;
         %
         m_ValVar = f_ValVar({'COMPR';'SHEAR'},true,c_ValVarMat,'Propiedades Materiales');
         e_DatMat = struct('modCompr',m_ValVar(1),'modCorte',m_ValVar(2));
      case 110  %J2 Plasticity   
         nhvart = 0;
         nhvare = 2+4;
         nauxvar = 1;
         %
         %La funci�n de m�dulo de endurecimiento isotr�pico es (p�g. 91 del Computational Elasticity - Simo):
         %K(alpha) = YIELD+TIT*HBA*alfa+(KIN-KCE)*(1-exp(-DEL*alfa))
         m_ValVar = f_ValVar({'YOUNG';'POISS';'YIELD';'HBA';'TIT';'KIN';'KCE';'DEL'},...
            true,c_ValVarMat,'Propiedades Materiales');
         e_DatMat = struct('young',m_ValVar(1),'poiss',m_ValVar(2),'yield',m_ValVar(3),...
            'hba',m_ValVar(4),'tit',m_ValVar(5),'kin',m_ValVar(6),'kce',m_ValVar(7),...
            'del',m_ValVar(8));
      otherwise
         error('Lectura de datos: Propiedades Materiales: Modelo constitutivo no definido.')
   end
   e_DatMat.conshyp = conshyp;
   %Se guarda estas matrices en la estructura de materiales por depende del modelo constitutivo, y
   %para no tener que pasar toda la estructura e_DatSet dentro del parfor, que implicar�a
   %posiblemente copiar conec, etc (todas matrices grandes, verificar esto).
   e_DatMat.sihvarpg = (nhvart*ntens+nhvare);
   e_DatMat.siavarpg = nauxvar;
   %Para no tener que verificar que el campo existe en las funciones de los elementos finitos (isfield) en
   %los modelos constitutivos que no tienen definido el implex, se fija que todos los modelos constitutivos
   %tengan definido esa propiedad. Como en los problemas multiescala, las funciones elementales se llaman
   %much�simas veces, se trata evitar operaciones simples pero que multiplicadas por la cantidad de veces
   %que se llama incremente el tiempo.
   if ~isfield(e_DatMat,'esImplex')
      %Por defecto se define que sea impl�cita.
      e_DatMat.esImplex = false;
   end
   
   %% Variables de la estructura principal (e_DatSet)
   e_DatSet(set).sihvare = e_DatMat.sihvarpg*e_DatElem.npg;
   %e_DatSet(set).sihvarpg = (nhvart*ntens+nhvare);
   %
   e_DatSet(set).siavare = nauxvar*e_DatElem.npg;
   %e_DatSet(set).siavarpg = nauxvar;
   %
   e_DatSet(set).sitvare = ntens*e_DatElem.npg;   %e_DatElem.ntens*e_DatElem.npg;
   %
   % Matrices de conectividades de cada set
   e_DatSet(set).conec = conec(m_IndElemSet,1:e_DatElem.npe);
   % Matrices de grados de libertad de cada elemento
   %Se realiza un prec�lculo para acelerar el c�lculo, aunque puede conllevar mayor transferencia al
   %paralelizar.
   e_DatSet(set).m_DofElem = reshape(f_DofElem(reshape(e_DatSet(set).conec',1,[]),ndn),e_DatElem.dofpe,[]);
   % Posici�n de los elementos en la matriz de conectividades global
   %(la que se ingresa en el archivo de datos) y es la numeraci�n interna dentro del programa.
   %Esta solo sirve en el caso que se utilice matrices donde una dimensi�n es la cantidad de
   %elementos (esta puede ser �til para simplificar ciertas operaciones pero ver si no queda muy
   %confuso y no conviene utilizar para todo lo que se refiere a matrices de elementos,
   %estructuras o celdas).
   %Esta matriz indica que posici�n tiene los elementos correspondiente al set dentro de la
   %matriz global de elementos (ordenados seg�n como fueron indicados en el archivo de datos).
   %Esta matriz global es la que se utiliza para indentificar a los elementos dentro del
   %programa.
   e_DatSet(set).m_IndElemSet = m_IndElemSet;
   % Matrices de numeraci�n de elementos indicada en el archivo de datos (solo para postproceso)
   %Esta es como se denomina los elementos en el archivo de datos (primera columna), esta puede
   %ser una lista no ordenada y no completa (puede faltar n�meros), pero tiene ser �nica (es
   %decir cada elemento tiene que ser identificado en forma �nica en la lista del archivo de
   %datos.
   %Notar que no se guarda m_SetElem, ya que todos los elementos de cada e_DatSet(set) tiene el
   %mismo n�mero de set.
   e_DatSet(set).m_NumElem = m_NumElemSet;
   % Almacenamiento de las estructuras principales
   e_DatSet(set).e_DatMat = e_DatMat;
   e_DatSet(set).e_DatElem = e_DatElem;
   
   seccion = f_ProxString(fid);
   f_VerifNom(seccion,'END_SET','Propiedades Materiales: Cada set debe terminar con "END_SET".')
   
end
seccion = f_ProxString(fid);
f_VerifNom(seccion,'END_SETS',...
   'Propiedades de materiales: Se debe definir el final de los sets, "END_SETS".')

%% FIN LECTURA DE DATOS GENERALES DEL MODELO
seccion = f_ProxString(fid);
f_VerifNom(seccion,'END_GENERAL_DATA',...
   'Datos generales: No est� definido el fin con "END_GENERAL_DATA".')

%*******************************************************************************
%* TIPO DE ELEMENTO FINITO (SOLO PUEDE DEFINIRSE UN TIPO DE ELEMENTO FINITO)   *
%*******************************************************************************
ndoft = nnod*ndn;

%*******************************************************************************
%* NUMERO DE PASOS DE TIEMPO                                                   *
%*******************************************************************************
seccion = f_ProxString(fid);
f_VerifNom(seccion,'INTERVAL_DATA','Datos del intervalo de tiempo: Se debe ingresar "INTERVAL_DATA".')
c_ValVarTiempo = f_ExtrVar(fid);
np = f_ValVar('NSTEP',true,c_ValVarTiempo,'Variables temporales');

%*******************************************************************************
%* INCREMENTO DE TIEMPO                                                        *
%*******************************************************************************
Dtime = f_ValVar('DTIME',true,c_ValVarTiempo,'Variables temporales');

%*******************************************************************************
%* FUNCION TEMPORAL PARA CONDICIONES DE BORDE                                  *
%*******************************************************************************
seccion = f_ProxString(fid);
f_VerifNom(seccion,'FUNCTION',['Funcion de evolucion temporal: Se debe empezar la funcion ',...
   'con "FUNCTION".'])
%Definir que variables si necesario ingresar variables para definir la funci�n.
c_ValVarFunTiempo = f_ExtrVar(fid);
%
formatFunTiempo = '%f %f';
funbc = textscan(fid,formatFunTiempo,'CollectOutput',1,'CommentStyle','$');
funbc = funbc{1};
seccion = f_ProxString(fid);
f_VerifNom(seccion,'END_FUNCTION',...
   'Funci�n de evoluci�n temporal: Se debe terminar la funci�n con "END_FUNCTION".')

%*******************************************************************************
%* FZAS ESTATICAS                                                              *
%*******************************************************************************
seccion = f_ProxString(fid);
f_VerifNom(seccion,'LOAD','Cargas aplicadas: Se debe definir las cargas, "LOAD".')
nomCarga = f_ProxString(fid);
tipoCarga = f_ProxString(fid);
f = zeros(ndoft,1);
switch tipoCarga
   case 'POINT_LOAD'
      c_NomVarFza = {'NODE','FX','FY','FZ'};
      if ndime==2
         m_ValExigCarga = [true,true,true,false];
      elseif ndime==3
         m_ValExigCarga = [true,true,true,true];
      else
         error('Lectura de datos: Cargas aplicadas: Cargas puntuales: Dimensi�n del espacio no definido.')
      end
      c_ValVarCarga = f_ExtrVar(fid);
      while ~strcmpi(c_ValVarCarga{1},'END_POINT_LOAD')
         m_ValFza = f_ValVar(c_NomVarFza,m_ValExigCarga,c_ValVarCarga,...
            'Cargas aplicadas: Cargas puntuales');
         %Con find(in==m_ValFza(1)) se convierte la numeraci�n de nodos a la interna.
         m_gdlNod = f_DofElem(find(in==m_ValFza(1)),ndn);
         %Se realiza una suma de fuerzas, considerando que puede poner varias l�neas con fuerzas en el
         %mismo nodo, o cuando se defina distintos tipos de carga.
         f(m_gdlNod) = f(m_gdlNod)+m_ValFza(2:ndn+1);
         c_ValVarCarga = f_ExtrVar(fid);
      end
   otherwise
      error('Lectura de datos: Cargas aplicadas: Tipo de carga no definida')
end
seccion = f_ProxString(fid);
f_VerifNom(seccion,'END_LOAD',...
   'Cargas aplicadas: Se debe definir el final de la secci�n de ingreso de cargas, "END_LOAD".')

%*******************************************************************************
%* CONDICIONES DE CONTORNO                                                     *
%*******************************************************************************
seccion = f_ProxString(fid);
f_VerifNom(seccion,'BOUNDARY','Condiciones de borde: Se debe definir las condiciones de borde, "BOUNDARY".')
% tipoBound = {{''}};
% while(~strcmpi(tipoBound{1}{1},'BOUNDARY'))
%    tipoBound = textscan(fid,'%s',1,'Delimiter',' :','MultipleDelimsAsOne',1,'CommentStyle','$');
% end
%tipoBound = tipoBound{2}{1};
%Se extrae las opciones de las condiciones de borde
c_ValVarCondBord = f_ExtrVar(fid);
if numel(c_ValVarCondBord)==1
   %En el caso que se ingrese solo el nombre de las condiciones de borde. En ese caso no se sigue la l�gica de
   %nombre de variable=valorVariable.
   s_nomCondBorde = c_ValVarCondBord{1};
   matCBFull = false;
else
   s_nomCondBorde = f_ValStrVar('NOM',true,c_ValVarCondBord,'Variables de condiciones de borde');
   matCBFull = f_ValDefecto(f_ValVar('MATCBFULL',false,c_ValVarCondBord,'Variables de condiciones de borde'),...
      false,[],'MATCBFULL');
end
%Capaz se podr�a usar %c para leer por separado las direcciones que est�n restringidas
%Se generaliza la lectura del valor de la restricci�n, se agrega una columna de opciones, que es optativo
%agregar en el archivo de datos. Esto permite agregar las coordenadas de la periodicidad en el caso de
%problemas de un grado de libertad.
%En el caso de poner la columna de opciones hay que poner un resultado en todas las filas.
%Por ahora no se pone una opci�n por grado de libertad, igual, en el caso de que querer poner opciones
%distintas para cada grado de libertad deber�a funcionar poner filas separadas para cada uno de ellos.
%formatBou = ['%f %f',repmat(' %f',1,ndn),repmat(' %s',1,ndn),'\n'];
%Cambiar que lea todo string entre par�ntesis, as� se puede agregar espacios dentro de los par�ntesis.
formatBou = ['%f %f',repmat(' %f',1,ndn),' %s\n'];
m_CondBord = textscan(fid,formatBou,'CollectOutput',true,'CommentStyle','$');
%Se ordena los datos ingresados seg�n la lista de nodos, para asegurar que los �ndices de los
%grados de libertad libre y fijos est�n bien calculados, y ordenados seg�n la matriz de rigidez
%global.
c_OpCondBord = m_CondBord{2};
[m_CondBord,m_IndOrd] = sortrows(m_CondBord{1});
c_OpCondBord = c_OpCondBord(m_IndOrd,:);
%Se cambia la numeraci�n de los nodos a la interna.
[m_CondBord(:,1),~] = find(bsxfun(@eq,in,m_CondBord(:,1)'));
%Lectura de comando de fin de boundary
seccion = f_ProxString(fid);
f_VerifNom(seccion,'END_BOUNDARY',['Condiciones de Borde: El pr�ximo string a las ',...
   'condiciones de borde tiene que ser "END_BOUNDARY".'])

%*******************************************************************************
%* PARAMETROS DEL ALGORITMO NO-LINEAL                                          *
%*******************************************************************************
seccion = f_ProxString(fid);
f_VerifNom(seccion,'STRATEGY','Estrategia: Se definir iniciar la misma con "STRATEGY".')
nomEstrat = f_ProxString(fid);

% BUSQUEDA LINEAL
seccion = f_ProxString(fid);
f_VerifNom(seccion,'LINE_SEARCH',['Estrategia: B�squeda lineal: Se debe indicar si activar la ',...
   'b�squeda lineal en "LINE_SEARCH".'])
lsearch = f_ProxString(fid);
if strcmp(lsearch,'OFF')
   lsearch = false;
else
   lsearch = true;
end

% SECCI�N DE DATOS DE CONVERGENCIA
seccion = f_ProxString(fid);
f_VerifNom(seccion,'CONVERGENCE',['Estrategia: Datos de convergencia: Se debe indicar los datos de ',...
   'convergencia con "CONVERGENCE".'])
c_ValVarConv = f_ExtrVar(fid);

% TOLERANCIA ESQUEMA DE NEWTON GLOBAL
tolnr = f_ValVar('TOL_NEWTON',true,c_ValVarConv,['Estrategia: Datos de convergencia: Tolerancia del ',...
   'esquema de newton global']);

% N�MERO DE ITERACIONES M�XIMAS
%El valor por defecto es 50 iteraciones.
iterMax = f_ValDefecto(f_ValVar('ITERMAX',false,c_ValVarConv,['Estrategia: Datos de convergencia: ',...
   'Iteraciones m�ximas del esquema de newton global']),50,[],'ITERMAX');


% TOLERANCIA NEWTON CONSTITUTIVO
toldg = f_ValVar('TOL_CONST_MODEL',false,c_ValVarConv,['Estrategia: Datos de convergencia: ',...
   'Tolerancia del esquema de newton constitutivo']);

% PARAMETROS BUSQUEDA LINEAL
m_ParamBL = f_ValVar({'ETA_MIN','ETA_MAX','ETA_PAR','TOLER'},lsearch,c_ValVarConv,...
   'Estrategia: B�squeda lineal: Par�metros');
eta_min = m_ParamBL(1);
eta_max = m_ParamBL(2);
eta_par = m_ParamBL(3);
tolls = m_ParamBL(4);

%*******************************************************************************
%* ESTRATEGIA DE CONTROL                                                       *
%*******************************************************************************
stratControl = f_ValVar('STRATEGY',true,c_ValVarConv,['Estrategia: Datos de convergencia: ',...
   'Estrategia de control']);
switch stratControl
   case 1
      lambda = 1.0;
      ds = nan;
   case 4
      lambda = 1.0e-08;
      ds = f_ValVar('DS',true,c_ValVarConv,['Estrategia: Datos de convergencia: ',...
         'Estrategia de control xx']);
      node1 = f_ValVar('NODE1',true,c_ValVarConv,['Estrategia: Datos de convergencia: ',...
         'NODE1']);
      dof1 = f_ValVar('DOF1',true,c_ValVarConv,['Estrategia: Datos de convergencia: ',...
         'DOF1']);
      node2 = f_ValVar('NODE2',false,c_ValVarConv,['Estrategia: Datos de convergencia: ',...
         'NODE2']);
      dof2 = f_ValVar('DOF2',false,c_ValVarConv,['Estrategia: Datos de convergencia: ',...
         'DOF2']);
   otherwise
      error('Lectura de datos: Datos de convergencia: Estrategia de control: Estrategia no definida.')
end

%*******************************************************************************
%* POST-PROCESO                                                                *
%*******************************************************************************
seccion = f_ProxString(fid);
f_VerifNom(seccion,'POST_PROCESS','PostProceso: Falta indicar los datos de postproceso,"POST_PROCESS".')
c_VarPostProc = f_ExtrVar(fid);
m_VarPostProc = f_ValVar({'STEP','MULTFILE'},[false,false],c_VarPostProc,'PostProceso');
postpro_impre_step = f_ValDefecto(m_VarPostProc(1),-1,[],'STEP');
isPostResMultFile = f_ValDefecto(m_VarPostProc(2),false,[],'MULTFILE');

%*******************************************************************************
%* RESULTADOS EN NODOS                                                         *
%*******************************************************************************
%Lectura de la pr�xima secci�n
seccion = f_ProxString(fid);
if strcmpi(seccion,'PLOT')
   %Seg�n el orden indicado en esta celda es como se codifica el tipo dato.
   c_NomDat = {'T','Dx','Dy','Dz','Fx','Fy','Fz','Bx','By','Bz','Tx','Ty','Tz',...
      'Exx','Eyy','Ezz','Exy','Eyx','Txx','Tyy','Tzz','Txy','Tyx','Efxx','Efyy','Efzz','Efxy',...
      'Da'};
   c_ValVarNGraf = f_ExtrVar(fid);
   nGraf = f_ValVar('NGRAF',true,c_ValVarNGraf,'Gr�ficos X-Y');
   if nGraf>0
      m_DatGrafXY = nan(nGraf,8);
      for iGraf =  1:nGraf
         c_ValVarGraf = f_ExtrVar(fid);
         %Se recupera el nodo, elemento y punto de gauss de los datos que se quiere graficar en los
         %ejes X e Y.
         m_DatGrafPos = f_ValVar({'Nx','Ny','Ex','Ey','PGx','PGy'},....
            [false,false,false,false,false,false],c_ValVarGraf,'Gr�ficos X-Y');
         c_DatGrafTipo = f_ValStrVar({'X','Y'},[true,true],c_ValVarGraf,'Gr�ficos X-Y');
         %Se recupera y se codifica que variable se quiere imprimir.
         m_DatGrafXY(iGraf,7) = find(strcmpi(c_DatGrafTipo{1},c_NomDat));
         m_DatGrafXY(iGraf,8) = find(strcmpi(c_DatGrafTipo{2},c_NomDat));
         % Se transforma seg�n la numeraci�n interna
         %Se utiliza la numeraci�n global interna (es decir la que se obtiene por el orden de los
         %elementos en la matriz de conectividad completa)
         %En el caso que no se ingrese alg�n nodo, elemento o punto de gauss, se deja la matriz
         %con los valores NaN. Esto puede hacer que alguna gr�fica le falta alguna informaci�n de
         %que dato se quiere graficar, por luego seg�n el modelo constitutivo y tipo de elemento
         %se verifica esto (notar que por ejemplo que si se grafica el tiempo, no es necesasrio
         %ingresar ning�n dato).
         for iEje = 1:2
            %X e Y de los datos nodales
            if ~isnan(m_DatGrafPos(iEje))
               m_DatGrafPos(iEje) = find(m_DatGrafPos(iEje)==in);
            end
            %X e Y de los datos elementales
            if ~isnan(m_DatGrafPos(iEje+2))
               m_DatGrafPos(iEje+2) = find(m_DatGrafPos(iEje+2)==m_NumElem);
            end
         end
         m_DatGrafXY(iGraf,1:6) = m_DatGrafPos;
      end
   else
      m_DatGrafXY = [];
   end
   %
   %Lectura de la pr�xima secci�n
   seccion = f_ProxString(fid);
else
   %En el caso que no se indique la secci�n PLOT
   m_DatGrafXY = [];
end

%*******************************************************************************
%* IMPRESION DE RESULTADOS                                                     *
%*******************************************************************************
if strcmpi(seccion,'PRINT')
   IRES = textscan(fid,'IRES=%f','CollectOutput',1,'CommentStyle','$');
   %Lectura de la pr�xima secci�n
   seccion = f_ProxString(fid);
else
   %Valores por defecto
   IRES = 1;
end

%*******************************************************************************
%* FINALIZACI�N SECCI�N ESTRATEGIA                                             *
%*******************************************************************************
f_VerifNom(seccion,'END_STRATEGY','Estrategias num�ricas: El string final tiene que ser "END_STRATEGY".')

%*******************************************************************************
%* SE CIERRA ARCHIVO                                                           *
%*******************************************************************************
fclose(fid);

%*******************************************************************************
%* OPERADORES TENSORIALES                                                      *
%*******************************************************************************
[SOIT,SSOIT,FOSIT,FOSPT,FODPT,FOAT1,FOAT2,SONT] = rmtens(struhyp,ntens);

%Se obtiene el nombre completo
[~,filename] = fileparts(file);
fileCompleto = fullfile(path_file,filename);
%Se genera el nombre del archivo que guarda los resultados para el GID
%filename_res = [fileCompleto,'.flavia.res'];

%% Variable usada posteriormente
max_time = Dtime*np;

%Variables generales o globales. La intenci�n es que estas variables no se modifiquen a partir de
%esta parte del c�digo, o por lo menos si se modifica en una funci�n a esta modificaci�n se la
%utilice en las funciones que llama esta (funciones internas), es decir que no sea necesario
%devolver como resultado e_VG.
%Para pasar en forma temporal el Set actual  a las funciones que est�n dentro del parfor del ensamblaje se
%agrega la variable iSet (esto solo sirve para debug).
K_GlobElast= sparse([]);
e_VG = struct('struhyp',struhyp,'conshyp',conshyp,...
   'ndime',ndime,'nnod',nnod,'ndoft',ndoft,'nElem',nElem,'nSet',nSet,'ndn',ndn,'ntens',ntens,...
   'np',np,'Dtime',Dtime,'max_time',max_time,...
   'tolnr',tolnr,'toldg',toldg,'tolls',tolls,'iterMax',iterMax,...
   'SOIT',SOIT,'SSOIT',SSOIT,'FOSIT',FOSIT,'FOSPT',FOSPT,'FODPT',FODPT,'FOAT1',FOAT1,...
   'FOAT2',FOAT2,'SONT',SONT,...
   'fileCompleto',fileCompleto,...
   'CONTROL_STRAT',stratControl,'ds',ds,'lambda',lambda,...
   'm_CondBord',m_CondBord,'c_OpCondBord',{c_OpCondBord},'m_ConecFront',m_ConecFront,'matCBFull',matCBFull,...
   'istep',[],'iter',[],'iSet',[],'iElem',[],'iElemSet',[],'iElemNum',[],'iPG',[],...
   'esME',0,'elast',0,'esKTSim',esKTSim,...
   'postpro_impre_step',postpro_impre_step,'isPostResMultFile',isPostResMultFile,...
   'IRES',IRES,'m_DatGrafXY',m_DatGrafXY,...
   'exist_CrackPath',exist_CrackPath,'K_GlobElast',K_GlobElast);

if stratControl==4
   e_VG.NODE1= node1;
   e_VG.DOF1 = dof1;
   if  ~isnan(node2)
      e_VG.NODE2 = node2;
   end
   if  ~isnan(dof2)
      e_VG.DOF2  = dof2;
   end
end

%% Prec�lculo de las matrices geom�tricas de los elementos finitos
%Ver si no cambiar la funci�n y meterla en el loop de los set de arriba (el problema que no est�
%definida la e_VG a ese nivel).
e_DatSet = f_MatBT(xx,e_DatSet,e_VG);

%% Prec�lculo de las matrices el�sticas y de las longitudes de regularizaci�n del modelo de
%%da�o regularizado
%Es necesario ponerlo al final del read_dat para que est� definido el e_VG (ver si no cambiar la
%funci�n c_elas para que no sea necesario pasarle e_VG, y directamente las variables que necesita)
ticIDPreCalc = tic;
for iSet = 1:nSet
   
   e_DatMatSet = e_DatSet(iSet).e_DatMat;
   e_DatElemSet = e_DatSet(iSet).e_DatElem;
   eltype = e_DatElemSet.eltype;
   
   %calculo de datos para elemento tipo 31 - 32 (BANDAS DE LOCALIZACION)
   nElem = e_DatSet(iSet).nElem;
   conec  = e_DatSet(iSet).conec;
   ksb = zeros(nElem,1);
   switch eltype
      case 31  %cuadrilatero
         for iElem = 1:nElem
            xx_elem=xx(conec(iElem,:),:);
            % decide which size is the larger one !
            x21= (xx_elem(2,1:2)-xx_elem(1,1:2));
            x41= (xx_elem(4,1:2)-xx_elem(1,1:2));
            if norm(x21) > norm(x41)
               normal= [x21(2) ; -x21(1) ]/norm([x21(2) ; -x21(1) ]) ;
               ksb(iElem)= abs(x41*normal);
            else
               normal= [x41(2) ; -x41(1) ]/norm([x41(2) ; -x41(1) ]);
               ksb(iElem)= abs(x21*normal);
            end
         end
         e_DatSet(iSet).e_DatElem.ksb = ksb;         
      case 32   %tri�ngulo de 3 nodos
         normal_micro = zeros(2,nElem);
         le_micro     = zeros(nElem,1);
         %Lectura de las normales micro de un archivo de datos.
         %El archivo est� organizado de la siguiente forma: NroElem AngNorm LongElem
         fId = fopen([fileCompleto,'.normElem'],'rt');
         if fId~=-1
            m_AngNormMic = fscanf(fId,'%f %f %f\n',[3,Inf])';
            fclose(fId);
         else
            m_AngNormMic = zeros(0,3);
         end
         m_NumElemSet = e_DatSet(iSet).m_NumElem;
         %Se espera el n�mero de elemento y el �ngulo.
         for iElem = 1:nElem
            xx_elem=xx(conec(iElem,:),:);
            % decide which size is the larger one !
            x(1,:)= (xx_elem(2,1:2)-xx_elem(1,1:2));
            len(1) = norm(x(1,:));
            x(2,:)= (xx_elem(3,1:2)-xx_elem(1,1:2));
            len(2) = norm(x(2,:));
            x(3,:)= (xx_elem(3,1:2)-xx_elem(2,1:2));
            len(3) = norm(x(3,:));
            ksb(iElem)= min(len);
            ind=find(len==ksb(iElem)) ;
            if ind==1
               n = [-x(2,2) x(2,1)];
               m = [-x(3,2) x(3,1)];
               p = [ x(1,1) x(1,2)];
            elseif ind==2
               n = [-x(1,2) x(1,1)];
               m = [-x(3,2) x(3,1)];
               p = [ x(2,1) x(2,2) ];
            elseif ind==3
               n = [-x(2,2) x(2,1)];
               m = [-x(1,2) x(1,1)];
               p = [ x(3,1) x(3,2)];
            end
            %
            if n*m' < 0
               n=-n;
            end
            %
            m_IndSet = m_NumElemSet(iElem)==m_AngNormMic(:,1);
            if any(m_IndSet)
               ang = m_AngNormMic(m_IndSet,2);
               normal_micro(:,iElem) = [cos(ang),sin(ang)];
               %Se divide la longitud en dos porque despu�s se hace una suma en todos los elementos para
               %el c�lculo de la longitud total y cuando se integra sobre S, y siempre hay dos tri�ngulos
               %superpuestos.
               le_micro(iElem) = m_AngNormMic(m_IndSet,3)/2;
            else
               %En el caso de no encontrarse en el archivo de datos se calcula.
               normal_micro(:,iElem) =(n+m)'/norm(n+m);
               %Se guarda la mitad de la longitud porque despu�s se hace una suma en todos los elementos para
               %el c�lculo de la longitud total y cuando se integra sobre S, y siempre hay dos tri�ngulos
               %superpuestos.
               le_micro(iElem)= 0.25*(norm(n)+norm(m));
            end
            %
            ksb(iElem)= abs(p*normal_micro(:,iElem));
         end
         e_DatSet(iSet).e_DatElem.ksb = ksb;
         e_DatSet(iSet).e_DatElem.normal_micro = normal_micro;
         e_DatSet(iSet).e_DatElem.le_Elem = le_micro;
      case 23
         %Definici�n de la energ�a de fractura de referencia usada para definir cuando inyectar el SD.
         if isnan(e_DatSet(iSet).e_DatElem.gfvRef)
            switch conshyp
               case 11
                  e_DatSet(iSet).e_DatElem.gfvRef = e_DatSet(iSet).e_DatMat.gfv;
               case {53,54}
                  %Seleccionar en duro el set material en la microescala donde se extrae la energ�a de
                  %fractura de referencia.
                  set = 3;
                  e_DatSet(iSet).e_DatElem.gfvRef = e_DatSet(iSet).e_DatMat.e_DatSet(set).e_DatMat.gfv;
                  warning(['Lectura de datos: Inicializaci�n variables: Elemento finito Q1-SD: ',...
                     'Energ�a de fractura de referencia: Se utiliza el set material %d de la microescala ',...
                     'para definir la energ�a de fractura de referencia del elemento Q1-SD.'],set)
               otherwise
                  error(['Lectura de datos: Inicializaci�n variables: Elemento finito SD-Q1: ',...
                     'Energ�a de fractura de referencia: Modelo constitutivo no definido.'])
            end
         end
   end
   
   %Inicializaci�n de variables seg�n el modelo constitutivo
   thickness = e_DatElemSet.thickness;
   switch e_DatMatSet.conshyp
      case {1,2,3,4,5,6,7,8,9}
         e_DatSet(iSet).e_DatMat.ce = c_elas(e_DatMatSet.young,e_DatMatSet.poiss,e_VG);
      case {10,11,12,13}
         e_DatSet(iSet).e_DatMat.ce = c_elas(e_DatMatSet.young,e_DatMatSet.poiss,e_VG);
         %
         m_NroIndHRegCalc = find(e_DatMatSet.m_hReg<=0);
         nHRegCalc = length(m_NroIndHRegCalc);
         for iCalc = 1:nHRegCalc
            nroIndHReg = m_NroIndHRegCalc(iCalc);
            tipoCalc = e_DatMatSet.m_hReg(nroIndHReg);
            switch tipoCalc
               case -1
                  %Se calcula segun un elemento de forma estandar con el mismo area.
                  switch eltype
                     case {2,10,32} %Tri�ngulo
                        %Para obtener una longitud del elemento se considera como un prisma recto
                        %de altura igual a espesor con bases de tri�ngulos rect�ngulos.
                        %Notar que esto es solo para problemas 2D.
                        %hElemReg = sqrt(2*volElem/espesor);
                        %Es importante usar e_DatSet(iSet).e_DatMat en lugar de e_DatMatSet para
                        %que los cambios queden guardados, ya que sino solo se modifica una copia
                        %propia de hReg que tiene e_DatMatSet, mientras que la hReg de
                        %e_DatSet(iSet) queda sin modificar.
                        e_DatSet(iSet).e_DatMat.m_hReg(nroIndHReg) = sqrt(...
                           2*e_DatSet(iSet).m_VolElem(nroIndHReg)/thickness);
                     case {4,8,20,21,22,23,31} %Cuadr�ngulo
                        %Para obtener una longitud del elemento se considera como un prisma recto
                        %de altura igual espesor con bases de cuadrados.
                        %hElemReg = sqrt(volElem/espesor);
                        e_DatSet(iSet).e_DatMat.m_hReg(nroIndHReg) = sqrt(...
                           e_DatSet(iSet).m_VolElem(nroIndHReg)/thickness);
                     otherwise
                        error(['Lectura de datos: Inicializaci�n variables: Modelo constitutivo ',...
                           'de da�o regularizado: C�lculo de longitud de regularizaci�n: ',...
                           'Tipo de c�lculo -1: Elemento finito no definido.'])
                  end
               case -2
                  %Se calcula en funci�n de la geometr�a del elemento.
                  switch eltype
                     case {31,32} %Elementos Bandas: Cuadr�ngulo y Tri�ngulo
                        %C�mo estos elementos calculan un espesor de regularizaci�n a partir de su geometr�a
                        %se utiliza en el caso de no imponerse explic�tamente ning�n valor.
                        %Ver que para estos elementos y material de da�o 11 no lo termina utilizando
                        %finalmente ya que usa directamente la matriz e_DatSet(iSet).e_DatElem.ksb.
                        %Para las otras funciones de da�os s� lo utilizar�a.
                        e_DatSet(iSet).e_DatMat.m_hReg(nroIndHReg) = ...
                                 e_DatSet(iSet).e_DatElem.ksb(nroIndHReg);
                     otherwise
                        error(['Lectura de datos: Inicializaci�n variables: Modelo constitutivo ',...
                           'de da�o regularizado: C�lculo de longitud de regularizaci�n: ',...
                           'Tipo de c�lculo -2: Elemento finito no definido.'])
                  end
               otherwise
                  error(['Lectura de datos: Inicializaci�n variables: Modelo constitutivo ',...
                     'de da�o regularizado: C�lculo Longitud de regularizaci�n: Tipo de c�lculo no ',...
                     'definido.'])
            end
         end
      case {50,51,53,55}
         %Se pone estos modelos constitutivos para indicar que est� definido el modelo
         %constitutivo y no salte error.
      case 54
         %Evaluaci�n de la matriz global del subdominio elastico
%          clear e_DatMat
%          e_VGMicro = e_DatMatSet.e_VG;
%          e_DatSetMicro = e_DatMatSet.e_DatSet;
%          [K_GlobElast,BTCe_Glob] = f_K_GlobElast(e_DatSetMicro,e_VGMicro);
%          
%          e_DatSet(iSet).e_DatMat.e_VG.K_GlobElast = K_GlobElast;
%          e_DatSet(iSet).e_DatMat.e_VG.BTCe_Glob   = BTCe_Glob;
      case 52
         %Se borra la estructura que antes se usa como variable temporal para transformar la estructura en una
         %correspondiente al modelo elastico.
         clear e_DatMat
         e_VGMicro = e_DatMatSet.e_VG;
         e_DatSetMicro = e_DatMatSet.e_DatSet;
         nElemMicro = e_VGMicro.nElem;
         [c_CTMicro,KTMicro] = f_TensTangElastCU(e_DatSetMicro,e_VGMicro);
         [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro,e_DatMatSet.xx,...
            e_DatSetMicro,e_VGMicro.m_ConecFront);
         e_DatMat.ce = f_ModTangHomog(KTMicro,c_CTMicro,m_LinCondMicro,doflMicro,doffMicro,...
            e_DatSetMicro,e_DatMatSet.omegaMicro,true(nElemMicro,1),true(nElemMicro,1),...
            e_VGMicro);
         e_DatMat.conshyp  = e_DatSet(iSet).e_DatMat.conshyp;
         e_DatMat.sihvarpg = e_DatSet(iSet).e_DatMat.sihvarpg;
         e_DatMat.siavarpg = e_DatSet(iSet).e_DatMat.siavarpg;
         e_DatMat.esImplex = e_DatSet(iSet).e_DatMat.esImplex;
         %Ac� se sobreescribe la estructura de e_DatMat que contiene todos los datos de la microcelda, y se
         %guarda la que tiene �nicamente los campos ce, conshyp, sihvarpg y siavarpg.
         e_DatSet(iSet).e_DatMat = e_DatMat;
      %
      %LARGE DEFORMATIONS
      case {100,110}
      otherwise
         error('Lectura de datos: Inicializaci�n variables: Modelo constitutivo no definido.')
   end

end
fprintf('Tiempo total de operaciones de pre-c�lculo en el preproceso: %f\n',toc(ticIDPreCalc))

if e_VG.exist_CrackPath
   e_VG.smooth_alpha = zeros(nnod,1);
   e_VG.smooth_dalpha = zeros(nnod,1);
   e_VG.smooth_dalpha_old = zeros(nnod,1);
   % *************************
   % * MATRIZ DE MASA LUMPED *
   % * ***********************
   [e_VG,e_DatSet] = Mmg_SDA(xx,e_VG,e_DatSet);
   %
   if exist('nRefSmoothing','var')&&~isnan(nRefSmoothing)
      e_VG.nRefSmoothing = [cos(nRefSmoothing);sin(nRefSmoothing)];
   end
   if exist('n_selec_mode','var')&&~isnan(n_selec_mode)
      e_VG.n_selec_mode = n_selec_mode;
   end
   if exist('angSmoothingImp','var')&&~isnan(angSmoothingImp)
      e_VG.angSmoothingImp = angSmoothingImp;
   else
      nomArch = [fileCompleto,'.AngSmooth'];
      if exist(nomArch,'file')
         fId = fopen(nomArch,'rt');
         e_VG.m_AngSmoothImpEl = fscanf(fId,'%f %f\n',[2,inf])';
         %Se ordena los �ngulos seg�n la numeraci�n interna de los elementos (la numeraci�n ingresada en el
         %archivo es la nominal, la que se ingresa en el archivo de datos en la primera de columna de las
         %conectividades).
         %Se usa matriz con NaN verificar que se hayan ingresado los �ngulos para todos los elementos
         m_AngSmoothImpEl = NaN(e_VG.nElem,1);
         for iElLec = 1:size(e_VG.m_AngSmoothImpEl,1)
            m_AngSmoothImpEl(m_NumElem==e_VG.m_AngSmoothImpEl(iElLec,1)) = e_VG.m_AngSmoothImpEl(iElLec,2);
         end
         e_VG.m_AngSmoothImpEl = m_AngSmoothImpEl;
         %
         fclose(fId);
      end
   end
end

end

function string = f_ProxString(fId)

   %Encuentra la pr�xima secci�n (busca el pr�ximo string ignorando comentarios y los signo : y =).
   string = '';
   while isempty(string)
      string = textscan(fId,'%s',1,'Delimiter',' :=','MultipleDelimsAsOne',1,'CommentStyle','$');
      string = string{1}{1};
   end

end

function f_VerifNom(texto1,texto2,mjeError)

   %Se verifica que el texto1 sea igual al texto2, sino se muestra el mensaje de error mjeError.
   if ~strcmpi(texto1,texto2)
      error(['Lectura de datos: ',mjeError])
   end

end

function c_Var = f_ExtrVar(fId)

   %Se extrae todas las valores que est� sobre una l�nea de las variables indicadas en c_ListVar.
   %Se considera una l�nea hasta que encuentra un salto l�nea, donde el s�mbolo "\" indica que la
   %l�nea continua en la siguiente.
   c_finLinea = {'NOVACIO'};
   textLeido = '';
   %textscan(fId,'%[:]',1,'CommentStyle','$');
   while ~isempty(c_finLinea{1})
      %Cuidado que el textscan considera los comentarios al inicio de un field, as� si la l�nea
      %completa (usando %[^/\n], se lee hasta que encuentra el car�cter / o salto
      %de l�nea \n, pero sin leerlos) tiene un comentario al final o en el medio, este se lee.
      %Lo que se hace es ignorar todo despu�s del s�mbolo de comentario
      c_textLeido = textscan(fId,'%[^/$\n]',1,'CommentStyle','$');
      textLeido = [textLeido,' ',c_textLeido{1}{1}]; %#ok<AGROW>
      %Se lee �nicamente el car�cter / en el primer elemento de la celda (primer field), luego hasta
      %final de l�nea (2do field). Esto se hace para ver que si se debe leer la pr�xima l�nea o no.
      %Se ignora todo lo que viene despu�s de / para que cuando se lea la siguiente, se empiece de
      %una l�nea nueva.
      c_finLinea = textscan(fId,'%[/] %[^\n]',1);
   end

   %Se busca los valores de la variables
   %El texto siguiente al igual, =, puede ser un n�mero o un string (si tiene espacios tiene que
   %estar entre comillas).
   c_Var = textscan(textLeido,'%s %q','Delimiter',' =','MultipleDelimsAsOne',1,...
      'CommentStyle','$');
   %Se transforma en una celda pura
   c_Var = [c_Var{1,1},c_Var{1,2}];

end

function m_ValNumVar = f_ValVar(c_NomVar,m_VarExig,c_Var,nomSecVar)

%Devuelve los valores de las variables definidas en c_NomVar. Esta funci�n est� pensada para
%recuperar valores num�ricos.
%En m_VarExig se indica si la variable es exigida (con 1 o true), en la cual si no est� definida
%se tira un error.
%En el caso que si el dato no se puede transformar a un valor num�rico (se utiliza la funci�n str2double),
%sin importar si es una variable exigida, se emitir� un valor de error (se asume que cuando dentro del
%programa se llama a esta funci�n, quiere decir que se exige que esa variable debe devolver un valor
%num�rico).

   if ischar(c_NomVar)
      nVar = 1;
      c_NomVar = {c_NomVar};
   else
      nVar = length(c_NomVar);
   end

   if nargin==3
      nomSecVar = '';
   else
      nomSecVar = [nomSecVar,': '];
   end
   mjeError = ['Lectura de datos: ',nomSecVar];

   if length(m_VarExig)==1
      m_VarExig = repmat(m_VarExig,nVar,1);
   end

   %Se transforma los valores de las variables en valores num�ricos en los casos que se pueda (en el caso que
   %no sea queda NaN en la matriz).
   m_ValNumTodasVar = str2double(c_Var(:,2));   %Ver si no poner esta transformaci�n dentro del loop.
   m_ValNumVar = zeros(nVar,1);
   for iVar = 1:nVar
      m_IndVarBuscada = strcmpi(c_NomVar{iVar},c_Var(:,1));
      if any(m_IndVarBuscada)
         if sum(m_IndVarBuscada)==1
            valVar = m_ValNumTodasVar(m_IndVarBuscada);
            if ~isnan(valVar)
               m_ValNumVar(iVar) = valVar;
            else
               error([mjeError,'La variable "%s" no tiene un valor num�rico v�lido.'],c_NomVar{iVar})
            end
         else
            error([mjeError,'La variable "%s" est� definida m�s de una vez.'],c_NomVar{iVar})
         end
      else
         if m_VarExig(iVar)
            error([mjeError,'La variable "%s" no est� definida en los datos.'],c_NomVar{iVar})
         else
            m_ValNumVar(iVar) = NaN;
         end
      end
   end

end

function c_ValStrVar = f_ValStrVar(c_NomVar,m_VarExig,c_Var,nomSecVar)

%Devuelve los valores de las variables definidas en c_NomVar. Esta funci�n est� pensada para
%recuperar valores tipo string.
%En m_VarExig se indica si la variable es exigida (con 1 o true), en la cual si no est� definida
%se tira un error.

   if ischar(c_NomVar)
      nVar = 1;
      c_NomVar = {c_NomVar};
   else
      nVar = length(c_NomVar);
   end

   if nargin==3
      nomSecVar = '';
   else
      nomSecVar = [nomSecVar,': '];
   end
   mjeError = ['Lectura de datos: ',nomSecVar];

   if length(m_VarExig)==1
      m_VarExig = repmat(m_VarExig,nVar,1);
   end

   c_ValStrVar = cell(nVar,1);
   for iVar = 1:nVar
      m_IndVarBuscada = strcmpi(c_NomVar{iVar},c_Var(:,1));
      if any(m_IndVarBuscada)
         if sum(m_IndVarBuscada)==1
            c_ValStrVar{iVar} = c_Var{m_IndVarBuscada,2};
         else
            error([mjeError,'La variable "%s" est� definida m�s de una vez.'],c_NomVar{iVar})
         end
      else
         if m_VarExig(iVar)
            error([mjeError,'La variable "%s" no est� definida en los datos.'],c_NomVar{iVar})
         else
            c_ValStrVar{iVar} = NaN;
         end
      end
   end

end

function valor = f_ValDefecto(valor,valorDef,iSet,var)

   if isnan(valor)
      valor = valorDef;
      if isnumeric(valorDef)||islogical(valorDef)
         valorDef = num2str(valorDef);
      end
      %warning('off','backtrace')
      if isempty(iSet)
         warning('Se adopta el valor por defecto %s=%s.',var,valorDef) %#ok<WNTAG>
      else
         warning('SET %d: Se adopta el valor por defecto %s=%s.',iSet,var,valorDef) %#ok<WNTAG>
      end
      %warning('on','backtrace')
   end

end