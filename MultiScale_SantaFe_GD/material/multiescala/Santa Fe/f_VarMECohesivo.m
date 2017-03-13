function e_DatMatME = f_VarMECohesivo(nomArchME,dirDat,m_ElemPGImpr,esImplex)

   %Se asume que el archivo de datos de la celda unitaria está en el mismo directorio el archivo de
   %datos de la malla macro.
   %Los datos de tiempo, ignora los leídos de este archivo (micro) y usa los leídos a nivel macro.
   nomArchDat = [nomArchME,'.mfl'];
   [in,xx,m_SetElem,f,funbc,e_DatSet,e_VG] = read_data(nomArchDat,dirDat);
   
   %Para decirle al programa que se está dentro del cálculo multiescala
   e_VG.esME = 1;
   %Para impresión y debug se guarda la iteración macro (notar que este campo solo se define en el
   %e_VG micro)
   e_VG.iterMacro = [];
   
   % Cálculo de la medida de la celda unitaria
   %Cómo se hace operaciones sobre el volumen, se guarda el volumen de los elementos en una sola matriz.
   %Ver si esto duplica los resultados y mejor obtener la matriz m_VolElem cada vez que se utilice.
   %Se la ordena según la numeración (orden) de la matriz de conectividad global (como fueron 
   %ingresados los elementos).
   m_VolElem = zeros(1,e_VG.nElem);
   m_VolElem([e_DatSet.m_IndElemSet]) = [e_DatSet.m_VolElem];
   %En caso que se ingrese en el script el área micro, se utiliza directamente esa.
   if ~exist('Omega_micro','var')
      %Esta expresión es correcta solo si la celda unitaria no tiene agujeros.
      Omega_micro = sum(m_VolElem);
   end
   
   % Matriz de elementos vecinos
   m_ElemVec = f_ElemVecino(e_DatSet,e_VG);
   
   e_DatMatME = struct('in',in,'xx',xx,'m_SetElem',m_SetElem,'f',f,'funbc',funbc,'e_DatSet',e_DatSet,...
      'm_ElemPGImpr',m_ElemPGImpr,'m_VolElem',m_VolElem,'omegaMicro',Omega_micro,...
      'm_ElemVec',m_ElemVec,'tipoCBL',[],'pasosPosBif',[],'esImplex',esImplex,'e_VG',e_VG);
   
   % Script con variables 
   %(ver si correlo acá, que significa que hay que llevar varias variables o en cada paso tiempo,
   %donde habría que analizar cuánto tiempo tarda)
   %Ocurre un error al usar usar run, no sé si pasa lo mismo con eval, ya que matlab no se da cuenta
   %que un script o función fue modificada para precompilarla (usa la precompilación anterior). Esto
   %hace que las modificaciones del script las ignora y usa por ejemplo valores de variable que son
   %de la versión del script anterior. Esto se arregla con clear all, pero eso puede traer muchos
   %problemas, además que se desaprovecha las precompilaciones previas. Se borra solo la
   %precompilación de la función (hay que usar el nombre de la función solo, sin camino).
   clear(nomArchME)
   %run corre el script sin estar en el path o que en el directorio activo
   run([e_VG.fileCompleto,'.m'])
   %En este script se define:
   %- m_ElemLoc: indicada como una lista de números de elementos del localizado impuesto.
   %- lMicro: es el ancho de la banda localizada en la celda unitaria. Zona donde se aplica la
   %deformación Bxn/lMicro.
   %- angBanda: el ángulo de la banda de localización, que sirve para indentificar cuál ángulo se
   %adopta de los dos devueltos del análisis de bifurcación.
   %- angNormal: el ángulo que se quiere forzar que la normal tenga, descartando la que se obtiene
   %del análisis de bifurcación.
   %- tipoCBL: el tipo de condición de borde del dominio localizado que se impone después de la
   %bifurcación. 0, no se impone CB, 1, se impone mínima restricción, 2, se impone periodicidad (en
   %en este caso habría que definir como imponerla, posiblemente crear una función que delvuelva las
   %condiciones de borde), 3, se impone CB lineal, 4, se impone CB de Taylor.
   %- tipoCBBif: se indica el tipo de CB de toda la celda unitaria después de la bifurcación (se
   %debe pasar una función que devuelva m_CondBord redefinido). No implementado todavía.
   %- pasosPosBif: es la cantidad de pasos después de la bifuración que se quiere imponer el SDA y
   %el modelo multiescala cohesivo.
   %- pasosPosBifCBExt y f_CambioCBExt: Si se quiere cambiar la condición externa de la CU después de la
   %bifurcación se debe indicar la cantidad de pasos después de la misma para hacer el cambio y el handle de
   %la función que devuelve la nuevas condiciones de borde.
   %- longFis: longitud de la curva central al dominio localizado (longitud de fisura).
   %- facNormMicro: Media a largo de la fisura S (línea central al dominio localizado) de la normales micro.

   if exist('m_ElemLoc','var')
      %Se recupera la matriz de denominación de los elementos (ver si no traerla directamente de
      %read_data).
      m_NumElem = zeros(1,e_VG.nElem);
      for iSet = 1:e_VG.nSet
         m_NumElem(e_DatSet(iSet).m_IndElemSet) = e_DatSet(iSet).m_NumElem;
      end
      %Elementos del dominio localizado impuesto.
      %Es importante poner explícitamente 1 en any, para el caso que se tenga un solo elemento que
      %localiza (ya que sino el vector lo analiza como una fila).
      e_DatMatME.m_ElemLocImp = any(bsxfun(@eq,m_ElemLoc',m_NumElem),1);
   end
   if exist('lMicro','var')
      %lMicro impuesto
      if isnumeric(lMicro)&&~isnan(lMicro)
         e_DatMatME.lMicroImp = lMicro;
      else
         %Al imponerse lMicro como un string o como NaN se lee las normales micro y las longitudes micro de un
         %archivo de datos.
         e_DatMatME.lMicroImp = NaN;
      end
   end
   if exist('angBanda','var')
      e_DatMatME.angBanda = angBanda;
   end
   if exist('angNormal','var')
      e_DatMatME.angNormalImp = angNormal;
   end
   if exist('tipoCBL','var')
      e_DatMatME.tipoCBL = tipoCBL;
   else
      %Se adopta que por defecto se impone la condición de borde lineal en el dominio localizado.
      e_DatMatME.tipoCBL = 3;
   end
   if exist('pasosPosBif','var')
      e_DatMatME.pasosPosBif = pasosPosBif;
   else
      %Por defecto son 0 pasos después de la bifurcación.
      e_DatMatME.pasosPosBif = 0;
   end
   %Para el cambio de la condición de borde despúes de bifurcación
   if exist('pasosPosBifCBExt','var')
      e_DatMatME.pasosPosBifCBExt = pasosPosBifCBExt;
   end
   %Función que devuelve la nueva condición de borde externa.
   %En la variable se debe pasar el handle de la función (funciones anónimas por ejemplo).
   if exist('f_CambioCBExt','var')
      e_DatMatME.f_CambioCBExt = f_CambioCBExt;
   else
      if exist('pasosPosBifCBExt','var')
         error(['Lectura de datos: Propiedades Materiales: Modelo multiescala cohesivo: ',...
            'Si se define la variable pasosPosBifCBExt se debe indicar la función que devuelve las nuevas',...
            'condiciones de borde.'])
      end
   end
   %Longitud total de la fisura en la microescala, se la define si se quiere imponer.
   if exist('longFis','var')
      e_DatMatME.longFisImp = longFis;
   end
   %Media de las normales micro a lo largo de la fisura
   if exist('facNormMicro','var')
      e_DatMatME.facNormMicroImp = facNormMicro;
   end
  
end