function e_DatMatME = f_VarME(nomArchME,dirDat,m_ElemPGImpr)

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

   % FUNCIÓN CONDICIONES DE BORDE
   %En es caso del modelo multiescala clásico se puede mantener las matrices de condiciones de borde
   %en la estructura de propiedades e_DatMat macro ya que no cambian durante el problema.
%    [m_LinCond,vfix,m_InvCRR,doff,dofl,doffCondCte] = f_CondBord(e_VG,xx,...
%      e_DatSet,e_VG.m_ConecFront);

   % Script con variables 
   %(ver si correlo acá, que significa que hay que llevar varias variables o en cada paso tiempo, donde habría
   %que analizar cuánto tiempo tarda)
   %Ocurre un error al usar usar run, no sé si pasa lo mismo con eval, ya que matlab no se da cuenta que un
   %script o función fue modificada para precompilarla (usa la precompilación anterior). Esto hace que las
   %modificaciones del script las ignora y usa por ejemplo valores de variable que son de la versión del
   %script anterior. Esto se arregla con clear all, pero eso puede traer muchos problemas, además que se
   %desaprovecha las precompilaciones previas. Se borra solo la precompilación de la función (hay que usar el
   %nombre de la función solo, sin camino).
   if exist([e_VG.fileCompleto,'.m'],'file')
      clear(nomArchME)
      %run corre el script sin estar en el path o que en el directorio activo
      run([e_VG.fileCompleto,'.m'])
   end
   
   % Cálculo de la medida de la celda unitaria
   %Cómo se hace operaciones sobre el volumen, se guarda el volumen de los elementos en una matriz.
   %Se la ordena según la numeración (orden) de la matriz de conectividad global (como fueron 
   %ingresados los elementos).
   m_VolElem = zeros(1,e_VG.nElem);
   m_VolElem([e_DatSet.m_IndElemSet]) = [e_DatSet.m_VolElem];
   %En caso que se ingrese en el script el área micro, se utiliza directamente esa.
   if ~exist('Omega_micro','var')
      %Esta expresión es correcta solo si la celda unitaria no tiene agujeros.
      Omega_micro = sum(m_VolElem);
   end
   
   e_DatMatME = struct('in',in,'xx',xx,'m_SetElem',m_SetElem,'f',f,'funbc',funbc,...
      'e_DatSet',e_DatSet,'m_ElemPGImpr',m_ElemPGImpr,'omegaMicro',Omega_micro,'e_VG',e_VG);
  
   %Imprensión de archivo de postprocesado de la malla y inicialización del archivo de datos
   %matlab2gid_mesh(in,xx,e_DatSet,e_VG)
   %f_InicArchDat(e_VG)
     
end