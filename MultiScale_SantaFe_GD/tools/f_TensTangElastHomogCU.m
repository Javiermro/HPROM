function m_Ce = f_TensTangElastHomogCU(archDat)

   %Esta programa está hecho para obtener el tensor tangente elástico de una celda unitaria sin entrar al
   %programa (en forma independiente).
   %Se llama directamente el archivo de datos de la celda unitaria.
   
   %Se asume que se está corriendo esta función con todos los directorios del programa agregado   
   %al path.
   
   if ~exist('archDat','var')
      [nomArchDat,pathArchDat] = uigetfile('*.mfl','Definir archivo de cálculo');
   else
      [pathArchDat,nomArchDat,ext] = fileparts(archDat);
      nomArchDat = [nomArchDat,ext];
   end
   %
   [~,m_Coord,~,~,~,e_DatSet,e_VG] = read_data(nomArchDat,pathArchDat);
   %
   [c_CT,KT] = f_TensTangElastCU(e_DatSet,e_VG);
   %
   nElem = e_VG.nElem;   
   omegaMicro = sum([e_DatSet.m_VolElem]);
   [m_LinCond,~,~,doff,dofl] = f_CondBord(e_VG,m_Coord,e_DatSet,e_VG.m_ConecFront);
   m_Ce = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,...
      e_DatSet,omegaMicro,true(nElem,1),true(nElem,1),e_VG);   
   
end