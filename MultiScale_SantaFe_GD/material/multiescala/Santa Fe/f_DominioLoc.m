function [m_ElemLoc,c_ProyIncrEps] = f_DominioLoc(e_VarEst_new,e_VarEst_old,m_DefLocal,e_DatSet,e_VG)

   nSet = e_VG.nSet;
   nElem = e_VG.nElem;
   ntens = e_VG.ntens;
   m_ElemLoc = false(1,nElem);
   %Determinación del incremento en la "dirección" de Bxn (tensor)
   %Se divide por 2 para la doble contracción 
   m_DefLocContr = e_VG.FOAT1*m_DefLocal;
   %Se normaliza el tensor proyección (m_DefLocContr ya tiene el término de corte dividido 2)
   m_DefLocContr = m_DefLocContr/sqrt((e_VG.SONT*m_DefLocContr)'*m_DefLocContr);
   %
   c_ProyIncrEps = cell(nSet,1);
   for iSet = 1:nSet
      npg = e_DatSet(iSet).e_DatElem.npg;
      %Incremento de deformaciones totales
      %m_IncrEps = e_VarHist(iSet).eps_new-e_VarHist(iSet).eps_old;
      % Incremento de deformaciones fluctuantes
      m_IncrEps = reshape(e_VarEst_new(iSet).eps_fluct-e_VarEst_old(iSet).eps_fluct,ntens,npg,[]);
      % Se adopta la media del incremento de las fluctuaciones de los 4 PGs.
      m_IncrEps = reshape(mean(m_IncrEps,2),ntens,[]);
      % Normalización del incremento de las fluctuaciones
      %m_IncrEps = bsxfun(@rdivide,m_IncrEps,sqrt(sum((e_VG.FOAT1*m_IncrEps).*m_IncrEps,1)));
      % Proyección sobre deformación de localización
      m_ProyIncrEps = sum(bsxfun(@times,m_IncrEps,m_DefLocContr),1);
%       % Condición de PG bifurcado
%       %Cuidado este máximo no tiene sentido porque es para cada set, y que tiene que ser máximo de todos los
%       %sets.
%       %m_PGBif = m_ProyIncrEps>0.01*max(m_ProyIncrEps(:));
%       m_PGBif = m_ProyIncrEps>0;
%       %Criterio de elemento localizado
%       %m_ElemLoc(e_DatSet(iSet).m_IndElemSet) = any(m_PGBif,2);  %all(m_PGBif,2)
%       m_ElemLoc(e_DatSet(iSet).m_IndElemSet) = m_PGBif;
      %
      %c_ProyIncrEps{iSet} = reshape(m_ProyIncrEps,npg,[]);
      c_ProyIncrEps{iSet} = m_ProyIncrEps;
   end
   
   %No interesa el orden con que se concatena ya que se realiza una media de los valores.
   m_ProyIncrEpsComp = [c_ProyIncrEps{:}];
   m_IndProyIncrEpsCompPos = m_ProyIncrEpsComp>0;
   % Se multiplica por el área del elemento para hacer una media ponderada.
   m_VolElemComp = [e_DatSet.m_VolElem];   
   %m_ProyIncrEpsComp = [c_ProyIncrEps{:}].*[e_DatSet.m_VolElem];
   % Media de las proyecciones positivas
   %medProyIncrEps = mean(m_ProyIncrEpsComp(m_ProyIncrEpsComp>0));
   medProyIncrEps = sum(m_ProyIncrEpsComp(m_IndProyIncrEpsCompPos).*m_VolElemComp(m_IndProyIncrEpsCompPos))/...
      sum(m_VolElemComp(m_IndProyIncrEpsCompPos));
   %
   for iSet = 1:nSet
      m_ProyIncrEps = c_ProyIncrEps{iSet};
      %Criterio de elemento localizado
      m_PGBif = m_ProyIncrEps>0.05*medProyIncrEps;
      m_ElemLoc(e_DatSet(iSet).m_IndElemSet) = m_PGBif;
   end
   
end