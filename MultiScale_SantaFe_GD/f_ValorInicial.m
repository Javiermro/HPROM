function uValIn = f_ValorInicial(xx,m_LinCond,dofl,doff,u,c_GdlCond,Fext,...
   e_VarEst_new,e_VarEst_old,e_VarAux,e_DatSet,DefMacro,e_VG)

   ndoft = e_VG.ndoft;
   
   %El residuo se calcula correctamente (considerando la posibilidad de comportamiento inelástico de
   %los materiales)
   [~,FintInelast] = f_Par_MatGlobales(xx,u,c_GdlCond,e_DatSet,e_VarEst_new,e_VarEst_old,...
      e_VarAux,DefMacro,e_VG);   
   %Se calcular la matriz de rigidez global forzando que los materiales tengan un comportamiento 
   %elástico (considerando la historia de daño, por ejemplo, que tenga el material en ese momento).
   e_VG.elast = 1;
   KT = f_Par_MatGlobales(xx,u,c_GdlCond,e_DatSet,e_VarEst_new,e_VarEst_old,e_VarAux,DefMacro,e_VG);
   e_VG.elast = 0;

   % ITERACION "0", COMPUTO DEL RESIDUO
   res = residuo(Fext,FintInelast,m_LinCond,dofl,doff,1);

   % INCREMENTO DE DESPLAZAMIENTO 
   %Como se resuelve elásticamente se exige que no use ninguna estrategia de control.
   e_VG.CONTROL_STRAT = 1;
   du_iterl = incremental_disp(res,KT,Fext,m_LinCond,dofl,doff,[],[],e_VG);
   du_iter = zeros(ndoft,1);
   du_iter(dofl) = du_iterl;
   du_iter(doff) = m_LinCond*du_iterl;
   
   % DESPLAZAMIENTO TOTAL
   uValIn = u+du_iter;

end