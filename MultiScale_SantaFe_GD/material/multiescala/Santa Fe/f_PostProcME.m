function f_PostProcME(iStep,m_NumElemMacro,m_DefMacro,e_DatMatSetMacro,hvar_newMacro,e_VGMacro)


   %Se recupera variables micro
   xx = e_DatMatSetMacro.xx;
   in = e_DatMatSetMacro.in;
   e_VG = e_DatMatSetMacro.e_VG;
   m_ElemPGImpr = e_DatMatSetMacro.m_ElemPGImpr;
   e_DatSet = e_DatMatSetMacro.e_DatSet;
   m_SetElem = e_DatMatSetMacro.m_SetElem;
   %Se guarda el paso tiempo a nivel micro.
   e_VG.istep = e_VGMacro.istep;
   %
   postpro_impre_step = e_VG.postpro_impre_step;
   IRES = e_VG.IRES;
   ndime = e_VG.ndime;
   nElem = e_VG.nElem;
   %Nombre completo del archivo micro
   nomArchCompl = e_VG.fileCompleto;
   %   
   for iPGGraf = 1:size(m_ElemPGImpr,2)
      
      iElemSetMacro = m_ElemPGImpr(1,iPGGraf);
      iPGMacro = m_ElemPGImpr(2,iPGGraf);
      %
      %m_iDefMacro = m_DefMacro(:,iPGMacro,iElemSetMacro);
      m_iDefMacro = m_DefMacro(:,iElemSetMacro);
      e_VarEst_new = hvar_newMacro(iPGMacro,iElemSetMacro).e_VarEst;
      e_VarAux = hvar_newMacro(iPGMacro,iElemSetMacro).e_VarAux;
      u = hvar_newMacro(iPGMacro,iElemSetMacro).u;
      c_GdlCond = hvar_newMacro(iPGMacro,iElemSetMacro).c_GdlCond;
      Fint = hvar_newMacro(iPGMacro,iElemSetMacro).Fint;
      m_DefMacroCU = hvar_newMacro(iPGMacro,iElemSetMacro).c_DefMacro;
      %
      m_ElemLoc = false(1,nElem);
      %Se cambia el nombre del archivo de postproceso de GiD de la celda unitaria, para
      %individualizar el elemento y el punto de Gauss que corresponde.
      e_VG.fileCompleto = [nomArchCompl,'_EM',num2str(m_NumElemMacro(iElemSetMacro)),'PGM',...
         num2str(iPGMacro)];

      % IMPRESIÓN DE RESULTADOS
      %Notar que se imprime este paso si a nivel macro este paso se imprime, y a la vez, que si a
      %nivel micro debe imprimirse.
      if IRES   
         if ~mod(iStep,postpro_impre_step)
            % DESPLAZAMIENTO TOTAL DE LA MICRO-CELDA
            uTotal = [m_iDefMacro(1),m_iDefMacro(4)/2;m_iDefMacro(4)/2,m_iDefMacro(2)]*xx(:,1:2)'...
               +reshape(u,ndime,[]);
            % POSTPROCESO GiD
            matlab2gid_res(iStep,in,u,c_GdlCond,e_DatSet,e_VarEst_new,e_VarAux,m_DefMacroCU,...
               uTotal,m_ElemLoc,e_VG)
         end
      end
      %Impresión gráficas
      f_SaveDatGraf(u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,e_DatSet,m_SetElem,[],[],e_VG) %JLM
      
   end
         
end