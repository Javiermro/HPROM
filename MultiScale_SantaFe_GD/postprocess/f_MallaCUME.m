function f_MallaCUME(m_NumElemMacro,e_DatMatSetMacro)

   %Se recupera variables micro
   xx = e_DatMatSetMacro.xx;
   in = e_DatMatSetMacro.in;
   e_VG = e_DatMatSetMacro.e_VG;
   m_ElemPGImpr = e_DatMatSetMacro.m_ElemPGImpr;
   e_DatSet = e_DatMatSetMacro.e_DatSet;
   %
   nomArchCompl = e_VG.fileCompleto;
   %   
   for iPGGraf = 1:size(m_ElemPGImpr,2)
      
      iElemSet = m_ElemPGImpr(1,iPGGraf);
      iPG = m_ElemPGImpr(2,iPGGraf);
      %Se cambia el nombre del archivo de malla del GiD de la celda unitaria, para
      %individualizar el elemento y el punto de Gauss que corresponde.
      e_VG.fileCompleto = [nomArchCompl,'_EM',num2str(m_NumElemMacro(iElemSet)),'PGM', num2str(iPG)];      
      matlab2gid_mesh(in,xx,e_DatSet,e_VG)
      
   end
   
end
