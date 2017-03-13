function GFMedia = f_GFMedia(volTotal,e_DatSet,e_VG)

   %% Media de la energía de fractura de la microcelda.

   c_GF = cell(nSet,1);
   for iSet = 1:nSet
      
      e_DatMatSet = e_DatSet(iSet).e_DatMat;
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      nPG = e_DatElemSet.npg;
      switch conshyp
         case 1
            m_Gf = zeros(nPG,nElem);            
         case 11
            m_Gf = e_DatMatSet.gfv*ones(nPG,nElem);   
         otherwise
            error('Energía de fractura media: Modelo constitutivo no definido.\n')         
      end
      c_GF{iSet} = m_Gf;
      
   end
   %
   GFMedia = f_HomogArea(c_GF,1,volTotal,{e_DatSet.m_DetJT},e_DatSet,e_VG);

end