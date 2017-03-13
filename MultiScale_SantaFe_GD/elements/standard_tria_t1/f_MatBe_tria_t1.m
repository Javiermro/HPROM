function [m_Be,m_DetJe] = f_MatBe_tria_t1(coord_n,e_DatElemSet,e_VG)

   ntens = e_VG.ntens;

   struhyp = e_VG.struhyp;
   
   dofpe = e_DatElemSet.dofpe;
   xg = e_DatElemSet.xg;   
   npg = e_DatElemSet.npg;

   m_Be = zeros(ntens,dofpe,npg);
   m_DetJe = zeros(npg,1);
   
   for iPG = 1:npg
      switch struhyp
         case {1,2}
            [m_Be(:,:,iPG),m_DetJe(iPG),~] = matrixB_tria_t1_epd(coord_n,e_DatElemSet,e_VG);
         case 3
            error('Elemento tria_t1: Matriz de deformación B Global: Estado 3D no disponible.');
         case 4
            [m_Be(:,:,iPG),m_DetJe(iPG)] = matrixB_tria_t1_axm(coord_n,xg(iPG,:),e_DatElemSet,e_VG);
         case 20   %Large deformations with plane deformation.
            [m_Be(:,:,iPG),m_DetJe(iPG)] = matrixB_tria_t1_epd_ld(coord_n,e_DatElemSet,e_VG);
         otherwise
            error('Elemento tria_t1: Matriz de deformación B Global: Estado no disponible.');
      end
   end
   
end