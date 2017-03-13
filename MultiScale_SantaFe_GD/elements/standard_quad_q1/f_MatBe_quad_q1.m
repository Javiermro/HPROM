function [m_Be,m_DetJe] = f_MatBe_quad_q1(coord_n,e_DatElemSet,e_VG)

   %xg = e_VG.xg;
   ntens = e_VG.ntens;
   %ndn = e_VG.ndn;
   struhyp = e_VG.struhyp;

   dofpe = e_DatElemSet.dofpe;
   xg = e_DatElemSet.xg;   
   npg = e_DatElemSet.npg;
      
   m_Be = zeros(ntens,dofpe,npg);
   m_DetJe = zeros(npg,1);
   
   for iPG = 1:npg
      switch struhyp
         case {1,2}
            [m_Be(:,:,iPG),m_DetJe(iPG)] = matrixB_quad_q1_epd(coord_n,xg(iPG,:),...
               e_DatElemSet,e_VG);
         case 3
            fprintf('Estado 3D no disponible');
         case 4
            [m_Be(:,:,iPG),m_DetJe(iPG)] = matrixB_quad_q1_axm(coord_n,xg(iPG,:),e_DatElemSet,e_VG);
         case 20   %Large deformations with plane deformation.
            [m_Be(:,:,iPG),m_DetJe(iPG)] = matrixB_quad_q1_epd_ld(coord_n,xg(iPG,:),e_DatElemSet,e_VG);
         otherwise
            error('Elemento quad_q1: Matriz de deformación B Global: Estado no disponible.');
      end
   end
   
end