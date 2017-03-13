function [m_Be,m_DetJe] = f_MatBe_bbar_q1(coord_n,e_DatElemSet,e_VG)

   %xg = e_VG.xg;
   ntens = e_VG.ntens;
   %dofpe = e_VG.dofpe;
   %npg = e_VG.npg;
   struhyp = e_VG.struhyp;
   
   dofpe = e_DatElemSet.dofpe;
   xg = e_DatElemSet.xg;   
   npg = e_DatElemSet.npg;

   m_Be = zeros(ntens,dofpe,npg);
   m_DetJe = zeros(npg,1);
   
   BH = matrixBh_bbar_q1(coord_n,e_DatElemSet,e_VG);
   for iPG = 1:npg
      switch struhyp
         case 1
            [m_Be(:,:,iPG),m_DetJe(iPG)] = matrixB_bbar_q1(coord_n,xg(iPG,:),e_DatElemSet,BH,e_VG);
         otherwise
            error('Elemento bbar_q1: Matriz de deformación B Global: Hipótesis o estado de cálculo no disponible.');
      end
   end
   
end