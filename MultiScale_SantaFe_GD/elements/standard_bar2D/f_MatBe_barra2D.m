function [m_Be,m_DetJe] = f_MatBe_barra2D(coord_n,e_VG)

   xg = e_VG.xg;
   ntens = e_VG.ntens;
   dofpe = e_VG.dofpe;
   npg = e_VG.npg;

   m_Be = zeros(ntens,dofpe,npg);
   m_DetJe = zeros(npg,1);
   
   for iPG = 1:npg
      [m_Be(:,:,iPG),m_DetJe(iPG)] = matrixB_barra2D (coord_n,xg(iPG,:),e_VG);
   end
   
end