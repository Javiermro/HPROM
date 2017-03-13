function e_DatSet = f_MatBT(xx,e_DatSet,e_VG)
   
   ntens = e_VG.ntens;
   ndime = e_VG.ndime;
   nSet = e_VG.nSet;
   %struhyp = e_VG.struhyp;
   
   for iSet = 1:nSet
      
      conec = e_DatSet(iSet).conec;
      nElem = e_DatSet(iSet).nElem;
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      eltype = e_DatElemSet.eltype;
      npe = e_DatElemSet.npe;
      dofpe = e_DatElemSet.dofpe;
      npg = e_DatElemSet.npg;
      wg = e_DatElemSet.wg;
      
      %La evaluación de las funciones de forma en los puntos de gauss no depende de la forma del elemento, por
      %lo que se evalúa para todos los puntos de gauss del elemento máster de todos los elementos, ahorrando
      %memoria.
      switch eltype
         case {2,10,32}
            m_FF = f_FF_tria_t1(e_DatElemSet,e_VG);
         case {4,8,20,21,22,23,31,108}
            m_FF = f_FF_quad_q1(e_DatElemSet,e_VG);
         otherwise
            error('Matrices de función de forma: Elemento no implementado.')
      end
      
      %
      m_BT = zeros(ntens,dofpe,npg,nElem);
      m_DetJT = zeros(npg,nElem);
      if eltype==31 
         dN_xy = zeros(ndime,npe,nElem);
      end
      
      %parfor iElem = 1:nElem
      for iElem = 1:nElem
         %Para que funcione correctamente con los índices, el parfor se hace en una función separada 
         %que resuelva los loops sobre los puntos de Gauss

         coord_n = f_CoordElem(xx,conec(iElem,:));         

         switch eltype
            case {2,10}
               [m_BT(:,:,:,iElem),m_DetJT(:,iElem)] = f_MatBe_tria_t1(...
                  coord_n,e_DatElemSet,e_VG);
            case {4,20,21,22,23,108}
               [m_BT(:,:,:,iElem),m_DetJT(:,iElem)] = f_MatBe_quad_q1(...
                  coord_n,e_DatElemSet,e_VG);
            case 5
               [m_BT(:,:,:,iElem),m_DetJT(:,iElem)] = f_MatBe_barra2D(coord_n,e_DatElemSet,e_VG);
            case 7
               error('Matrices de deformación: Falta implementar el caso de hexaedros de 8 nodos')
            case 8
               [m_BT(:,:,:,iElem),m_DetJT(:,iElem)] = f_MatBe_bbar_q1(coord_n,e_DatElemSet,e_VG);
            case 31
               [m_BT(:,:,:,iElem),m_DetJT(:,iElem)] = ...
                   f_MatBe_quad_q1(coord_n,e_DatElemSet,e_VG);
               E = 0; n =0;
               dN_E = [(1/4*n-1/4)  (-1/4*n+1/4)  (1/4*n+1/4)  (-1/4*n-1/4)];
               dN_n = [(1/4*E-1/4)  (-1/4*E-1/4)  (1/4*E+1/4)  (-1/4*E+1/4)];
               dN_En = [dN_E ; dN_n];
               J11 = coord_n(1,:)*dN_E';
               J12 = coord_n(2,:)*dN_E';
               J21 = coord_n(1,:)*dN_n';
               J22 = coord_n(2,:)*dN_n';
               J = [J11 J12 ; J21 J22];
               dN_xy(:,:,iElem) = J\dN_En;
            case 32
               [m_BT(:,:,:,iElem),m_DetJT(:,iElem)] = f_MatBe_tria_t1(coord_n,e_DatElemSet,e_VG);
               dN_E = [-1 1 0];
               dN_n = [-1 0 1];
               dN_En = [dN_E ; dN_n];
               J11 = coord_n(1,:)*dN_E';
               J12 = coord_n(2,:)*dN_E';
               J21 = coord_n(1,:)*dN_n';
               J22 = coord_n(2,:)*dN_n';
               J = [J11 J12 ; J21 J22];
               dN_xy(:,:,iElem) = J\dN_En;
            otherwise
               error('Matrices de deformación: Elemento no implementado.')
         end
      end
      %
      e_DatSet(iSet).m_BT = m_BT;
      e_DatSet(iSet).m_DetJT = m_DetJT;
      e_DatSet(iSet).m_FF = m_FF;
      %Volumen de los elementos
      %El punto de gauss ya viene multiplicado por el espesor.
      %Para considerar que el elemento SDA_tria_t1 tiene dos PG ubicados en la misma posición para
      %cálculo del volumen del elemento se debe considerar solo de ellos (ver si no poner que el
      %punto de gauss singular tenga peso nulo).
      %Parecido ocurre con el elemento MixStrInj_quad_q1, donde el PG 5 tiene el peso como si fuera un único
      %PG del elemento, por lo que para determinar el volumen se utiliza los 4 primeros PGs.
      switch eltype
          case {2,4,7,8,108}
              e_DatSet(iSet).m_VolElem = wg'*m_DetJT;
              e_DatSet(iSet).m_VolElemGP = bsxfun(@times,wg,m_DetJT);
          case 10
              e_DatSet(iSet).m_VolElem = wg(1)'*m_DetJT(1,:);
              e_DatSet(iSet).m_VolElemGP = bsxfun(@times,wg(1),m_DetJT(1,:));
           case {20,21,22,23}
              e_DatSet(iSet).m_VolElem = wg(1:4)'*m_DetJT(1:4,:);
              e_DatSet(iSet).m_VolElemGP = bsxfun(@times,wg(1:4),m_DetJT(1:4,:));
          case {31,32}
              e_DatSet(iSet).m_VolElem = wg'*m_DetJT;
              e_DatSet(iSet).m_VolElemGP = bsxfun(@times,wg,m_DetJT);
              e_DatSet(iSet).dN_xy = dN_xy;
          otherwise
            error('Volumen del elemento: Elemento no implementado.')
      end
      
   end

end
