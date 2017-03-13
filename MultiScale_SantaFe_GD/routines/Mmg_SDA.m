function [e_VG,e_DatSet] = Mmg_SDA(xx,e_VG,e_DatSet)

M = zeros(e_VG.nnod,1);
for iSet = 1:e_VG.nSet
   nElem     = e_DatSet(iSet).nElem;
   eltype    = e_DatSet(iSet).e_DatElem.eltype;
   npg       = e_DatSet(iSet).e_DatElem.npg;
   npe       = e_DatSet(iSet).e_DatElem.npe;
   wg        = e_DatSet(iSet).e_DatElem.wg;
   xg        = e_DatSet(iSet).e_DatElem.xg;
   conec     = e_DatSet(iSet).conec;
   m_DetJT   = e_DatSet(iSet).m_DetJT;
   m_VolElem = e_DatSet(iSet).m_VolElem;
   switch eltype
      case 4
         m_corrMASA = zeros(npe,npe,nElem);
         N_vector   = zeros(npe,nElem);
         dN_xy      = zeros(e_VG.ndime,npe,nElem);
         for iElem=1:nElem
            coord_n = xx(conec(iElem,:),:)';
            m_pesoPG = m_DetJT(:,iElem).*wg(:);
            m_masa   = zeros(npe);
            for i = 1:npg
               E = xg(i,1);
               n = xg(i,2);
               N = [(1/4)*(1-E)*(1-n) (1/4)*(1+E)*(1-n) (1/4)*(1+E)*(1+n) (1/4)*(1-E)*(1+n)];
               m_masa = m_masa + N'*N*m_pesoPG(i);
               N_vector(:,iElem) = N_vector(:,iElem) + N'*m_pesoPG(i);
            end
            m_corrMASA(:,:,iElem) = inv(m_masa).*(m_VolElem(iElem)/npe);
            for j = 1:length(conec(iElem,:))
               M(conec(iElem,j)) = M(conec(iElem,j))+ (m_VolElem(iElem)/npe);
            end
            %
            E = 0;
            n = 0;
            dN_E = [(1/4*n-1/4)  (-1/4*n+1/4)  (1/4*n+1/4)  (-1/4*n-1/4)];
            dN_n = [(1/4*E-1/4)  (-1/4*E-1/4)  (1/4*E+1/4)  (-1/4*E+1/4)];
            dN_En = [dN_E ; dN_n];
            %
            J11 = coord_n(1,:)*dN_E';
            J12 = coord_n(2,:)*dN_E';
            J21 = coord_n(1,:)*dN_n';
            J22 = coord_n(2,:)*dN_n';
            J = [J11 J12 ; J21 J22];
            dN_xy(:,:,iElem) = J\dN_En;
         end
         e_DatSet(iSet).MElemInv = m_corrMASA;
         e_DatSet(iSet).N_vector = N_vector;
         e_DatSet(iSet).dN_xy    = dN_xy;
      case {21,22,23}
         m_corrMASA = zeros(npe,npe,nElem);
         N_vector   = zeros(npe,nElem);
         dN_xy      = zeros(e_VG.ndime,npe,nElem);
         for iElem=1:nElem
            coord_n = xx(conec(iElem,:),:)';
            m_pesoPG = m_DetJT(:,iElem).*wg(:);
            m_masa   = zeros(npe);
            for i = 1:(npg-2);
               E = xg(i,1); n = xg(i,2);
               N = [(1/4)*(1-E)*(1-n) (1/4)*(1+E)*(1-n) (1/4)*(1+E)*(1+n) (1/4)*(1-E)*(1+n)];
               m_masa = m_masa + N'*N*m_pesoPG(i);
               N_vector(:,iElem) = N_vector(:,iElem) + N'*m_pesoPG(i);
            end
            m_corrMASA(:,:,iElem) = inv(m_masa).*(m_VolElem(iElem)/npe);
            for j = 1:length(conec(iElem,:))
               M(conec(iElem,j)) = M(conec(iElem,j))+ (m_VolElem(iElem)/npe);
            end
            %
            E = 0; n = 0;
            dN_E = [(1/4*n-1/4)  (-1/4*n+1/4)  (1/4*n+1/4)  (-1/4*n-1/4)];
            dN_n = [(1/4*E-1/4)  (-1/4*E-1/4)  (1/4*E+1/4)  (-1/4*E+1/4)];
            dN_En = [dN_E ; dN_n];
            %
            J11 = coord_n(1,:)*dN_E';
            J12 = coord_n(2,:)*dN_E';
            J21 = coord_n(1,:)*dN_n';
            J22 = coord_n(2,:)*dN_n';
            J = [J11 J12 ; J21 J22];
            dN_xy(:,:,iElem) = J\dN_En;
         end
         e_DatSet(iSet).MElemInv = m_corrMASA;
         e_DatSet(iSet).N_vector = N_vector;
         e_DatSet(iSet).dN_xy    = dN_xy;
      otherwise
         error('Lectura de datos: Cálculo de matrices de masas: Elemento finito no definido.')
   end
end
e_VG.MGlobInv=1./M;

