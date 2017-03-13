function e_VarAux =...
   SDA_Properties(xx,e_VG,e_VarEst_new,e_DatSet,e_VarAux)

% Calculo del gradiente de Phi y seleccion de los nodos solitarios
% ****************************************************************

%NORMAL VECTOR SMOOTHING
for iSet = 1:e_VG.nSet
   e_DatElemSet  = e_DatSet(iSet).e_DatElem;
   eltype        = e_DatElemSet.eltype;
   
   switch eltype
      case {21,22,23}
         conec            =  e_DatSet(iSet).conec;
         nElem            =  e_DatSet(iSet).nElem;
         npe              =  e_DatElemSet.npe;
         pointersVAE      =  e_DatElemSet.pointersVAE ;
         p_phi_grad       =  pointersVAE.p_phi_grad  ;
         p_n_tens         =  pointersVAE.p_n_tens    ;
         p_fii            =  pointersVAE.p_fii       ;
         i_indST          =  e_DatElemSet.pointersVHE.i_indST    ;
         m_VarAuxElem     =  e_VarAux(iSet).VarAuxElem;
         m_VarHistElemNew =  e_VarEst_new(iSet).VarHistElem;
         dN_MACRO         =  e_DatSet(iSet).dN_xy  ;
         
         %  p_EFsidesCutCPF = pointersVAE.p_EFsidesCutCPF ;
         %  EF_sidesCutCPF  = m_VarAuxElem(p_EFsidesCutCPF,:);
         
         EF_sidesCutCPF   =  e_DatSet(iSet).e_DatElem.EF_sidesCutCPF;
         
         
         for iElem = 1:nElem
            
            m_phi_grad      =  m_VarAuxElem(p_phi_grad,iElem) ;
            m_phi_grad      =  reshape(m_phi_grad,4,2);
            m_n_tens        =  m_VarAuxElem(p_n_tens,iElem) ;
            m_n_tens        =  reshape(m_n_tens,4,2);
            m_fii           =  m_VarAuxElem(p_fii,iElem) ;
            m_indSTmacro    =  m_VarHistElemNew (i_indST ,iElem)  ;
            coord_n         =  xx(conec(iElem,:),:)';
            
            [m_phi_grad,m_fii] = ...
               node_selection_quad_q1_SDA2...
               (coord_n,npe,m_n_tens,m_phi_grad,dN_MACRO(:,:,iElem),...
               m_indSTmacro,m_fii,EF_sidesCutCPF(iElem));
            
            m_phi_grad                      = reshape(m_phi_grad,8,1);
            m_VarAuxElem(p_phi_grad,iElem)  = m_phi_grad  ;
            m_VarAuxElem(p_fii,iElem)       = m_fii ;
            
         end
         %
         e_VarAux(iSet).VarAuxElem =  m_VarAuxElem;
   end
end
