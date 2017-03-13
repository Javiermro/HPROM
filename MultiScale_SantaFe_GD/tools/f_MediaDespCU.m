function m_uMedio = f_MediaDespCU(u,omegaMicro,e_DatSet,e_VG)

   nSet = e_VG.nSet;
   ndn = e_VG.ndn;
   %El desplazamiento es un campo nodal, por lo que para intregrar sobre el dominio se lo lleva a los punto de
   %gauss.
   c_DespElem = cell(nSet,1);

   %
   for iSet = 1:nSet
      %
      nElem = e_DatSet(iSet).nElem;
      nPG = e_DatSet(iSet).e_DatElem.npg;
      m_DofElem = e_DatSet(iSet).m_DofElem;
      m_FF = e_DatSet(iSet).m_FF;
      %
      m_uElemSet = reshape(u(m_DofElem),[],nElem);
      %
      m_DespPG = zeros(ndn,nPG,nElem);
      for iElem = 1:nElem
         %squeeze llama reshape, así que no es más rápido que esta si se conoce cuál es la dimensión de la
         %matriz con valor 1.
         m_DespPG(:,:,iElem) = squeeze(sum(bsxfun(@times,m_FF(:,:,:),m_uElemSet(:,iElem)'),2));
      end
      %
      c_DespElem{iSet} = m_DespPG;
   end
   
   m_uMedio = f_HomogArea(c_DespElem,ndn,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VG);

end