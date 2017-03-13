%function m_CTHomog = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,Omega_micro,e_VG)
function m_CTHomog = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,Omega_micro,...
   m_ElemLocHomog,m_ElemLocPert,e_VG)
   %Ver que si no conviene por velocidad tener un función separada cuando se homogeniza en el
   %dominio localizado.
   
   % DETERMINACIÓN DEL MODULO TANGENTE HOMOGENEIZADO
   ntens = e_VG.ntens;
   nSet = e_VG.nSet;
   %ndn = e_VG.ndn;
   
   % Variación de las fluctuaciones de los desplazamientos micros
   %Creo varias operaciones pueden ser evitadas, si se junta la siguiente función con esta, 
   %principalmente si se vectorizan las operaciones (ambas son integraciones sobre los elementos, 
   %pero que se debe realizar en forma separada, ya que hay resolver un sistema para obtener primero
   %los incrementos de las fluctuaciones m_VarFluc).
   m_VarFluc = f_VarFluct(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,m_ElemLocPert,e_VG);

   % Homogenización
   m_CTHomog = zeros(ntens,ntens);
   %Matriz identidad tensorial (se usa para ver si acelera el cálculo, al sacar como factor común
   %m_CT(:,:,iPG,iElem)).
   m_IndTens = eye(ntens);
   %m_IndTens = m_TensProy;
   for iSet = 1:nSet
      
      nElem = e_DatSet(iSet).nElem;
      %conec = e_DatSet(iSet).conec;
      m_DofElem = e_DatSet(iSet).m_DofElem;
      m_DetJT = e_DatSet(iSet).m_DetJT;
      m_BT = e_DatSet(iSet).m_BT;
      m_CT = c_CT{iSet};
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      npg = e_DatElemSet.npg;
      wg = e_DatElemSet.wg;
      
      m_ElemLocSet = m_ElemLocHomog(e_DatSet(iSet).m_IndElemSet);
      
      for iElem = 1:nElem
         %Se ensambla el elemento si está en el dominio localizado.
         if m_ElemLocSet(iElem)
            %dofElem = f_DofElem(conec(iElem,:),ndn);
            dofElem = m_DofElem(:,iElem);
            m_pesoPG = m_DetJT(:,iElem).*wg;
            m_VarFlucElem = m_VarFluc(dofElem,:);
            for iPG = 1:npg
               m_CTHomog = m_CTHomog+m_CT(:,:,iPG,iElem)*(m_IndTens+...
                  m_BT(:,:,iPG,iElem)*m_VarFlucElem)*m_pesoPG(iPG);
            end
         end
      end
      
   end
   
   m_CTHomog = m_CTHomog/Omega_micro;
   
end