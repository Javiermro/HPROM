function m_VarFluc = f_VarFluct(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,m_ElemLocPert,e_VG)

   %Número de grados de libertad por nodo.
   %ndn = e_VG.ndn;
   ndoft = e_VG.ndoft;
   ntens = e_VG.ntens;
   nSet = e_VG.nSet;
   
   m_FCHom = zeros(ndoft,ntens);

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
      dofpe = e_DatElemSet.dofpe;


      m_ElemLocSet = m_ElemLocPert(e_DatSet(iSet).m_IndElemSet);

      for iElem = 1:nElem        
         
         %Se ensambla el elemento si está en el dominio localizado.
         if m_ElemLocSet(iElem)
            m_pesoPG = m_DetJT(:,iElem).*wg;

            m_FCHomElem = zeros(dofpe,ntens);
            for iPG = 1:npg
               %No se coloca el signo menos porque en incremental_disp cuando se resuelve el sistema ya
               %se lo pone.
               m_FCHomElem = m_FCHomElem+m_BT(:,:,iPG,iElem)'*m_CT(:,:,iPG,iElem)*m_pesoPG(iPG);
            end
            %dofElem = f_DofElem(conec(iElem,:),ndn);
            dofElem = m_DofElem(:,iElem);
            m_FCHom(dofElem,:) = m_FCHom(dofElem,:)+m_FCHomElem;
         end
         
      end
      
   end
   
%    m_FCHom = sparse(m_Fil(:),m_Col(:),m_FCHom(:),ndoft,ntens);   
   %No se usa la función residuo porque no maneja el caso del vector de fuerzas con muchas columnas.
   %Este problema es directo, es decir no se está calculando diferencia de fuerzas (residuo), sino
   %que son "fuerzas" aplicadas al problema (condiciones de borde de tipo fuerzas).
   m_FCHom = m_FCHom(dofl,:)+m_LinCond'*m_FCHom(doff,:);
   
   %Determinación del incremento o variación de las fluctuaciones.
   m_VarFluc = zeros(ndoft,ntens);
   m_VarFluc(dofl,:) = incremental_disp(m_FCHom,KT,0,m_LinCond,dofl,doff,[],[],e_VG);
   %En los grados de libertad restringidos no hay variación de las fluctuaciones, es decir como si
   %las condiciones de borde constantes para este problema fuera siempre nulas. No aporta a ecuación
   %siguiente al ser C = 0 (uR = L*uP+C).
   m_VarFluc(doff,:) = m_LinCond*m_VarFluc(dofl,:);   
   
end