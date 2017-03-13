function [sigma,eps,eps_fluct,hvar,m_VarAuxElem] = f_OperPosConvSDA_tria_t1_dano(...
   sigma,eps,eps_fluct,hvar,m_VarAuxElem,m_CT,nElem,m_BT,e_DatMatSet,e_DatElemSet,m_NumElem,e_VG)

   ntens = e_VG.ntens;
   %Se considera que este elemento suele puede tener 2 puntos de gauss, uno para la parte
   %regular y otra para la singular.
   %m_NodSolitElem = [0,0,2,3,0,0];
   sihvarpg = e_DatMatSet.sihvarpg; 
   for iElem = 1:nElem
      %
      % Operaciones que se realiza hasta que se detecta la bifurcación.
      condBif = m_VarAuxElem(1,iElem);
      if ~condBif
         %Para el análisis de bifurcación se toma el tensor de punto de gauss 
         %regular (justamente es único que interesa).
         [condBif,m_angBif] = f_CondBifct(m_CT(:,:,1,iElem),e_VG);
         if condBif==1
            %Para determinar el ángulo que corresponde con la normal de los dos
            %posibles que devuelve el análisis de bifurcación de la discontinuidad
            %fuerte se utiliza la dirección de la tensión principal 1 (mayor valor
            %con signo), que en el caso que la bifurcación se produzca en comprensión
            %sería incorrecto (ver qué hacer?).
            %Parece que no sirve este criterio porque las direcciones devueltas por
            %el análisis de bifurcación son generalmente simétricas con respecto a la
            %dirección principal 1.
            %m_angTensPpal = cribis_sp(sigma(1:ntens,iElem));
            %m_angTensPpal = -15/180*pi;
            m_angTensPpal = -0.05;
            %
            vecBif = [cos(m_angBif);sin(m_angBif)];
            vecDirPpal = [cos(m_angTensPpal);sin(m_angTensPpal)];
            [~,ind] = max(abs(vecBif'*vecDirPpal));
            m_angBif = m_angBif(ind);
            m_Normal = [cos(m_angBif);sin(m_angBif)];
            %m_Normal = [1;0];
            %Se toma el dominio omega_phi+ como el que que contiene el nodo único
            %del elemento triangular. La normal de la fisura se cambia de sentido
            %para apuntar a ese nodo.
            [m_Normal,m_GradPhi] = get_solitary_node(m_Normal,m_BT(:,:,1,iElem),...
               e_DatElemSet,e_VG);
%             B = m_BT(:,:,1,iElem);
%             dN_xy = zeros(e_VG.ndime,e_DatSet(iSet).e_DatElem.npe);
%             dN_xy(1,:) = B(1,1:e_VG.ndn:e_DatSet(iSet).e_DatElem.dofpe);
%             dN_xy(2,:) = B(2,2:e_VG.ndn:e_DatSet(iSet).e_DatElem.dofpe);
%             nodSolit = m_NodSolitElem(e_DatSet(iSet).m_IndElemSet(iElem));
%             m_GradPhi = dN_xy(:,nodSolit);
%             if m_Normal'*m_GradPhi<0
%                m_Normal = -m_Normal;
%             end
            %Se guarda en notación de Voigt el vector normal.
            m_VarAuxElem(1,iElem) = condBif;
            m_VarAuxElem(2:9,iElem) = reshape(f_Vec2Voigt2D(m_Normal,e_VG),[],1);
            m_VarAuxElem(10:17,iElem) = reshape(f_Vec2Voigt2D(m_GradPhi,e_VG),[],1);
            m_VarAuxElem(18,iElem) = e_VG.istep;
%             if e_DatElemSet.simetrico
%                m_VarAuxElem(19:26,iElem) = reshape(f_Vec2Voigt2D(...
%                   m_GradPhi/norm(m_GradPhi),e_VG),[],1);
%             else
%                m_VarAuxElem(19:26,iElem) = m_VarAuxElem(2:9,iElem);
%             end
            %
            fprintf('** Se detectó bifurcación en el elemento %d.\n',m_NumElem(iElem))
         end
      end
      %
      %Se separa los condicionales por si se quiere implementar la discontinuidad
      %fuerte desde el paso que se detectó la bifurcación.
      % Operaciones que ocurren una vez cuando se activa el SDA, en el mismo paso o ciertos pasos
      % después que se detecta la bifurcación.
      condBif = m_VarAuxElem(1,iElem);
      pasoBif = m_VarAuxElem(18,iElem);
      noActBif = m_VarAuxElem(21,iElem);
      if condBif==1&&e_VG.istep==pasoBif+1&&~noActBif
         condBif = 2;
         m_VarAuxElem(1,iElem) = condBif;
         %Se copia las variables históricas del PG regular al singular.
         sigma(ntens+1:2*ntens,iElem) = sigma(1:ntens,iElem);
         eps(ntens+1:2*ntens,iElem) = eps(1:ntens,iElem);
         eps_fluct(ntens+1:2*ntens,iElem) = eps_fluct(1:ntens,iElem);
         hvar(sihvarpg+1:2*sihvarpg,iElem) = hvar(1:sihvarpg,iElem);
      end
      %
      % Operaciones que ocurren desde que se impone el SDA.
      %if condBif==2

      %end
   end

end