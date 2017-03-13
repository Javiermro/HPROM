function c_NormalesMicro = f_NormalesMicro(m_ElemLoc,e_VarEst,e_DatSet,m_NormalMacroVoigt,e_VG)

   %Leer las normales y longitudes micro impuestas en f_VarMECohesivo y imponerlas al  

   nSet = e_VG.nSet;
   ntens = e_VG.ntens;
   nDime = e_VG.ndime;
   
   % Normales geométrica
   %Nombre de archivo de normales geométricas (en realidad es el vector perpendicular a las normales
   %geométricas).
   s_Arch = [e_VG.fileCompleto,'.NormalMicro'];
   fId = fopen(s_Arch);
   m_NormalGeom = fscanf(fId,'%d %f %f %f',[4,Inf])';
   fclose(fId);
   m_NormalGeom(:,2:3) = bsxfun(@rdivide,m_NormalGeom(:,2:3),sqrt(sum(m_NormalGeom(:,2:3).^2,2)));
   m_NormalGeom(:,2:3) = m_NormalGeom(:,3:-1:2);
   m_NormalGeom(:,3) = -m_NormalGeom(:,3);
   
   % Normal macro
   m_NormalMacro = [m_NormalMacroVoigt(1,1);m_NormalMacroVoigt(2,2)];
   
   %Se guarda en la segunda columna de c_NormalesMicro las longitudes micros. Se guarda un valor por elemento
   %finito, y se aplica la misma a todos los puntos de gauss.
   c_NormalesMicro = cell(nSet,2);
   for iSet = 1:nSet
      
      nElem = e_DatSet(iSet).nElem;
      e_DatMat = e_DatSet(iSet).e_DatMat;
      e_DatElem = e_DatSet(iSet).e_DatElem;
      esImplex = e_DatMat.esImplex;      
      nPG = e_DatElem.npg;
      conshyp = e_DatMat.conshyp;
      
      hvar = reshape(e_VarEst(iSet).hvar,[],nPG,nElem);
      
      m_ElemLocSet = m_ElemLoc(e_DatSet(iSet).m_IndElemSet);
      
      %m_NormMicro = zeros(ntens,nDime,nPG,nElem);
      m_NormMicro = zeros(ntens,nDime,nElem);
      %m_LongMicro = zeros(1,nElem);
      %Se utiliza Inf para que en la operaciones de inserción de deformación y homegeneización de las
      %tracciones no resulte en NaN si su usa cero o NaN. Esto ocurre aún en el caso de que el lmicro
      %corresponda a un elemento finito que no pertenece al dominio localizado (debido la deformación fuera
      %del dominio que queda NaN, al ser 0/NaN ó 0/0, y no cero como quedaba antes).
      m_LongMicro = Inf(1,nElem);
      for iElem = 1:nElem
         if m_ElemLocSet(iElem)
            switch conshyp
               %case 1
                  %Si el material elástico se adopta la normal macro, pero capaz hay que adoptar la dirección
                  %de la tensión principal 1.
                  %m_NormMicro(:,:,:,iElem) = repmat(m_NormalMacroVoigt,[1,1,nPG]);
               case {1,10,12}
                  iElemNum = e_DatSet(iSet).m_NumElem(iElem);
                  m_indElemSet = m_NormalGeom(:,1)==iElemNum;
                  if isempty(m_indElemSet)
                     error('Normales micro: La normal y longitud micro del elemento %d no está definida',...
                        iElemNum)
                  end
                  m_NormalGeomElem = m_NormalGeom(m_indElemSet,2:3);
                  m_LongMicroElem = m_NormalGeom(m_indElemSet,4);                     
%                   if esImplex
%                      m_AngCrit = hvar(4:5,:,:);
%                   else
%                      m_AngCrit = hvar(3:4,:,:);
%                   end
                  %for iPG = 1:nPG                     
%                      m_NormalCrit = [cos(m_AngCrit(:,iPG,iElem)');sin(m_AngCrit(:,iPG,iElem)')];
%                      [~,ind] = min(abs(m_NormalGeomElem*m_NormalCrit));
%                      m_NormalCrit = m_NormalCrit(:,ind);
                     m_NormalCrit = m_NormalGeomElem';
                     %Se fuerza que la normal micro apunte a la normal macro.
                     if m_NormalCrit'*m_NormalMacro<0
                        m_NormalCrit = -m_NormalCrit;
                     end
                     %m_NormMicro(:,:,iPG,iElem) = f_Vec2Voigt2D(m_NormalCrit,e_VG);
                     m_NormMicro(:,:,iElem) = f_Vec2Voigt2D(m_NormalCrit,e_VG);
                     %
                     m_LongMicro(iElem) = m_LongMicroElem;
                  %end
               otherwise
                  error('Normales micro: Modelo constitutivo no definido.')
            end
         %else
            %En el caso que no pertenezca al dominio localizado no interesa la normal micro, por lo que se la
            %impone igual a cero.
            %No es necesario si las matrices están inicializados en cero.
            %m_NormMicro(:,:,:,iElem) = zeros(ntens,nDime,nPG);
            %m_NormMicro(:,:,iElem) = zeros(ntens,nDime);
         end
      end
      c_NormalesMicro{iSet,1} = m_NormMicro;
      c_NormalesMicro{iSet,2} = m_LongMicro;
      
   end
   
end