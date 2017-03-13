function [e_DatSet,e_VarEst_new,e_VarAux,e_VG] =...
   f_OperPosConv_CPF(u,xx,e_VarEst_new,e_VarAux,e_DatSet,e_VG)

%NORMAL VECTOR SMOOTHING  &  indicator field SMOOTHING
V = zeros(e_VG.nnod,1);
for iSet = 1:e_VG.nSet
   e_DatElemSet  = e_DatSet(iSet).e_DatElem;
   conshyp       = e_DatSet(iSet).e_DatMat.conshyp;
   eltype        = e_DatElemSet.eltype;
   %       npg           = e_DatElemSet.npg;
   conec         = e_DatSet(iSet).conec;
   dofElemSet    = e_DatSet(iSet).m_DofElem(:);
   nElem         = e_DatSet(iSet).nElem;
   m_IndElemSet = e_DatSet(iSet).m_IndElemSet;
   uElemSet      = reshape(u(dofElemSet),[],nElem);
   m_VarAuxElem  = e_VarAux(iSet).VarAuxElem;
   %       sihvarpg      = e_DatSet(iSet).e_DatMat.sihvarpg;
   %       hvar_new      = e_VarEst_new(iSet).hvar;
   m_VarHistElemNew = e_VarEst_new(iSet).VarHistElem;
   switch eltype
      case {21,22,23}
         
         switch conshyp
            case {11,53,54}               
               i_vectVHElem   =  e_DatElemSet.pointersVHE.i_vectVHElem;
               p_varHistSmooth= i_vectVHElem(4);
               
               p_nSmoothing  = e_DatElemSet.pointersVAE.p_nSmoothing;
               dN_xy     =  e_DatSet(iSet).dN_xy  ;
               MElemInv  =  e_DatSet(iSet).MElemInv;
               NT_MACRO  =  e_DatSet(iSet).N_vector;
               for iElem =  1:nElem
                  %Ver si no sacar este condicional afuera del loop de elementos.
                  if ~isfield(e_VG,'angSmoothingImp')&&~isfield(e_VG,'m_AngSmoothImpEl')
                     % Selection of normal vector for smoothing
                     nSmoothing = ...
                        Smoothing_normal_vector...
                        (xx,conec(iElem,:),uElemSet(:,iElem),dN_xy(:,:,iElem),eltype,e_VG);
                     %ang = -70*pi/180;
                     %ang = 0;
                     %nSmoothing=[cos(ang),0;0,sin(ang);0,0;sin(ang),cos(ang)];
                  elseif isfield(e_VG,'angSmoothingImp')
                     ang = e_VG.angSmoothingImp;
                     nSmoothing = [cos(ang),0;0,sin(ang);0,0;sin(ang),cos(ang)];
                  else
                     ang = e_VG.m_AngSmoothImpEl(m_IndElemSet(iElem));
                     nSmoothing = [cos(ang),0;0,sin(ang);0,0;sin(ang),cos(ang)];
                  end
                  m_VarAuxElem(p_nSmoothing,iElem) = nSmoothing(:);
                  % indicator field SMOOTHING
                  HV_smooth = m_VarHistElemNew(p_varHistSmooth,iElem);
                  V_alpha = MElemInv(:,:,iElem)*NT_MACRO(:,iElem).*HV_smooth;
                  %                     for j = 1:length(conec(iElem,:))
                  %                         V(conec(iElem,j)) = V(conec(iElem,j)) + V_alpha(j);
                  %                     end
                  V(conec(iElem,:)) = V(conec(iElem,:)) + V_alpha(:);
               end
               e_VarAux(iSet).VarAuxElem = m_VarAuxElem;
            otherwise
               error(['Operaciones de fin de paso: Normal de Smoothing: Cuadrángulo SDA bilineal ',...
                  'del modelo multiescala cohesivo: Modelo constitutivo no definido.'])
         end
         
      case 4
         p_nSmoothing  = e_DatElemSet.pointersVAE.p_nSmoothing ;
         dN_xy         =  e_DatSet(iSet).dN_xy  ;
         for iElem =  1:nElem
            % Selection of normal vector for smoothing
            %Ver si no sacar este condicional afuera del loop de elementos.
            if ~isfield(e_VG,'angSmoothingImp')
               nSmoothing = ...
                  Smoothing_normal_vector...
                  (xx,conec(iElem,:),uElemSet(:,iElem),dN_xy(:,:,iElem),eltype,e_VG);       
               %ang = -70*pi/180;
               %ang = 0;
               %nSmoothing=[cos(ang),0;0,sin(ang);0,0;sin(ang),cos(ang)];           
            elseif isfield(e_VG,'angSmoothingImp')
               ang = e_VG.angSmoothingImp;
               nSmoothing = [cos(ang),0;0,sin(ang);0,0;sin(ang),cos(ang)];
            else
               ang = e_VG.m_AngSmoothImpEl(m_IndElemSet(iElem));
               nSmoothing = [cos(ang),0;0,sin(ang);0,0;sin(ang),cos(ang)];
            end
            %
            m_VarAuxElem(p_nSmoothing ,iElem) = nSmoothing(:);
         end
         e_VarAux(iSet).VarAuxElem = m_VarAuxElem;
      otherwise
         error('Operaciones de fin de paso: Normal de Smoothing: Elemento Finito no definido.')
   end
end
smooth_alpha = e_VG.MGlobInv.*V;
e_VG.smooth_alpha = smooth_alpha;
%
% SMOOTHING of the norma derivative of the indicator field
V = zeros(e_VG.nnod,1);
for iSet = 1:e_VG.nSet
   e_DatElemSet  = e_DatSet(iSet).e_DatElem;
   eltype        = e_DatElemSet.eltype;
   conec         = e_DatSet(iSet).conec;
   nElem         = e_DatSet(iSet).nElem;
   m_VarAuxElem  = e_VarAux(iSet).VarAuxElem;
   switch eltype
      case {21,22,23}
         p_nSmoothing  = e_DatElemSet.pointersVAE.p_nSmoothing ;
         M_elcorr_MASA_MACRO  =  e_DatSet(iSet).MElemInv;
         NT_MACRO             =  e_DatSet(iSet).N_vector;
         dN_xy                =  e_DatSet(iSet).dN_xy  ;
         for iElem = 1:nElem
            nSmoothing = reshape( m_VarAuxElem(p_nSmoothing ,iElem) , 4 ,2) ;
            d_alpha = [nSmoothing(1,1) nSmoothing(2,2)]*(dN_xy(:,:,iElem)*smooth_alpha(conec(iElem,:)));
            V_alpha = M_elcorr_MASA_MACRO(:,:,iElem)*NT_MACRO(:,iElem).*d_alpha;
            for j = 1:length(conec(iElem,:))
               V(conec(iElem,j)) = V(conec(iElem,j)) + V_alpha(j);
            end
         end
      case 4
         p_nSmoothing  = e_DatElemSet.pointersVAE.p_nSmoothing ;
         M_elcorr_MASA_MACRO  =  e_DatSet(iSet).MElemInv;
         NT_MACRO             =  e_DatSet(iSet).N_vector;
         dN_xy                =  e_DatSet(iSet).dN_xy  ;
         for iElem = 1:nElem
            nSmoothing = reshape( m_VarAuxElem(p_nSmoothing ,iElem) , 4 ,2) ;
            d_alpha = [nSmoothing(1,1) nSmoothing(2,2)]*(dN_xy(:,:,iElem)*smooth_alpha(conec(iElem,:)));
            V_alpha = M_elcorr_MASA_MACRO(:,:,iElem)*NT_MACRO(:,iElem).*d_alpha;
            for j = 1:length(conec(iElem,:))
               V(conec(iElem,j)) = V(conec(iElem,j)) + V_alpha(j);
            end
         end
   end
end
smooth_dalpha = e_VG.MGlobInv.*V;
e_VG.smooth_dalpha = smooth_dalpha;

%---------------------------------------
%Freezing elements belong the crack path
%---------------------------------------
FREEZING_FACT = 0.85;
for iSet = 1:e_VG.nSet
   e_DatElemSet     = e_DatSet(iSet).e_DatElem;
   conshyp          = e_DatSet(iSet).e_DatMat.conshyp;
   eltype           = e_DatElemSet.eltype;
   conec            = e_DatSet(iSet).conec;
   m_VarHistElemNew = e_VarEst_new(iSet).VarHistElem;
   switch eltype
      case {21,22,23}
         FRAC_ENER = e_DatElemSet.gfvRef;
%          switch conshyp
%             case 11
%                FRAC_ENER = e_DatSet(iSet).e_DatMat.gfv;    % set = logical(e_VG.iel==iSET);
%             case {53,54}
%                %aaa=1;
%                %%%   Set= 1;  % modificar esto !!!!!
%                Set= 3;  % modificar esto !!!!!
%                FRAC_ENER = e_DatSet(iSet).e_DatMat.e_DatSet(Set).e_DatMat.gfv;
%                %                    FRAC_ENER = 0;
%          end
         
         i_vectVHElem   =  e_DatElemSet.pointersVHE.i_vectVHElem ;
         p_varDissip= i_vectVHElem(1);
         elem_freeze = find(m_VarHistElemNew(p_varDissip,:)>FREEZING_FACT*FRAC_ENER);
         if find(elem_freeze>0)
            for iElem_freeze = 1:length(elem_freeze)
               e_VG.smooth_dalpha(conec(elem_freeze(iElem_freeze),:)) = ...
                  e_VG.smooth_dalpha_old(conec(elem_freeze(iElem_freeze),:));
            end
            %else
            %    aaaa=1;
         end
   end
end


e_VG.smooth_dalpha_old= e_VG.smooth_dalpha;

% smooth_alpha = M_INV_MACRO*V;

% Puntos de corte de los vertices de los elementos finitos
% **** Para la seleccion de la normal de referencia ******
for iSet = 1:e_VG.nSet
   e_DatElemSet   = e_DatSet(iSet).e_DatElem;
   eltype         = e_DatElemSet.eltype;
   nElem         = e_DatSet(iSet).nElem;
   switch eltype
      case {21,22,23}
         % EF_sidesCutCPF = e_DatSet(iSet).e_DatElem.EF_sidesCutCPF;
         EF_sidesCutCPF = zeros(nElem,1);
         NumsidesCutCPF = zeros(nElem,1);
         conec          = e_DatSet(iSet).conec;
         for iElem = 1:nElem
            % Seleccion de los lados por donde pasa la discontinuidad
            if (smooth_dalpha(conec(iElem,1))*smooth_dalpha(conec(iElem,2))<0) && (smooth_dalpha(conec(iElem,2))*smooth_dalpha(conec(iElem,3))<0)
               EF_sidesCutCPF(iElem) = 12;
            elseif (smooth_dalpha(conec(iElem,1))*smooth_dalpha(conec(iElem,2))<0) && (smooth_dalpha(conec(iElem,3))*smooth_dalpha(conec(iElem,4))<0)
               EF_sidesCutCPF(iElem) = 13;
            elseif (smooth_dalpha(conec(iElem,1))*smooth_dalpha(conec(iElem,2))<0) && (smooth_dalpha(conec(iElem,4))*smooth_dalpha(conec(iElem,1))<0)
               EF_sidesCutCPF(iElem) = 14;
            elseif (smooth_dalpha(conec(iElem,2))*smooth_dalpha(conec(iElem,3))<0) && (smooth_dalpha(conec(iElem,3))*smooth_dalpha(conec(iElem,4))<0)
               EF_sidesCutCPF(iElem) = 23;
            elseif (smooth_dalpha(conec(iElem,2))*smooth_dalpha(conec(iElem,3))<0) && (smooth_dalpha(conec(iElem,4))*smooth_dalpha(conec(iElem,1))<0)
               EF_sidesCutCPF(iElem) = 24;
            elseif (smooth_dalpha(conec(iElem,3))*smooth_dalpha(conec(iElem,4))<0) && (smooth_dalpha(conec(iElem,4))*smooth_dalpha(conec(iElem,1))<0)
               EF_sidesCutCPF(iElem) = 34;
               %                     else
               %                         EF_sidesCutCPF(iElem) = 0;
            end
            
            cont=0;
            if (smooth_dalpha(conec(iElem,1))*smooth_dalpha(conec(iElem,2))<0)
               cont=cont+1;
            end
            if (smooth_dalpha(conec(iElem,2))*smooth_dalpha(conec(iElem,3))<0)
               cont=cont+1;
            end
            if (smooth_dalpha(conec(iElem,3))*smooth_dalpha(conec(iElem,4))<0)
               cont=cont+1;
            end
            if (smooth_dalpha(conec(iElem,4))*smooth_dalpha(conec(iElem,1))<0)
               cont=cont+1;
            end
            
            NumsidesCutCPF(iElem) = cont;
            if cont > 2
               EF_sidesCutCPF(iElem) = 0;
            end
         end
         e_DatSet(iSet).e_DatElem.EF_sidesCutCPF = EF_sidesCutCPF;
         e_DatSet(iSet).e_DatElem.NumsidesCutCPF = NumsidesCutCPF;
   end
end
end
