function [vectVHElem_new]  = ...
   Regularization_parameters_MICRO(xx_MICRO,e_VarEst_new,...
   e_VarEst_old,e_VarAux,delta_uMICRO,e_VG,e_DatSet,vectVHElem_old,...
   kinf,omegaMicro,ind_ActState_old)

% Calculo de la distancia S_mu en la microcelda
% *********************************************
length_fractMICRO_new = 0;
inc_dissip_CELL       = 0;
m_valHomog_old_5      = 0;
m_valHomog_new_5      = 0;
m_valHomog_old        = 0;
m_valHomog_new        = 0;
m_IncrRHomog = 0;

damage_micro          = 0;
area_damage           = 0;
areaElemDano = 0;

nSet = e_VG.nSet;
vectVHElem_new = vectVHElem_old;
for iSet = 1:nSet
   eltype            = e_DatSet(iSet).e_DatElem.eltype;
   switch eltype
      case 31 %sï¿½lo hacer el calculo para materiales de la banda de falla
         nElem_MICRO       = e_DatSet(iSet).nElem;
         m_DofElem         = e_DatSet(iSet).m_DofElem;
         e_DatElemSet      = e_DatSet(iSet).e_DatElem;
         conec_MICRO       = e_DatSet(iSet).conec;
         wg                = e_DatSet(iSet).e_DatElem.wg;
         npg               = e_DatSet(iSet).e_DatElem.npg;
         m_DetJT           = e_DatSet(iSet).m_DetJT;
         hvar_old          = e_VarEst_old(iSet).hvar;
         hvar_new          = e_VarEst_new(iSet).hvar;
         dN_xy             = e_DatSet(iSet).dN_xy  ;
         m_n_vector        = zeros(nElem_MICRO,e_VG.ndn);
         area_elem         = e_DatSet(iSet).m_VolElem;
         p_IntDissip       =  e_DatElemSet.pointersVHE.i_IntDissip  ;
         m_VarHistElemNew  =  e_VarEst_new(iSet).VarHistElem ;
         for iElem_MICRO = 1:nElem_MICRO
            hvar_OLD        = reshape(hvar_old(:,iElem_MICRO),[],npg);
            hvar_NEW        = reshape(hvar_new(:,iElem_MICRO),[],npg);
            m_pesoPG        = m_DetJT(:,iElem_MICRO).*wg;
            %
            m_valHomog_old_5  = m_valHomog_old_5+hvar_OLD(5,:)*m_pesoPG;
            m_valHomog_new_5  = m_valHomog_new_5+hvar_NEW(5,:)*m_pesoPG;
            m_valHomog_old  = m_valHomog_old+ hvar_OLD(7,:)*m_pesoPG;
            m_valHomog_new  = m_valHomog_new+ hvar_NEW(7,:)*m_pesoPG;
            damage_micro = damage_micro+hvar_NEW(1,:)*m_pesoPG;
            delta_hvar_test = hvar_NEW(4,:)-hvar_OLD(4,:);
            % material que dania y ablanda con deformacion - ielem_material==1 la segunda condicion revisa que en los 4 PGS se este presentando
            % carga, de lo contrario no tiene en cuenta el elemento finito para el calculo de la longitud de fallo
            m_IntDissip     =  m_VarHistElemNew (p_IntDissip   ,iElem_MICRO)   ;
            inc_dissip_CELL = inc_dissip_CELL + m_IntDissip;
            %% if (iSet==1) && (sum(delta_hvar_test(11:data_micro.e_VG_MICRO.sihvarpg:end,iElem_MICRO)>0)==4)
            if (sum(delta_hvar_test>0)==npg)
               %%%       dof           = f_DofElem(conec_MICRO(iElem_MICRO,:),ndn);
               dof           = m_DofElem(:,iElem_MICRO);
               delta_u       = delta_uMICRO(dof);
               conec_Elem    = conec_MICRO(iElem_MICRO,:);
               m_ConecFrontElem = [conec_Elem(1) conec_Elem(2);conec_Elem(2) conec_Elem(3);...
                  conec_Elem(3) conec_Elem(4);conec_Elem(4) conec_Elem(1)];
               m_n_vector(iElem_MICRO,:) = matrix_grad_u_quad_q1_epd(delta_u,dN_xy(:,:,iElem_MICRO));
               tvector_Elem = xx_MICRO(m_ConecFrontElem(:,2),:)-xx_MICRO(m_ConecFrontElem(:,1),:);
               le_Elem = sqrt(sum(tvector_Elem.^2,2));
               tvector_Elem = bsxfun(@rdivide,tvector_Elem,le_Elem);
               n_side_Elem = [tvector_Elem(:,2),-tvector_Elem(:,1)];
               dot_pr = abs(n_side_Elem*m_n_vector(iElem_MICRO,:)');
               for i = 1:4
                  if dot_pr(i)<0.997
                     le_Elem(i)=0;
                  end
               end
               % LONGITUD TOTAL DE FISURACION EN LA MICROESTRUCTURA
               length_fractMICRO_new = length_fractMICRO_new +  sum(le_Elem)/2;
               % AREA TOTAL DE LA BANDA FISURADA ACTIVA EN LA MICROESTRUCTURA
               %  area_FractMICRO_new = area_FractMICRO_new + m_VolElem(iElem_MICRO);
               % INFORMACION DE VARIABLES INTERNAS EN LOS SEGMENTOS DE BANDAS ACTIVOS
               % hvar_assessment_q(:,iElem_MICRO) = hvar_new(:,iElem_MICRO);
            end
            area_damage = area_damage + area_elem(iElem_MICRO);
         end
      case 32
         conshyp  = e_DatSet(iSet).e_DatMat.conshyp;
         switch  conshyp
            case {11,12}
               
               nElem_MICRO       = e_DatSet(iSet).nElem;
               e_DatElemSet      = e_DatSet(iSet).e_DatElem;
               wg                = e_DatSet(iSet).e_DatElem.wg;
               npg               = e_DatSet(iSet).e_DatElem.npg;
               m_DetJT           = e_DatSet(iSet).m_DetJT;
               hvar_old          = e_VarEst_old(iSet).hvar;
               hvar_new          = e_VarEst_new(iSet).hvar;
               area_elem         = e_DatSet(iSet).m_VolElem;
               p_IntDissip       =  e_DatElemSet.pointersVHE.i_IntDissip  ;
               m_VarHistElemNew  =  e_VarEst_new(iSet).VarHistElem ;
               %le_Elem = e_DatSet(iSet).e_DatElem.le_Elem;
               for iElem_MICRO = 1:nElem_MICRO
                  hvar_OLD        = reshape(hvar_old(:,iElem_MICRO),[] , npg);
                  hvar_NEW        = reshape(hvar_new(:,iElem_MICRO),[] , npg);
                  m_pesoPG        = m_DetJT(:,iElem_MICRO).*wg;
                  %
                  %% Criterio de descarga del elemento
                  %Se utiliza el incremento implícito de r.
                  %vectVHElem_new(5) = hvar_new(3,6);
                  %m_valHomog_old_5  = m_valHomog_old_5+hvar_OLD(5,:)*m_pesoPG;
                  %m_valHomog_new_5  = m_valHomog_new_5+hvar_NEW(5,:)*m_pesoPG;
                  %Se utiliza el incremento explícito de r.
                  m_IncrRHomog = m_IncrRHomog+hvar_NEW(8,:)*m_pesoPG;
                  %
                  m_valHomog_old  = m_valHomog_old+ hvar_OLD(7,:)*m_pesoPG;
                  m_valHomog_new  = m_valHomog_new+ hvar_NEW(7,:)*m_pesoPG;
                  damage_micro = damage_micro +  hvar_NEW(1,:)*m_pesoPG;
                  %delta_hvar_test = hvar_NEW(4,:)-hvar_OLD(4,:);
                  % material que dania y ablanda con deformacion - ielem_material==1 la segunda condicion revisa que en los 4 PGS se este presentando
                  % carga, de lo contrario no tiene en cuenta el elemento finito para el calculo de la longitud de fallo
                  %Se realiza una media de la disipación solo en los elementos micro que dañaron.
                  %Se empieza a sumar a partir que bifurca (sólo durante el estado 1).
                  if any(hvar_NEW(8,:)>0)&&ind_ActState_old>0
                     areaElemDano = areaElemDano+area_elem(iElem_MICRO);
                     m_IntDissip = m_VarHistElemNew(p_IntDissip,iElem_MICRO);
                     inc_dissip_CELL = inc_dissip_CELL+m_IntDissip;
                  end
                  %Esto parece que se utiliza únicamente con el modelo Barcelona.
%                   if (sum(delta_hvar_test>0)==npg)% LONGITUD TOTAL DE FISURACION EN LA MICROESTRUCTURA
%                      length_fractMICRO_new = length_fractMICRO_new +  sum(le_Elem(iElem_MICRO));
%                   end
                  area_damage = area_damage + area_elem(iElem_MICRO);
               end
         end
   end
end

%m_valHomog_old_5  = m_valHomog_old_5 /omegaMicro ;
%m_valHomog_new_5  = m_valHomog_new_5 /omegaMicro;
m_valHomog_old    = m_valHomog_old /omegaMicro ;
m_valHomog_new    = m_valHomog_new /omegaMicro;
%vectVHElem_new(7) = m_valHomog_new_5-m_valHomog_old_5; % CRITERIO PARA RECUPERAR EL ELEMENTO FINITO ESTANDAR EN MODELOS MULTIESCALA.
%% Criterio de descarga homogeneizado
%%vectVHElem_new(5) = m_valHomog_new_5-m_valHomog_old_5; 
vectVHElem_new(5) = m_IncrRHomog/omegaMicro;
%fprintf('Criterio de carga %f del elemento %d.\n',vectVHElem_new(5),e_VG.iElemNumMacro)
%
%vectVHElem_new(5) = vectVHElem_old(5)+ kinf*(m_valHomog_new-m_valHomog_old); % VARIABLE DE SUAVIZADO MULTIESCALA
vectVHElem_new(4) = vectVHElem_old(4)+ kinf*(m_valHomog_new-m_valHomog_old); % VARIABLE DE SUAVIZADO MULTIESCALA

% if length_fractMICRO_new==0
%    length_fractMICRO_new=1e-16;
% end
% vectVHElem_new(7) = length_fractMICRO_new ;     % S_mu
% vectVHElem_new(6) = omegaMicro/length_fractMICRO_new; %long caracteristica del modelo multiescala
% 
vectVHElem_new(8) = damage_micro/area_damage;

%if area_FractMICRO_new==0; area_FractMICRO_new=1e-16; end

% Es necesario almacenar el valor de la longitud de fallo de la microcelda
%VPG_new.length_fractMICRO  = length_fractMICRO_new ;
%VPG_new.char_length_MICRO = data_micro.area/length_fractMICRO_new;


%m_DsigmaHomog     = f_HomogArea(c_Tens_new,ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VG);

% NEW ENERGY CRITERIA
% % *******************
% for iElem_MICRO = 1:nElem_MICRO
%     ielem_material = data_micro.e_VG_MICRO.iel(iElem_MICRO);
%     % material que dania y ablanda con deformacion - ielem_material==1
%     if (ielem_material==1)
%         inc_dissip_CELL = inc_dissip_CELL + incr_diss(iElem_MICRO);
%     end
% end

%if (length_fractMICRO_new==1e-16); FACT=0; else FACT=(kinf/data_micro.area); end
%if (length_fractMICRO_new==1e-16)||areaElemDano==0
if areaElemDano==0 
   FACT=0;
else
   %FACT=(kinf/omegaMicro);
   FACT=(kinf/areaElemDano);
end
% Variable de disipacion acumulada en la microcelda
%alf VPG_new.dissip_CELL  = VPG_old.dissip_CELL + FACT*inc_dissip_CELL;
vectVHElem_new(1)  = vectVHElem_old(1) + FACT*inc_dissip_CELL;
%vectVHElem_new(1) = vectVHElem_old(1)+inc_dissip_CELL/omegaMicro;


%alf % *******************************************************
%alf % ACTUALIZACION VARIABLES HOMOGENIZADAS - MICROESTRUCTURA
%alf % *******************************************************
%alf hvar_old_test = hvar_old;
%alf Hbar = (data_micro.Eprop_MICRO(1,6)^2)/(2*data_micro.Eprop_MICRO(1,4)*data_micro.Eprop_MICRO(1,14));
%alf
%alf for j = 1:nElem_MICRO
%alf     if sum(hvar_old_test(5:data_micro.e_VG_MICRO.sihvarpg:end,j)==0)==4 && ...
%alf             sum(hvar_old_test(6:data_micro.e_VG_MICRO.sihvarpg:end,j)==0)==4
%alf         %
%alf         hvar_old_test(5:data_micro.e_VG_MICRO.sihvarpg:end,j) = ...
%alf             hvar_new(5:data_micro.e_VG_MICRO.sihvarpg:end,j); % Problema de asignacion en la primera iteracion de r
%alf        %
%alf        hvar_old_test(6:data_micro.e_VG_MICRO.sihvarpg:end,j) = ...
%alf             hvar_new(6:data_micro.e_VG_MICRO.sihvarpg:end,j); % Problema de asignacion en la primera iteracion de q
%alf     %
%alf     end
%alf end

%m_hvar_Homog_old = f_HomogArea(hvar_old,data_micro.area,m_DetJT,data_micro.e_VG_MICRO);
%m_hvar_Homog     = f_HomogArea(hvar_new,data_micro.area,m_DetJT,data_micro.e_VG_MICRO);
%%m_hvar_Homog(10,1) = m_hvar_Homog(11,1)-m_hvar_Homog_old(11,1); % ro_new - ro_old (delta_ro) para inyectar la Strong Discontinuity

% delta_hvar = hvar_new - hvar_old_test;
% m_delta_hvar_Homog = f_HomogArea(delta_hvar,data_micro.area,m_DetJT,data_micro.e_VG_MICRO);

%alfredo   m_delta_hvar_Homog_varq = f_HomogArea(delta_hvar,area_FractMICRO_new,m_DetJT,data_micro.e_VG_MICRO);
%if length_fractMICRO_new==1e-16; m_delta_hvar_Homog(5,1)=0; end
%m_delta_hvar_Homog_varq2 = f_HomogArea(hvar_assessment_q,area_FractMICRO_new,m_DetJT,data_micro.e_VG_MICRO);




% variable 11 no es necesaria
% m_hvar_Homog(11,1) = hvar_Homog_old_MACRO(11,1) + VPG_new.char_length_MICRO*m_delta_hvar_Homog(5,1);

% VARIABLE DE INYECCION
%%%m_hvar_Homog(15,1) = hvar_Homog_old_MACRO(15,1) + Hbar*VPG_new.char_length_MICRO*m_delta_hvar_Homog(5,1); %HOMOGENIZED Q CRITERIA % 30/05/2013
%m_hvar_Homog(15,1) = VPG_new.dissip_CELL; % FRACTURE ENERGY CRITERIA

% VARIABLES AUXILIARES POSTPROCESADAS
%m_hvar_Homog(16,1) = norm(m_DefMACRO);
%%% m_hvar_Homog(17,1) = hvar_Homog_old_MACRO(17,1) + m_delta_hvar_Homog_varq(6,1);
%m_hvar_Homog(17,1) = m_delta_hvar_Homog_varq2(6,1);


