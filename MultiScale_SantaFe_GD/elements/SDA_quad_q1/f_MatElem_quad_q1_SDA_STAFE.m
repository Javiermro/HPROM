function [kuu,kbu_SDA,kbu_MSL,kbbIN_SDA,kbbIN_MSL,fint,res_beta_SDA,...
    res_beta_MSL,sigma_new,hvar_new,eps_new,m_TensorTang,...
    ind_state_new,ind_ActState_new,elem_type_new,ksd,vectVHElem_new,sigmaTilde_new] =...
 f_MatElem_quad_q1_SDA_STAFE...
     (delta_u,eps_old,Dbeta_e_SDA,Dbeta_e_MSL,...
     m_VarAuxPG,BI_n,leq_elem,Phi_Grad,n_tensor,sigma_old,...
     hvar_old,e_DatElemSet, e_DatMatSet,m_Be,m_DetJe,ind_state_old,ind_ActState_old,...
     e_VG,vectVHElem_old,fact_inyect,CPI_n,sigmaTilde_old)
  
 TWO_SCALE      = 0 ;  % comentado para la unificacion de codigos aeh ago_2013
 
 ntens = e_VG.ntens;
 ndime = e_VG.ndime;
 dofpe = e_DatElemSet.dofpe;
 npg = e_DatElemSet.npg;
 wg = e_DatElemSet.wg;
 
 % Propiedades materiales
 sihvarpg = e_DatMatSet.sihvarpg;
 siavarpg = e_DatMatSet.siavarpg;
 conshyp  = e_DatMatSet.conshyp;
 esImplex = e_DatMatSet.esImplex;
 
 if conshyp  == 54
     TWO_SCALE      = 1 ;
 end
 
 % Inicializaciones
 % ****************
 kuu = zeros(dofpe,dofpe);
 fint = zeros(dofpe,1);
 
 % SDA ELEMENT
 kbu_SDA      = zeros(ndime,dofpe);
 kbbIN_SDA    = zeros(ndime,ndime);
 res_beta_SDA = zeros(ndime,1);
 
 % MSL ELEMENT
 kbu_MSL      = zeros(ntens,dofpe);
 kbbIN_MSL    = zeros(ntens,ntens);
 res_beta_MSL = zeros(ntens,1);
 
 sigma_new     = zeros(ntens,npg);                     % Included Gauss Points (5 & 6) in the center of the finite element
 sigmaTilde_new = zeros(ntens,npg);
 eps_new       = zeros(ntens,npg);                       % Included Gauss Points (5 & 6) in the center of the finite element
 Btotal        = zeros(ntens,dofpe);
 
 %Para la homogenizacion del tensor tangente
 if esImplex
     m_TensorTang = zeros(ntens,ntens,2*npg);
 else
     m_TensorTang = zeros(ntens,ntens,npg);
 end
 
 % Redimensionado de matrices
 sigma_old        = reshape(sigma_old,ntens,[]);
 sigmaTilde_old   = reshape(sigmaTilde_old,ntens,[]);
 eps_old          = reshape(eps_old,ntens,[]);
 hvar_old         = reshape(hvar_old,sihvarpg,[]);
 aux_var          = reshape(m_VarAuxPG,siavarpg,npg);
 hvar_new         = f_DefinicionhVar(conshyp,sihvarpg,npg);
 vectVHElem_new   = vectVHElem_old;
 
 ind_state_new    = 0;
 elem_type_new    = 0;
 ind_ActState_new = 0;
 
 % SDA PARAMETERS: definicion del ancho de banda
 % *********************************************
 area      = m_DetJe(1:4)'*wg(1:4);
 h_eq_zero = sqrt(area);
 
 % Problema cuando no se tiene l_equivalente en el elemento finito
 
 h_eq = leq_elem; if h_eq ==0; h_eq = h_eq_zero; end
 
 ksd = h_eq;
 
 if TWO_SCALE
     kinf = e_DatElemSet.kinf;
     le   = e_DatElemSet.le ;
 else
     kinf = e_DatElemSet.kinf*h_eq ;
     le   = e_DatElemSet.le*h_eq  ;
 end

 % Parametros iniciales para el elemento finito.
 gamma = 0;
 tau = 1;
  
 elem_type_new = 0;
 
 SLI_n = vectVHElem_old(5);
 % RLI_n = vectVHElem_old(8);
  
 if (BI_n == 1)

     if conshyp==2
     elseif (conshyp==4 || conshyp==44 || conshyp==11 || conshyp==53|| conshyp==54)
         ro     = vectVHElem_old(1) ;  %  hvar_old(15,6);
         ro_inj = vectVHElem_old(2) ;  %   hvar_old(8,6);
         ro_sda = vectVHElem_old(3) ;  %  hvar_old(9,6);
     end     
     
     % FROM (STATE==0) TO (STATE==1 || STATE==2)
     if ind_ActState_old == 0
         % CASE 1
         %if ((RLI_n>0) && (ro<ro_inj)) || (((RLI_n>0) && (ro>ro_inj)) && (CPI_n==0 || CPI_n==4))
         if ((SLI_n>0) && (ro<ro_inj)) || (((SLI_n>0) && (ro>ro_inj)) && (CPI_n==0 || CPI_n==4))    
             gamma = 1;  tau = 0.;
             ind_state_new = 1;
             ind_ActState_new = 1 ;
             ksd = kinf;
             
             % CASE 2
         %elseif (((SLI_n>0) && (ro>ro_inj)) || ((RLI_n>0) && (ro>ro_inj))) && (CPI_n==2)
         elseif (((SLI_n>0) && (ro>ro_inj)) || ((SLI_n>0) && (ro>ro_inj))) && (CPI_n==2)
             gamma = 1;  tau = 0.;
             ind_state_new = 2;
             ind_ActState_new = 2 ;
             ksd = kinf;
         end
     end
     
     % FROM (STATE==1) TO (STATE==0 || STATE==2)
     if ind_ActState_old == 1

         % CASE 3
         if (SLI_n<=0)
             gamma = 0;  tau = 1;
             ind_state_new = 1;
             ind_ActState_new = 0 ;
             ksd = 0.99999*h_eq;
             
             % CASE 4
         elseif (SLI_n>0) && (ro>ro_inj) && (CPI_n==2)
             gamma = 1;  tau = 0;
             ind_state_new = 2;
             ind_ActState_new = 2;
             ksd = kinf;
         
             
         else
             ind_ActState_new = ind_ActState_old;
             gamma = 1;  tau = 0;
             ind_state_new = 1;
             ksd = kinf;
             
         end
         
     end
     
     % FROM (STATE==2) TO (STATE==0 || STATE==1)
     if ind_ActState_old == 2
%          % CASE 5
% %alf_enero14
%          if (SLI_n<=0)
%              gamma = 0;  tau = 1;
%              ind_state_new = 2;
%              ind_ActState_new = 0 ;
%              ksd = 0.99999*h_eq;
%              % CASE 6
%          elseif (SLI_n>0) && (CPI_n==0 || CPI_n==4)
%              gamma = 1;  tau = 0;
%              ind_state_new = 2;             
%              ind_ActState_new = 1;
%              ksd = kinf;
%          else
%              ind_ActState_new = ind_ActState_old;
%              gamma = 1;  tau = 0;
%              ind_state_new = 2;
%              ksd = kinf;
%          end
%alf_enero14
%alf_enero14
              ind_ActState_new = ind_ActState_old;
              gamma = 1;  tau = 0;
              ind_state_new = 2;
              ksd = kinf;
%alf_enero14
     end    
 end
 
 if ind_state_new<ind_state_old
     ind_state_new = ind_state_old;
 end

 
C_inyect_SD=1;
if ind_ActState_new==2 && ind_ActState_old<=1
   C_inyect_SD=0;
   gamma = 1;  tau = 0.;
   ksd = kinf;
end

 if TWO_SCALE == 1
     if ksd>10; ksd=kinf; end % VERIFICAR PORQUE Leq SI PUEDE SER >10 DEPENDIENDO DEL GROSOR DEL MALLADO
     data_micro.Eprop_MICRO(:,15)=ksd; data_micro.Eprop_MICRO(:,17)=gamma;
 else
     if gamma == 0; ksd=0.99999*h_eq; end
     Eprop_qSD.ksd=ksd; Eprop_qSD.gamma=gamma;
 end
 
if ind_ActState_new==2 && ind_ActState_old==2 && C_inyect_SD
    kinf= 2.e-4; ksd=kinf;Eprop_qSD.ksd=ksd;
end

 
 
  m_pesoPG = m_DetJe.*wg;
 
  % B MATRIX ELEMENTAL INTEGRATION
  % ******************************
  for jPG = 1:(npg-2)
      Btotal = Btotal + m_Be(:,:,jPG)*m_pesoPG(jPG);
  end
  
  psi=ksd/h_eq;
  if gamma == 1 % PESOS DE GAUSS INTEGRACION MATRIZ K
      %m_pesoPG(1:4) = 0;
      m_pesoPG(5) = (1.0-ksd/h_eq)*area;
      m_pesoPG(6) = (ksd/h_eq)*area;

      if ind_ActState_new == 2 && C_inyect_SD  %%alfredo 
         m_pesoPG(5) = area;
         m_pesoPG(6) = 0;
         psi=0;
      end
  end
  
  % *************************************************************
  % *************************************************************
  % REGULAR GAUSS POINT (GP 5)
  % STRAIN & STRESS TENSORS ON THE 5th GAUSS POINT (UNLOAD POINT)
  % *************************************************************
  % *************************************************************
  B               = m_Be(:,:,5);
%alf14  m_DefMACRO      = eps_old(:,5) + B*delta_u;
  m_DefMACRO      =  B*delta_u;
  
  % ELEMENT SELECTION
  if ind_ActState_new == 1 || (ind_ActState_new == 2 && ~C_inyect_SD)  
      % MSL ELEMENT
      m_DefMACRO = B*delta_u - gamma*(1/(h_eq-ksd))*Dbeta_e_MSL;
  elseif ind_ActState_new == 2 && C_inyect_SD
      % SDA_ELEMENT
      m_DefMACRO = B*delta_u - gamma*Phi_Grad*Dbeta_e_SDA;
  end
  
  d_eps_new(:,5) = m_DefMACRO;
  eps_new(:,5)   = eps_old(:,5) + d_eps_new(:,5) ;

  switch conshyp % Modelo constitutivo
      case 11   % modelo da�o monoescala %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if gamma == 0
              [ctR,sigma_new(:,5),sigma_new_impl(:,5),hvar_new(:,5)] = ...
               rmap_damage_IMPLEX (eps_new(:,5),hvar_old(:,5),e_DatMatSet,...
                       Eprop_qSD);
          else
              [ctR,sigma_new(:,5),sigma_new_impl(:,5),hvar_new(:,5)] ...
                  = rmap_damage_elastic(eps_new(:,5),hvar_old(:,5),e_DatMatSet) ;
          end
      case 54   %Modelo multiescala cohesivo SANTA FE PG 5
          e_VG.iPG = 5;
          if gamma   == 0
              e_VG.elast = 0;
          else
              e_VG.elast = 1;
          end
          condBifR = 0;
          m_IDefR = d_eps_new(:,5);
          m_IDefS = zeros(ntens,1);
          %
         % [ctR,sigma_new(:,5),hvar_new(:,5),~] =...
         %     f_RMap_MEBcna(m_IDefR,hvar_old(:,5),...
         %     e_DatMatSet,condBifR,e_VG,vectVHElem_old,kinf);
          [ctR,sigma_new(:,5),hvar_new(:,5),~, ~] =...
              f_RMap_MEStafe(m_IDefR,m_IDefS,hvar_old(:,5),...
                  e_DatMatSet,condBifR,ind_ActState_new,C_inyect_SD,...
                  e_VG,vectVHElem_old,kinf);
  end
  c_tilde5              = ctR.Implex;
  m_TensorTang(:,:,5)   = c_tilde5 ;
  m_TensorTang(:,:,11)  = ctR.Impli;
  delta_sigma(:,5) = sigma_new(:,5)-sigma_old(:,5);
  %                  sigma_tilde_new(:,5) = sigma_tilde_old(:,5)+delta_sigma(:,5);
  %                  hvar_new(1:ntens,5)        =   sigma_tilde_new(:,5) ;   %aeh_agosto_2013
  if (gamma == 1) % SD Range
      kuu = kuu  + (B'* c_tilde5 *B)*m_pesoPG(5);
  end
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % SINGULAR POINT
  % STRAIN & STRESS TENSORS ON THE 6th GAUSS POINT (LOAD POINT)
  % ***********************************************************
  B               = m_Be(:,:,6);
%alf14  m_DefMACRO      = eps_old(:,6) + B*delta_u;
  m_DefMACRO      =   B*delta_u;

  % ELEMENT SELECTION
    if ind_ActState_new == 1 || (ind_ActState_new == 2 && ~C_inyect_SD)       % MSL ELEMENT
%alf      m_DefMACRO = eps_old(:,6) + B*delta_u + gamma*(1/ksd)*Dbeta_e_MSL;
      m_DefMACRO = B*delta_u + gamma*(1/ksd)*Dbeta_e_MSL;
  elseif ind_ActState_new == 2 && C_inyect_SD
      % SDA ELEMENT
%alf      m_DefMACRO = eps_old(:,6) + B*delta_u + gamma*((1/ksd)*n_tensor-Phi_Grad)*Dbeta_e_SDA;
%%      m_DefMACRO = B*delta_u + gamma*((1/ksd)*n_tensor-Phi_Grad)*Dbeta_e_SDA;
      m_DefMACRO = B*delta_u + gamma*(-Phi_Grad)*Dbeta_e_SDA;
  end
  
  if (BI_n == 1); data_micro.e_VG_MICRO.bandwidth  = ksd ; end
  switch conshyp % Modelo constitutivo
      case 11
          m_DefMACRO = m_DefMACRO + gamma/ksd*n_tensor*Dbeta_e_SDA;
          d_eps_new(:,6) = m_DefMACRO;
          eps_new(:,6)   = eps_old(:,6) + d_eps_new(:,6) ;
          % ONE SCALE MODEL
          [ctR,sigma_new(:,6),sigma_new_impl(:,6),hvar_new(:,6)] = ...
          rmap_damage_IMPLEX (...
              eps_new(:,6),hvar_old(:,6),e_DatMatSet,Eprop_qSD);
          %VARIABLE DE INYECCION - DISIPACION DE ENERGIA
          vectVHElem_new(7) = 0.5*((1-hvar_new(1,6))*eps_new(:,6)'*e_DatMatSet.ce*eps_new(:,6));  % energia interna
          val_1 = sigma_new(:,6)'*(eps_new(:,6) - eps_old(:,6));
          val_2 =(vectVHElem_new(7)-vectVHElem_old(7));
          if val_1 >= val_2
              vectVHElem_new(1) = vectVHElem_old(1) + ksd* (val_1 - val_2);
          else
              vectVHElem_new(1) = vectVHElem_old(1);
          end
          vectVHElem_new(4) = hvar_new(4,6);
          %vectVHElem_new(5) = hvar_new(4,6);   %hvar_new(5,6)-hvar_old(5,6);
          vectVHElem_new(5) = hvar_new(3,6);   %hvar_new(5,6)-hvar_old(5,6);
          % *******************************************
          delta_sigma(:,6) = sigma_new(:,6)-sigma_old(:,6);
          m_d_TraccNew     = n_tensor'*delta_sigma(:,6);
          
          if ind_ActState_new == 2  && C_inyect_SD
              m_CTHomog_ub= n_tensor'*ctR.Implex;
              m_CTHomog_bb= n_tensor'*ctR.Implex*n_tensor/ksd;
              ctR.Implex(1:2,1:2)   = m_CTHomog_bb;
              ctR.Implex(3:4,:) = m_CTHomog_ub;
          end
      case 54   % TWO SCALE MODEL  SANTA FE PG 6
          e_VG.iPG   = 6;
          e_VG.elast = 0;
%alf14          m_IDefR = eps_new(:,6);
          d_eps_new(:,6) = m_DefMACRO;
  
          m_IDefR = m_DefMACRO;;
          m_IDefS = Dbeta_e_SDA;
          %
          [ctR,sigma_new(:,6),hvar_new(:,6),vectVHElem_new,m_d_TraccNew] =...
              f_RMap_MEStafe(m_IDefR,m_IDefS,hvar_old(:,6),...
                  e_DatMatSet,BI_n,ind_ActState_new,C_inyect_SD,...
                  e_VG,vectVHElem_old,kinf);
          % *******************************************
%          eps_new(:,6)   = hvar_new(:,6).c_DefMacro ;

          delta_sigma(:,6) = sigma_new(:,6)-sigma_old(:,6);
      %    m_d_TraccNew   = delta_sigma([1:2],6);
          if ind_ActState_new == 2  && C_inyect_SD
%               m_CTHomog_ub= n_tensor'*ctR.Implex;
%               m_CTHomog_bb= n_tensor'*ctR.Implex*n_tensor/ksd;
%               ctR.Implex(1:2,1:2)   = m_CTHomog_bb;
%               ctR.Implex(3:4,:)     = m_CTHomog_ub;
          end
  end
  
  c_tilde6              = ctR.Implex;
  m_TensorTang(:,:,6)   = c_tilde6 ;
  m_TensorTang(:,:,12)  = ctR.Impli;
  
  if gamma ==1
      kuu = kuu  + (B'*c_tilde6*B)*m_pesoPG(6);
  end

  
  % REVISION DE LOS PARAMETROS DE INYECCION DE LA SD
  % ************************************************
  if (BI_n == 0) % SD Range
      if TWO_SCALE ==1;
          iSet= 1;  % modificar esto !!!!!
          FRAC_ENER = e_DatMatSet.e_DatSet(iSet).e_DatMat.gfv;
      else
          %% FRAC_ENER = Eprop(1,14);
          FRAC_ENER       = e_DatMatSet.gfv;
      end
      %%%          hvar_new(8,6)= 0.22*FRAC_ENER;  % FRACTURE ENERGY CRITERIA
      vectVHElem_new(2) = fact_inyect*FRAC_ENER;  % FRACTURE ENERGY CRITERIA
    %  vectVHElem_new(2) = 1e9*FRAC_ENER;  % FRACTURE ENERGY CRITERIA
      vectVHElem_new(3) =  100000.0*fact_inyect*vectVHElem_new(2); % RANGO NO USADO (POR SI SE QUIERE USAR 3 TIPOS DE ELEMENTOS FINITOS)
      % hvar_new(8,6)= 0.02*FRAC_ENER;  % FRACTURE ENERGY CRITERIA
      % hvar_new(9,6)= 100000.0*hvar_new(8,6); % RANGO NO USADO (POR SI SE QUIERE USAR 3 TIPOS DE ELEMENTOS FINITOS)
  else
      
      vectVHElem_new(2) = vectVHElem_old(2)   ; %injection  % en el PG=6 (singular) almaceno las variables elementales
      vectVHElem_new(3) = vectVHElem_old(3)   ; %Posibilidad de un tercer tipo de EF
      %  hvar_new(8,6)=hvar_old(8,6);   %injection  % en el PG=6 (singular) almaceno las variables elementales
      %  hvar_new(9,6)=hvar_old(9,6);   %Posibilidad de un tercer tipo de EF
      
  end
  
  %------------------
  Eprop_qSD.ksd=h_eq;
  %------------------
  
  % STRAIN & TENSION TENSORS ON THE STANDARD (REGULAR) GAUSS POINTS (1-4)
  % *********************************************************************
  for iPG = 1:(npg-2) % swap the initial four gauss points - regular points
      B = m_Be(:,:,iPG);
%alf14      m_DefMACRO  = eps_old(:,iPG) + B*delta_u;
      m_DefMACRO  = B*delta_u;

      % ELEMENT SELECTION
    if ind_ActState_new == 1 || (ind_ActState_new == 2 && ~C_inyect_SD)  
          % MSL ELEMENT
          m_DefMACRO = B*delta_u - gamma*(1/(h_eq-ksd))*Dbeta_e_MSL;
  elseif ind_ActState_new == 2 && C_inyect_SD
          % SDA_ELEMENT
          m_DefMACRO = B*delta_u - gamma*Phi_Grad*Dbeta_e_SDA;
      end
      
%alf14      eps_new(:,iPG)  = m_DefMACRO;
      d_eps_new(:,iPG)  = m_DefMACRO;
      eps_new(:,iPG)    = eps_old(:,iPG) + d_eps_new(:,iPG);
      switch conshyp % Modelo constitutivo
          case 11
              [ctR,sigma_new(:,iPG),sigma_new_impl(:,iPG),...
                  hvar_new(:,iPG)] = rmap_damage_IMPLEX (...
                  eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,Eprop_qSD);
          case 54   % TWO SCALE MODEL  PG : iPG
              e_VG.iPG   = iPG;
          %    if gamma   == 0
                  e_VG.elast = 0;
          %    else
          %        e_VG.elast = 1;
          %    end
              
              m_IDefR = d_eps_new(:,iPG);
              m_IDefS = zeros(ntens,1);
              %
              [ctR,sigma_new(:,iPG),hvar_new(:,iPG),~, ~] =...
                  f_RMap_MEStafe(m_IDefR,m_IDefS,hvar_old(:,iPG),...
                  e_DatMatSet,condBifR,ind_ActState_new,C_inyect_SD,...
                  e_VG,vectVHElem_old,kinf);
      end
      c_tilde                   = ctR.Implex;
      m_TensorTang(:,:,iPG)     = c_tilde ;
      m_TensorTang(:,:,iPG+npg) = ctR.Impli;

      if (gamma == 0) % SD Range
          kuu = kuu  + tau*(B'* c_tilde *B)*m_pesoPG(iPG);
      end
      
      
      
      delta_sigma(:,iPG) = (1-gamma)*(tau*(sigma_new(:,iPG)-sigma_old(:,iPG))+...
          (1-tau)*(sigma_new(:,6)-sigma_old(:,6)))...
          + ...
          gamma*(psi*(sigma_new(:,6)-sigma_old(:,6)) + (1-psi)*(sigma_new(:,5)-sigma_old(:,5)));
      
      sigmaTilde_new(:,iPG) = sigmaTilde_old(:,iPG) + delta_sigma(:,iPG);

      % Calculo de fint
      % ***************
      % fint = fint+B'*delta_sigma(:,iPG)*m_pesoPG(iPG);
      fint = fint+B'*sigmaTilde_new(:,iPG)*m_pesoPG(iPG);
  end
 
  
 % VARIABLE AUXIILAR PARA RECUPERACION DE EF INYECTADO 
% vectVHElem_new(8) = hvar_new(3,1:4)*hvar_new(3,1:4)';

%  if hvar_new(3,1)>0 && hvar_new(3,2)>0 && hvar_new(3,3)>0 && hvar_new(3,4)>0
%      vectVHElem_new(8) = hvar_new(3,1:4)*hvar_new(3,1:4)';
%  else
%      vectVHElem_new(8) = 0;
%  end
        
 if (gamma == 1)
      s=area/h_eq;
      % ELEMENT SELECTION
      if ind_ActState_new == 1  || (ind_ActState_new == 2 && ~C_inyect_SD) 
          % MSL ELEMENT
          % ***********
          c6_menos_c5 = (c_tilde6-c_tilde5);
          %%           kub_MSL = m_Be(:,:,5)'*(c_tilde5*(1/h_eq)*m_pesoPG(5) + c_tilde6*(-1/ksd)*m_pesoPG(6));
          kub_MSL = m_Be(:,:,5)'*c6_menos_c5 *s ;
          kbu_MSL = s*c6_menos_c5 *m_Be(:,:,5);
          kbbIN_MSL = inv(s*( 1/(h_eq-ksd) *c_tilde5 + 1/ksd*c_tilde6 ));
          
          % TRACTION CONTINUITY OVER THE LOCALIZATION BAND
          % **********************************************
          res_beta_MSL = s*(delta_sigma(:,6)-delta_sigma(:,5));
          
          % K MATRIX AND RESIDUAL CONDENSED FORCE
          % *************************************
          kuu = kuu-kub_MSL*kbbIN_MSL*kbu_MSL;
          fint = fint-kub_MSL*kbbIN_MSL*res_beta_MSL;
      elseif ind_ActState_new == 2  && C_inyect_SD
           % SDA ELEMENT
           % ***********
%alf14           kub_SDA =m_Be(:,:,5)'*(-c_tilde5*Phi_Grad*m_pesoPG(5) + c_tilde6*((1/ksd)*n_tensor-Phi_Grad)*m_pesoPG(6));
           kub_SDA = -m_Be(:,:,5)'*c_tilde5*Phi_Grad*m_pesoPG(5) ;
%alf14           kbu_SDA = s*n_tensor'*(c_tilde6*m_Be(:,:,6) - c_tilde5 *m_Be(:,:,5));
           m_Dctn  = c_tilde6(ndime+1:2*ndime,:) +   n_tensor'*( -c_tilde5);
           kbu_SDA = s*m_Dctn *m_Be(:,:,5) ;
%alf14           kbbIN_SDA = inv(s*n_tensor'*(c_tilde6*(-Phi_Grad + (1/ksd)*n_tensor) + c_tilde5 *Phi_Grad));
           kbbIN_SDA = (s* (c_tilde6(1:ndime, 1:ndime)- m_Dctn*Phi_Grad ) )\eye(ndime);
 
           % TRACTION CONTINUITY OVER THE LOCALIZATION BAND
           % **********************************************
%alf14           res_beta_SDA = s*n_tensor'*(delta_sigma(:,6) - delta_sigma(:,5));
           
%%          res_beta_SDA = s*(m_d_TraccNew - n_tensor'*delta_sigma(:,5)) ;
          res_beta_SDA = s*(m_d_TraccNew - n_tensor'*sigma_new(:,5) ) ;

           
           % K MATRIX AND RESIDUAL CONDENSED FORCE
           % *************************************
           m_k1 = kub_SDA*kbbIN_SDA;
           kuu  = kuu -  m_k1*kbu_SDA;
           
           fint = fint - m_k1*res_beta_SDA;
       
       end
     
   end
   
   %Se ordena las matrices como vectores columnas (es mas rapido usar Mat(:), que reshape(Mat,[],1);).
   sigma_new = sigma_new(:);
   sigmaTilde_new = sigmaTilde_new(:);
   eps_new = eps_new(:);
   hvar_new = hvar_new(:);
%   sigma_tilde_new = sigma_tilde_new(:);
   
end