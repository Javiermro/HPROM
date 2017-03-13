function [e_VarEst_new,e_VarAux] = f_OperPosConv_BNS(u,xx,e_VarEst_new,...
   e_VarEst_old,e_VarAux,e_DatSet,c_CT,e_VG)

for iSet = 1:e_VG.nSet
   e_DatElemSet = e_DatSet(iSet).e_DatElem;
   e_DatMatSet = e_DatSet(iSet).e_DatMat;
   eltype = e_DatElemSet.eltype;
   conshyp = e_DatMatSet.conshyp;
   conec = e_DatSet(iSet).conec;
   m_DofElem = e_DatSet(iSet).m_DofElem;
   dofElemSet = m_DofElem(:);
   nElem = e_DatSet(iSet).nElem;
   uElemSet = reshape(u(dofElemSet),[],nElem);
   %
   sigma_new = e_VarEst_new(iSet).sigma;
   eps_new = e_VarEst_new(iSet).eps;
   eps_fluct_new = e_VarEst_new(iSet).eps_fluct;
   hvar_new = e_VarEst_new(iSet).hvar;
   hvar_old = e_VarEst_old(iSet).hvar;
   m_VarAuxGP = e_VarAux(iSet).VarAuxGP;
   m_VarAuxElem = e_VarAux(iSet).VarAuxElem;
   m_VarHistElem =  e_VarEst_new(iSet).VarHistElem;
   
   switch eltype
      case 2
         switch conshyp
            case 8
               e_VarEst_new(iSet).hvar = f_ReCalcNormal(sigma_new,hvar_new,e_DatMatSet,e_VG);
            case 12
               %Determinación del ángulo de bifurcación (por ahora no se realiza ningún análisis, solo se
               %devuelve el ángulo asumiendo que está bifurcado).
               e_VarEst_new(iSet).hvar = f_AnalisisBif(eps_new,hvar_new,e_DatSet(iSet),e_VG);
            case 55  %MULTIESCALA MODELO CLÁSICO CON ANÁLISIS DE BIFURCACIÓN
               m_CT = c_CT{iSet};
               f_OperPosConv_MEClasBif(eps_new,m_CT,e_DatSet(iSet),e_VG);
         end
      case 4
         switch conshyp
            case 12
               e_VarEst_new(iSet).hvar = f_AnalisisBif(eps_new,hvar_new,e_DatSet(iSet),e_VG);
            case 55  %MULTIESCALA MODELO CLÁSICO CON ANÁLISIS DE BIFURCACIÓN
               m_CT = c_CT{iSet};
               f_OperPosConv_MEClasBif(m_CT,e_DatSet(iSet),e_VG);
            otherwise
               %%                  error('error f_OperPosConv_BNS: Bif.Modelo constitutivo no definido.')
         end
      case 10
         switch conshyp
            case {10,11,12}
               m_BT = e_DatSet(iSet).m_BT;
               m_CT = c_CT{iSet};
               [sigma_new,eps_new,eps_fluct_new,hvar_new,m_VarAuxElem] = ...
                  f_OperPosConvSDA_tria_t1_dano(...
                  sigma_new,eps_new,eps_fluct_new,hvar_new,m_VarAuxElem,m_CT,nElem,m_BT,...
                  e_DatMatSet,e_DatElemSet,e_DatSet(iSet).m_NumElem,e_VG);
               e_VarEst_new(iSet).sigma = sigma_new;
               e_VarEst_new(iSet).eps = eps_new;
               e_VarEst_new(iSet).eps_fluct = eps_fluct_new;
               e_VarEst_new(iSet).hvar = hvar_new;
               %e_VarAux(iSet).VarAuxGP = m_VarAuxGP;
               e_VarAux(iSet).VarAuxElem = m_VarAuxElem;
            case 51
               m_BT = e_DatSet(iSet).m_BT;
               m_CT = c_CT{iSet};
               m_ElemPGImpr = e_DatMatSet.m_ElemPGImpr;
               [sigma_new,eps_new,eps_fluct_new,hvar_new,m_VarHistElem,m_VarAuxGP,m_VarAuxElem] = ...
                  f_OperPosConvSDA_tria_t1_MECohesivo(...
                  sigma_new,eps_new,eps_fluct_new,hvar_new,hvar_old,m_VarHistElem,m_VarAuxGP,...
                  m_VarAuxElem,m_CT,nElem,m_BT,e_DatMatSet,e_DatElemSet,...
                  e_DatSet(iSet).m_NumElem,m_ElemPGImpr,e_VG);
               %
               e_VarEst_new(iSet).sigma = sigma_new;
               e_VarEst_new(iSet).eps = eps_new;
               e_VarEst_new(iSet).eps_fluct = eps_fluct_new;
               e_VarEst_new(iSet).hvar = hvar_new;
               e_VarAux(iSet).VarAuxGP = m_VarAuxGP;
               e_VarAux(iSet).VarAuxElem = m_VarAuxElem;
         end
      case 20
         m_CT = c_CT{iSet};
         m_VarAuxElem = f_OperPosConv_MixStrInj_quad_q1(m_VarAuxElem,m_VarAuxGP,e_DatMatSet,m_CT,e_VG);
         e_VarAux(iSet).VarAuxElem = m_VarAuxElem;
         %%%%%%%%%%%%%%%%%%%%%
      case {21,22,23}
         
         m_VarHistElem_old =  e_VarEst_old(iSet).VarHistElem;
         EF_sidesCutCPF = e_DatSet(iSet).e_DatElem.EF_sidesCutCPF ;
         
         switch conshyp
            case {11,53}
               m_CT = c_CT{iSet};
               [m_VarAuxElem,m_VarHistElem] = ...
                  f_OperPosConvSDA_QuadQ1...
                  (uElemSet,xx,conec,m_VarHistElem,m_VarHistElem_old,EF_sidesCutCPF,...
                  m_VarAuxElem,m_CT,e_DatSet(iSet),e_VG,hvar_new,hvar_old);
               
               e_VarAux(iSet).VarAuxElem = m_VarAuxElem;
               e_VarEst_new(iSet).VarHistElem = m_VarHistElem;
            case 54
               m_CT     = c_CT{iSet};
               hvar_old          =  e_VarEst_old(iSet).hvar;
               [hvar_new,m_VarHistElem,m_VarAuxGP,m_VarAuxElem] = ...
                  f_OperPosConvSDA_QSFe(uElemSet,xx,conec,EF_sidesCutCPF,...
                  hvar_new,hvar_old,...
                  m_VarHistElem,m_VarHistElem_old,m_VarAuxGP,...
                  m_VarAuxElem,m_CT,nElem,e_DatMatSet,e_DatElemSet,...
                  e_DatSet(iSet), e_VG);
               %
               e_VarEst_new(iSet).hvar   = hvar_new;
               e_VarAux(iSet).VarAuxGP   = m_VarAuxGP;
               e_VarAux(iSet).VarAuxElem = m_VarAuxElem;
               e_VarEst_new(iSet).VarHistElem = m_VarHistElem;
            case 51
               nElem = e_DatSet(iSet).nElem;
               m_BT = e_DatSet(iSet).m_BT;
               m_CT = c_CT{iSet};
               hvar_old = e_VarEst_old(iSet).hvar;
               m_ElemPGImpr = e_DatMatSet.m_ElemPGImpr;
               [sigma_new,eps_new,eps_fluct_new,hvar_new,m_VarHistElem,m_VarAuxGP,m_VarAuxElem] = ...
                  f_OperPosConvSDA_tria_t1_MECohesivo(...
                  sigma_new,eps_new,eps_fluct_new,hvar_new,hvar_old,m_VarHistElem,m_VarAuxGP,...
                  m_VarAuxElem,m_CT,nElem,m_BT,e_DatMatSet,e_DatElemSet,...
                  e_DatSet(iSet).m_NumElem,m_ElemPGImpr,e_VG);
               %
               e_VarEst_new(iSet).sigma = sigma_new;
               e_VarEst_new(iSet).eps   = eps_new;
               e_VarEst_new(iSet).eps_fluct = eps_fluct_new;
               e_VarEst_new(iSet).hvar = hvar_new;
               e_VarEst_new(iSet).VarHistElem = m_VarHistElem;
               e_VarAux(iSet).VarAuxGP = m_VarAuxGP;
               e_VarAux(iSet).VarAuxElem = m_VarAuxElem;
         end
      %%%%%%%%%%%%%%%%%%%%%
      %LARGE DEFORMATIONS
      case 108  %Cuadrángulo bilineal FBar
         switch conshyp
            case 55  %MULTIESCALA MODELO CLÁSICO CON ANÁLISIS DE BIFURCACIÓN
               m_CT = c_CT{iSet};
               f_OperPosConv_MEClasBif(eps_new,m_CT,e_DatSet(iSet),e_VG);
         end
         
   end
   
end