function [kt,fint,sigma_new,eps_new,eps_fluct,hvar_new,aux_var,m_TensorTang] = ...
   f_MatElem_quad_q1(u,eps_old,hvar_old,aux_var,e_DatElemSet,e_DatMatSet,...
   m_Be,m_DetJe,DefMacro,sigma_old,m_VolElem,e_VG)

% Variable globales
%conshyp  = e_VG.conshyp;
%npe      = e_VG.npe;
%dofpe    = e_VG.dofpe;
%nPG      = e_VG.nPG;
%wg       = e_VG.wg;
ntens    = e_VG.ntens;
dofpe = e_DatElemSet.dofpe;
nPG = e_DatElemSet.npg; 
wg = e_DatElemSet.wg;

sihvarpg = e_DatMatSet.sihvarpg;
siavarpg = e_DatMatSet.siavarpg;
conshyp  = e_DatMatSet.conshyp;
esImplex = e_DatMatSet.esImplex;

% Inicializaciones
kt            = zeros(dofpe,dofpe);
fint          = zeros(dofpe,1);
sigma_new     = zeros(ntens,nPG);
eps_new       = zeros(ntens,nPG);
eps_fluct     = zeros(ntens,nPG);
%hvar_new     = zeros(sihvarpg,nPG);
hvar_new      = f_DefinicionhVar(conshyp,sihvarpg,nPG);
%aux_var      = zeros(siavarpg,nPG);
if esImplex
   m_TensorTang = zeros(ntens,ntens,2*nPG);
else
   m_TensorTang = zeros(ntens,ntens,nPG);
end
hvar_old = reshape(hvar_old,sihvarpg,[]);
%hvar_new = reshape(hvar_new,sihvarpg,[]);
aux_var = reshape(aux_var,siavarpg,nPG);
%sigma_old = reshape(sigma_old,ntens,nPG);
%eps_old = reshape(eps_old,ntens,[]);
%DefMacro = reshape(DefMacro,ntens,[]);
%En el peso de gauss ya viene multiplicado el espesor.
m_pesoPG = m_DetJe.*wg;

for iPG = 1:nPG
   
   e_VG.iPG = iPG;

   B = m_Be(:,:,iPG);
   
   eps_fluct(:,iPG) = B*u;

   % Deformacion aplicada a cada punto de Gauss
   %eps_new(:,iPG) = DefMacro(:,iPG)+eps_fluct(:,iPG);
   eps_new(:,iPG) = DefMacro+eps_fluct(:,iPG);

   % Modelo constitutivo
   switch conshyp
      case 1   %Elasticidad
         [ct,sigma_new(:,iPG)] = rmap_elast(eps_new(:,iPG),e_DatMatSet.ce);
      case 2   %ELASTO - PLASTICIDAD J2: HARDENING-SOFTENING ISOTROPIC
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_plasJ2(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
      case 4
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG); 
      case 7
         [sigma_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmapfi_danio_plasdp_new(...
             eps_new(:,iPG),hvar_old(:,iPG),Eprop,ce,e_VG);
         ct = rmapct_danio_plasdp(eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),Eprop,ce,e_VG);
      case 10  %Daño isotrópico regularizado
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage_reg(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
      case 11  %Daño isotrópico solo tracción regularizado
          if e_VG.elast                           
              [ct,sigma_new(:,iPG),~,hvar_new(:,iPG)] ...
                          = rmap_damage_elastic(eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet);
          else
              Eprop_qSD.ksd = e_DatMatSet.m_hReg(e_VG.iElemSet);
              [ct,sigma_new(:,iPG),~,hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage_IMPLEX (eps_new(:,iPG),...
                  hvar_old(:,iPG),e_DatMatSet,Eprop_qSD);
          end
      case 12  %Daño isotrópico solo tracción regularizado
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] =...
             RMapDanoSTraccReg (...
               eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),e_DatMatSet,e_VG);
      case 13 %Daño isotrópico solo tracción regularizado
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = ...
            RMapDanoRankine(...
            eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),e_DatMatSet,e_VG);
      case {50,55}  %Modelo multiescala clásico
         %fprintf('**-- Inicio de return mapping del modelo multiescala\n')
         %cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
         %listeners = cmdWinDoc.getDocumentListeners;
         %jFxCommandArea = listeners(3);
         %set(jFxCommandArea,'Background','red');
         [ct,sigma_new(:,iPG),hvar_new(:,iPG)] = f_RMap_ME(eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
         %set(jFxCommandArea,'Background','yellow');
         %fprintf('**-- Fin de return mapping del modelo multiescala\n')
      case 52  %Elasticidad usando el tensor elástico obtenido de una homogenización de una microcelda.
         [ct,sigma_new(:,iPG)] = rmap_elast(eps_new(:,iPG),e_DatMatSet.ce);
      case 100 %Elastic Material neo-Hookean
         [ct,sigma_new(:,iPG)] = f_RMap_NeoHook_DefPl(eps_new(:,iPG),e_DatMatSet,e_VG);
         %Determinación numérica del tensor tangente constitutivo
         dF = 1e-8;
         m_DPTF = zeros(ntens,ntens);
         m_F = eps_new(:,iPG);
         m_Fdf = m_F;
         for i = 1:ntens
            m_Fdf(i) = m_F(i)+dF;
            [~,m_PT1] = f_RMap_NeoHook_DefPl(m_Fdf,e_DatMatSet,e_VG);
            m_Fdf(i) = m_F(i)-dF;
            [~,m_PT0] = f_RMap_NeoHook_DefPl(m_Fdf,e_DatMatSet,e_VG);
            %Diferencia finita centrada
            m_DPTF(:,i) = (m_PT1-m_PT0)/2/dF;
            %Se recupera el valor original para realizar los siguientes incrementos.
            m_Fdf(i) = m_F(i);
         end
         ct = m_DPTF;
      case 110 %Large deformations J2 Plasticity
         [ct,sigma_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG),strsg,dmatx] = ...
            f_RMapPlastJ2LD(eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),e_DatMatSet,e_VG);
      otherwise
         error('Matrices Elementales Quad_q1: Modelo constitutivo no definido.')
   end

   %Se almacena para los tensor tangente constitutivos para realizar homogeneización y análisis de
   %bifurcación
   if esImplex
      % IMPLEX: En el caso de tener activado el implex se almacena las matrices constitutivas
      %tangentes implícitas para el análisis de bifurcación.
      %Por ahora se asume que el modelo constitutivo siempre devuelve una estructura cuando está seleccionado
      %el implex.
      m_TensorTang(:,:,iPG) = ct.Implex;
      %Se almacena los tensores tangentes constitutivo implícitos para análisis de bifurcación como si fuera
      %PG adicionales, tantos como nPG. Se almacena en los índices (:,:,nPG+1:2*nPG).
      %Hay que tener cuidado cuando se hace un loop por la componente de PG sea hasta nPG, y no usando la
      %tercera dimensión de m_TensorTang (size(m_TensorTang,3)).
      m_TensorTang(:,:,iPG+nPG) = ct.Impli;
      %En los cálculos para el ensamblaje se utiliza el implex.
      ct = ct.Implex;
   else
      if isfield(ct,'Impli')  
         ct = ct.Impli;
      end
      m_TensorTang(:,:,iPG) = ct;
   end
   
   % Cálculo de fint
   %delta_sigma = sigma_new(:,iPG)-sigma_old(:,iPG);   
   %fint = fint+B'*delta_sigma*m_pesoPG(iPG);
   fint = fint+B'*sigma_new(:,iPG)*m_pesoPG(iPG);   
      
   % Cálculo de Kep
   kt = kt+B'*ct*B*m_pesoPG(iPG);
   
end
   
% Se ordena las matrices como vectores columnas
sigma_new  = sigma_new(:);
hvar_new   = hvar_new(:);
aux_var    = aux_var(:);
eps_new    = eps_new(:);
eps_fluct  = eps_fluct(:);
   
end