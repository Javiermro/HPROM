function [kt,fint,sigma_new,eps_new,eps_fluct,hvar_new,aux_var,m_TensorTang] = ...
   f_MatElem_bbar_q1(u,hvar_old,aux_var,e_DatElemSet,e_DatMatSet,...
   m_Be,m_DetJe,DefMacro,sigma_old,e_VG)

% Variable globales
ntens = e_VG.ntens;
dofpe = e_DatElemSet.dofpe;
nPG = e_DatElemSet.npg; 
wg = e_DatElemSet.wg;

sihvarpg = e_DatMatSet.sihvarpg;
siavarpg = e_DatMatSet.siavarpg;
conshyp  = e_DatMatSet.conshyp;
esImplex = e_DatMatSet.esImplex;
      
% Inicializaciones
kt           = zeros(dofpe,dofpe);
fint         = zeros(dofpe,1);
sigma_new    = zeros(ntens,nPG);
eps_new      = zeros(ntens,nPG);
eps_fluct    = zeros(ntens,nPG);
%hvar_new     = zeros(sihvarpg,nPG);
hvar_new = f_DefinicionhVar(conshyp,sihvarpg,nPG);
%aux_var      = zeros(siavarpg,nPG);
if esImplex
   m_TensorTang = zeros(ntens,ntens,2*nPG);
else
   m_TensorTang = zeros(ntens,ntens,nPG);
end
hvar_old     = reshape(hvar_old,sihvarpg,[]);
%hvar_new = reshape(hvar_new,sihvarpg,[]);
aux_var = reshape(aux_var,siavarpg,nPG);
%sigma_old    = reshape(sigma_old,ntens,nPG);
%DefMacro     = reshape(DefMacro,ntens,[]);
%m_pesoPG     = m_DetJe.*wg*thickness;
%En el peso de gauss ya viene multiplicado el espesor.
m_pesoPG = m_DetJe.*wg;
%volElem = sum(m_pesoPG);
 
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
      case 2   %ELASTO - PLASTICIDAD J2: HARDENING-SOFTENING ISOTROPICO
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_plasJ2(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
      case 4   %Daño isótropico
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG); 
      case 9
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_plasJ2bl(...
            ntens,eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
      case 10  %Daño isotrópico regularizado
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage_reg(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
      case 11  %Daño isotrópico solo tracción regularizado
          if e_VG.elast                           
              [ct,sigma_new(:,iPG),~,hvar_new(:,iPG)] ...
                          = rmap_damage_elastic(eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet);
          else
              Eprop_qSD.ksd = e_DatMatSet.m_hReg(e_VG.iElemSet);
              [ct,sigma_new(:,iPG),~,hvar_new(:,iPG)] = rmap_damage_IMPLEX (eps_new(:,iPG),...
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
      case 50  %Modelo multiescala clásico          
%          %fprintf('**-- Inicio del return mapping del modelo multiescala\n')
%          cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
%          listeners = cmdWinDoc.getDocumentListeners;
%          jFxCommandArea = listeners(3);
%          %colorTextArea = get(jTextArea,'Background');
%          set(jFxCommandArea,'Background','red');
         [ct,sigma_new(:,iPG),hvar_new(:,iPG)] = f_RMap_ME(eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
%          set(jFxCommandArea,'Background','yellow');
%          %fprintf('**-- Fin del return mapping del modelo multiescala\n')
      case 52  %Elasticidad usando el tensor elástico obtenido de una homogenización de una microcelda.
         [ct,sigma_new(:,iPG)] = rmap_elast(eps_new(:,iPG),e_DatMatSet.ce);
      otherwise
         error('Matrices Elementales BBar_q1: Modelo constitutivo no definido.')
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
   %fint = fint+B'*(sigma_new(:,iPG)-sigma_old(:,iPG))*m_pesoPG(iPG);
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