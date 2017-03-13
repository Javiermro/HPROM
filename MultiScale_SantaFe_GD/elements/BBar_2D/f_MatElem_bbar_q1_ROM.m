function [kt,fint,ROMI,sigma_new,eps_new,eps_fluct,hvar_new,aux_var,m_TensorTang,IntEnergyNew] = ...
 f_MatElem_bbar_q1_ROM(kin_Var,eta,ModoPHI_EPS,ROMI,eps_old,hvar_old,aux_var,e_DatElemSet,e_DatMatSet,...
      m_Be,m_DetJe,DefMacro,sigma_old,leq,IntEnergyOld,e_VG)

% Variable globales
%conshyp  = e_VG.conshyp;
%npe      = e_VG.npe;
%dofpe    = e_VG.dofpe;
%npg      = e_VG.npg;
%wg       = e_VG.wg;
ntens = e_VG.ntens;
%sihvarpg = e_VG.sihvarpg;
%siavarpg = e_VG.siavarpg;
dofpe = e_DatElemSet.dofpe;
nPG = e_DatElemSet.npg; 
wg = e_DatElemSet.wg;

sihvarpg = e_DatMatSet.sihvarpg;
siavarpg = e_DatMatSet.siavarpg;
conshyp  = e_DatMatSet.conshyp;
esImplex = e_DatMatSet.esImplex;

if e_VG.MOD_TYPE == 1
    kin_Var =  reshape(kin_Var,ntens,[]);   
end

% Inicializaciones
kt           = zeros(dofpe,dofpe);
fint         = zeros(dofpe,1);
sigma_new    = zeros(ntens,nPG);
eps_new      = zeros(ntens,nPG);
eps_fluct    = zeros(ntens,nPG);
IntEnergyNew  = zeros(nPG,1) ;
%hvar_new     = zeros(sihvarpg,npg);
hvar_new = f_DefinicionhVar(conshyp,sihvarpg,nPG);
%aux_var      = zeros(siavarpg,npg);
if esImplex
   m_TensorTang = zeros(ntens,ntens,2*nPG);
else
   m_TensorTang = zeros(ntens,ntens,nPG);
end
hvar_old     = reshape(hvar_old,sihvarpg,[]);
%hvar_new = reshape(hvar_new,sihvarpg,[]);
aux_var = reshape(aux_var,siavarpg,nPG);

sigma_old    = reshape(sigma_old,ntens,nPG);
eps_old    = reshape(eps_old,ntens,nPG);

%DefMacro     = reshape(DefMacro,ntens,[]);
%m_pesoPG     = m_DetJe.*wg*thickness;
Reg_PAR.leq = leq;

%En el peso de gauss ya viene multiplicado el espesor.
m_pesoPG = m_DetJe.*wg;
 
for iPG = 1:nPG
   
   e_VG.iPG = iPG;

   if e_VG.MOD_TYPE == 1
       eps_fluct(:,iPG) = kin_Var(:,iPG) ;
   else
       B = m_Be(:,:,iPG);
       eps_fluct(:,iPG) = B*kin_Var;
   end

   % Deformacion aplicada a cada punto de Gauss
   eps_new(:,iPG) = DefMacro(:,iPG) + eps_fluct(:,iPG);
   %eps_new(:,iPG) = DefMacro+eps_fluct(:,iPG);
            
   % Modelo constitutivo
   switch conshyp
      case 1   %Elasticidad
         [ctR,sigma_new(:,iPG)] = rmap_elast(eps_new(:,iPG),e_DatMatSet.ce);
      case 2   %ELASTO - PLASTICIDAD J2: HARDENING-SOFTENING ISOTROPICO
         [ctR,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_plasJ2(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
      case 4   %Danio isotropico
         [ctR,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG); 
      case 9
         [ctR,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_plasJ2bl(...
            ntens,eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
      case 10  %Danio isotropico regularizado
         %Para obtener una longitud del elemento se considera como un prisma recto de altura 1 con
         %bases de cuadrados.
         %hElemReg = sqrt(volElem/1);
         [ctR,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage_reg(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
       case 11  %Danio isotropico solo tracci�n regularizado
           [ctR,sigma_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage_BCN(...
               eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_DatElemSet,e_VG,Reg_PAR);
       case 12  %Danio isotropico solo tracci�n regularizado
           [ctR,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] =RMapDanoSoloTraccionReg(...
               eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_DatElemSet,e_VG,Reg_PAR);
       case 13 %Danio isotropico solo tracci�n regularizado
           [ctR,sigma_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage_RANKINE_BCN(...
               eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_DatElemSet,e_VG,Reg_PAR);
      case 50  %Modelo multiescala clasico          
%          %fprintf('**-- Inicio del return mapping del modelo multiescala\n')
%          cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
%          listeners = cmdWinDoc.getDocumentListeners;
%          jFxCommandArea = listeners(3);
%          %colorTextArea = get(jTextArea,'Background');
%          set(jFxCommandArea,'Background','red');
         [ctR,sigma_new(:,iPG),hvar_new(:,iPG)] = f_RMap_ME(eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
%          set(jFxCommandArea,'Background','yellow');
%          %fprintf('**-- Fin del return mapping del modelo multiescala\n')
      case 52  %Elasticidad usando el tensor elastico obtenido de una homogenizacion de una microcelda.
         [ctR,sigma_new(:,iPG)] = rmap_elast(eps_new(:,iPG),e_DatMatSet.ce);
      otherwise
         error('Matrices Elementales BBar_q1: Modelo constitutivo no definido.')
   end
      
   %Se almacena para los tensor tangente constitutivos para realizar homogeneizacion y analisis de
   %bifurcacion
   if esImplex
      % IMPLEX: En el caso de tener activado el implex se almacena las matrices constitutivas
      %tangentes implicitas para el analisis de bifurcacion.
      %Por ahora se asume que el modelo constitutivo siempre devuelve una estructura cuando est� seleccionado
      %el implex.
      m_TensorTang(:,:,iPG) = ctR.Implex;
      %Se almacena los tensores tangentes constitutivo implicitos para analisis de bifurcacion como si fuera
      %PG adicionales, tantos como nPG. Se almacena en los indices (:,:,nPG+1:2*nPG).
      %Hay que tener cuidado cuando se hace un loop por la componente de PG sea hasta nPG, y no usando la
      %tercera dimension de m_TensorTang (size(m_TensorTang,3)).
      m_TensorTang(:,:,iPG+nPG) = ctR.Impli;
      %En los calculos para el ensamblaje se utiliza el implex.
      ct = ctR.Implex;
   else
      ct = ctR;
      m_TensorTang(:,:,iPG) = ct;
   end
   
   % Calculo de la energia interna
   switch conshyp
       case 1 % Elastic material
           IntEnergyNew(iPG) = 0;
       case 2 % Plastic material W = 0.5*(e_e*C*e_e)
           %IntEnergyNew(iPG) = 0.5*((eps_new(:,iPG)-hvar_new(1:4,iPG))'*e_DatMatSet.ce*(eps_new(:,iPG)-hvar_new(1:4,iPG)));
           IntEnergyNew(iPG) = IntEnergyOld(iPG) + 0.5*(sigma_new(:,iPG)+sigma_old(:,iPG))'*((eps_new(:,iPG)-hvar_new(1:4,iPG))-(eps_old(:,iPG)-hvar_old(1:4,iPG))) ...
               - 0.5*(hvar_old(12,iPG)+hvar_new(12,iPG))*(hvar_new(5,iPG)-hvar_old(5,iPG));
       case 11 % Continuum damage model
           IntEnergyNew(iPG) = 0.5*((1-aux_var(1,iPG))*eps_new(:,iPG)'*...
               e_DatMatSet.ce*eps_new(:,iPG));  % energia interna
       case 50 % Multiscale classic model
           IntEnergyNew(iPG) = 0.5*(sigma_new(:,iPG)'*eps_new(:,iPG));  % energia interna
       otherwise
           IntEnergyNew(iPG) = 0.5*((1-aux_var(1,iPG))*eps_new(:,iPG)'*...
               e_DatMatSet.ce*eps_new(:,iPG));  % energia interna
   end
   
   iPG_Val = (iPG-1)*ntens+1:iPG*ntens ;
   if e_VG.MOD_TYPE ==1
       
       eta_voigt = e_VG.E_MATRIX*[eta(1) eta(2) 0 eta(3)]';
       
       % Calculo de fint
       ROMI.Fint_ROM = ROMI.Fint_ROM + ModoPHI_EPS(iPG_Val,:)'*(sigma_new(:,iPG)+eta_voigt)*m_pesoPG(iPG);
       
       % Calculo de Kep implex
       ROMI.KT_ROM = ROMI.KT_ROM + ModoPHI_EPS(iPG_Val,:)'*ct*ModoPHI_EPS(iPG_Val,:)*m_pesoPG(iPG);
       
       if esImplex % Implex
           
           % Calculo de Kep implex
           ROMI.KT_ROM_IMP  = ROMI.KT_ROM_IMP + ModoPHI_EPS(iPG_Val,:)'*ctR.Impli*ModoPHI_EPS(iPG_Val,:)*m_pesoPG(iPG);
           
           % right hand vector for Homogenized Tangent tensor (implex)
           ROMI.Q_ROM = ROMI.Q_ROM + ModoPHI_EPS(iPG_Val,:)'*ct*m_pesoPG(iPG);
           
           % right hand vector for Homogenized Tangent tensor (implicit)
           ROMI.Q_ROM_IMP = ROMI.Q_ROM_IMP + ModoPHI_EPS(iPG_Val,:)'*ctR.Impli*m_pesoPG(iPG);
           
       else % Caso implicito
           
           % Calculo de Kep implex
           ROMI.KT_ROM_IMP  = ROMI.KT_ROM_IMP + ModoPHI_EPS(iPG_Val,:)'*ct*ModoPHI_EPS(iPG_Val,:)*m_pesoPG(iPG);
           
           % right hand vector for Homogenized Tangent tensor (implex)
           ROMI.Q_ROM = ROMI.Q_ROM + ModoPHI_EPS(iPG_Val,:)'*ct*m_pesoPG(iPG);
           
           % right hand vector for Homogenized Tangent tensor (implicit)
           ROMI.Q_ROM_IMP = ROMI.Q_ROM_IMP + ModoPHI_EPS(iPG_Val,:)'*ct*m_pesoPG(iPG);
           
       end
       
       % % DISSIPATION ASSESSMENT
       % IntEnergyNew(iPG) = 0.5*((1-aux_var(1,iPG))*eps_new(:,iPG)'*...
       %      e_DatMatSet.ce*eps_new(:,iPG)) + eta_voigt'*kin_Var(:,iPG);  % energia interna expandida
       
   else
       
       % Calculo de fint
       fint = fint+B'*sigma_new(:,iPG)*m_pesoPG(iPG);
       
       % Calculo de Kep
       kt = kt+B'*ct*B*m_pesoPG(iPG);
   end
   
   
end

% Se ordena las matrices como vectores columnas
sigma_new  = sigma_new(:);
hvar_new   = hvar_new(:);
aux_var    = aux_var(:);
eps_new    = eps_new(:);
eps_fluct  = eps_fluct(:);
   
end
