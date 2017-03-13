function [kt,fint,sigma_new,eps_new,eps_fluct,hvar_new,aux_var,m_TensorTang] = f_MatElem_MixStrInj_quad_q1(...
   u,hvar_old,m_VarAuxPG,m_VarAuxElem,e_DatElemSet,e_DatMatSet,m_Be,m_DetJe,DefMacro,sigma_old,e_VG)

   %% ESTE ELEMENTO FINITO NO ESTÁ FUNCIONANDO CORRECTAMENTE

   % Variable globales
   ntens = e_VG.ntens;
   %
   dofpe = e_DatElemSet.dofpe;
   nPG = e_DatElemSet.npg; 
   wg = e_DatElemSet.wg;

   % Propiedades materiales
   sihvarpg = e_DatMatSet.sihvarpg;
   siavarpg = e_DatMatSet.siavarpg;
   conshyp  = e_DatMatSet.conshyp;
   esImplex = e_DatMatSet.esImplex;
   
   % Variable auxiliar del elemento
   %Condición de bifurcación
   %Se determina mediante el análisis de bifurcación del PG central (5), una vez que este bifurcó todo el
   %elemento se considera bifurcado (y en consecuencia los restantes PGs).   
   condBif = m_VarAuxElem(1);

   % Inicializaciones
   kt = zeros(dofpe,dofpe);
   fint = zeros(dofpe,1);
   sigma_new = zeros(ntens,nPG);
   eps_new = zeros(ntens,nPG);
   eps_fluct = zeros(ntens,nPG);
   %hvar_new = zeros(sihvarpg,nPG);
   hvar_new = f_DefinicionhVar(conshyp,sihvarpg,nPG);
   if esImplex
      m_TensorTang = zeros(ntens,ntens,2*nPG);
   else
      m_TensorTang = zeros(ntens,ntens,nPG);
   end
   
   % Redimensionado de matrices
   hvar_old = reshape(hvar_old,sihvarpg,[]);
   aux_var = reshape(m_VarAuxPG,siavarpg,nPG);
   sigma_old = reshape(sigma_old,ntens,[]);
   
   %En el peso de gauss ya viene multiplicado el espesor.
   m_pesoPG = m_DetJe.*wg;
   
   %Variable condBif, por ahora se sigue la siguiente convención
   %condBif=1: Pertenece al dominio de elementos mixtos.
   %condBif=0: No pertenece al dominio de elementos mixtos.
   
   % Parámetros de estabilización por PG
   %Se multiplica los parámetros de estabilizacón por los pesos de los puntos de gauss directamente, así tener
   %en cuenta este parámetro en la integración. 
   if condBif>0
      %Si condBif es igual a 1, significa que se detectó bifurcación.
      %Ver si no vale la pena activar la subintegración unos pasos después que bifurcó (condBif=2), o eso se
      %deja para la SD.
      paramEstab = e_DatElemSet.estabBif;
   elseif condBif==0
      %Si condBif es cero, significa que todavía no se detectó la bifurcación, o que descargó.
      paramEstab = e_DatElemSet.estabNoBif;
   end
   paramEstabRest = 1-paramEstab;
   
   % Punto de gauss 5
   %Se considera que el PG 5 es el central.
   iPG = 5;
   e_VG.iPG = iPG;
   B = m_Be(:,:,iPG);
   eps_fluct(:,iPG) = B*u;
   % Deformacion aplicada a cada punto de Gauss
   eps_new(:,iPG) = DefMacro+eps_fluct(:,iPG);
   % Modelo constitutivo
   switch conshyp
      case 2   %ELASTO - PLASTICIDAD J2: HARDENING-SOFTENING ISOTROPIC
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_plasJ2(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
      case 4
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG); 
      case 10  %Daño isotrópico regularizado
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage_reg(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
      case 12  %Daño isotrópico solo tracción regularizado
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = RMapDanoSTraccReg(...
            eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),e_DatMatSet,e_VG);
      case 50  %Modelo multiescala clásico
         %fprintf('**-- Inicio de return mapping del modelo multiescala\n')
         %cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
         %listeners = cmdWinDoc.getDocumentListeners;
         %jFxCommandArea = listeners(3);
         %set(jFxCommandArea,'Background','red');
         [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = f_RMap_ME(...
            eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
         %set(jFxCommandArea,'Background','yellow');
         %fprintf('**-- Fin de return mapping del modelo multiescala\n')
      otherwise
         error(['Matrices Elementales Mixed Strain Injection Quad_q1: Punto de gauss 5: Modelo ',...
            'constitutivo no definido.'])
   end
   %Se almacena para los tensor tangente constitutivos para realizar homogeneización y análisis de
   %bifurcación
   if esImplex
      % IMPLEX: En el caso de tener activado el implex se almacena las matrices constitutivas tangentes 
      %implícitas para el análisis de bifurcación.
      %Por ahora se asume que el modelo constitutivo siempre devuelve una estructura cuando está
      %seleccionado el implex.
      m_TensorTang(:,:,iPG) = ct.Implex;
      %Se almacena los tensores tangentes constitutivo implícitos para análisis de bifurcación como si
      %fuera PG adicionales, tantos como nPG. Se almacena en los índices (:,:,nPG+1:2*nPG).
      %Hay que tener cuidado cuando se hace un loop por la componente de PG sea hasta nPG, y no usando la
      %tercera dimensión de m_TensorTang (size(m_TensorTang,3)).
      m_TensorTang(:,:,iPG+nPG) = ct.Impli;
      % Aporte al tensor tangente constitutivo del PG 5 a los otros puntos de gauss
      ct5Estab = paramEstabRest*ct.Implex;
   else
      m_TensorTang(:,:,iPG) = ct;
      % Aporte al tensor tangente constitutivo del PG 5 a los otros puntos de gauss
      ct5Estab = paramEstabRest*ct;
   end
   % Aporte de las tensiones del PG 5 a los puntos de gauss
   dSigma5Estab = paramEstabRest*(sigma_new(:,iPG)-sigma_old(:,iPG));
   %fint = fint+B*dSigma5Estab;
   B5 = B;

   %Se considera que el puntos de gauss 5 es el central y que peso nulo, por lo que directamente se lo
   %descarta.
   %No se considera como optimización que cuando paramEstab es cero no interesa la determinación de las
   %tensiones para los PG 1 a 4, ya que no solo se debería utilizar las tensiones con el valor del PG 5, sino
   %que como también se debe actualizar las variables internas no se puede evitar los llamados a las funciones
   %constitutivas.
   for iPG = 1:nPG-1

      e_VG.iPG = iPG;

      B = m_Be(:,:,iPG);

      eps_fluct(:,iPG) = B*u;

      % Deformacion aplicada a cada punto de Gauss
      eps_new(:,iPG) = DefMacro+eps_fluct(:,iPG);

      % Modelo constitutivo
      switch conshyp
         case 2   %ELASTO - PLASTICIDAD J2: HARDENING-SOFTENING ISOTROPIC
            [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_plasJ2(...
               eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
         case 4
            [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage(...
               eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG); 
         case 10  %Daño isotrópico regularizado
            [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmap_damage_reg(...
               eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
         case 12  %Daño isotrópico solo tracción regularizado
            [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = RMapDanoSTraccReg(...
               eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),e_DatMatSet,e_VG);
         case 50  %Modelo multiescala clásico
            %fprintf('**-- Inicio de return mapping del modelo multiescala\n')
            %cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
            %listeners = cmdWinDoc.getDocumentListeners;
            %jFxCommandArea = listeners(3);
            %set(jFxCommandArea,'Background','red');
            [ct,sigma_new(:,iPG),eps_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = f_RMap_ME(...
               eps_new(:,iPG),hvar_old(:,iPG),e_DatMatSet,e_VG);
            %set(jFxCommandArea,'Background','yellow');
            %fprintf('**-- Fin de return mapping del modelo multiescala\n')
         otherwise
            error('Matrices Elementales Quad_q1: Modelo constitutivo no definido.')
      end

      %Se almacena para los tensor tangente constitutivos para realizar homogeneización y análisis de
      %bifurcación
      if esImplex
         % IMPLEX: En el caso de tener activado el implex se almacena las matrices constitutivas tangentes 
         %implícitas para el análisis de bifurcación.
         %Por ahora se asume que el modelo constitutivo siempre devuelve una estructura cuando está
         %seleccionado el implex.
         m_TensorTang(:,:,iPG) = ct.Implex;
         %Se almacena los tensores tangentes constitutivo implícitos para análisis de bifurcación como si
         %fuera PG adicionales, tantos como nPG. Se almacena en los índices (:,:,nPG+1:2*nPG).
         %Hay que tener cuidado cuando se hace un loop por la componente de PG sea hasta nPG, y no usando la
         %tercera dimensión de m_TensorTang (size(m_TensorTang,3)).
         m_TensorTang(:,:,iPG+nPG) = ct.Impli;
         %En los cálculos para el ensamblaje se utiliza el implex.
         ctEstab = paramEstab*ct.Implex;
      else
         m_TensorTang(:,:,iPG) = ct;
         ctEstab = paramEstab*ct;         
      end

      % Incremento de tensiones
      dSigmaiEstab = paramEstab*(sigma_new(:,iPG)-sigma_old(:,iPG));
      dSigma = dSigmaiEstab+dSigma5Estab;

      % Actualización de las tensiones estabilizadas
      sigma_new(:,iPG) = sigma_old(:,iPG)+dSigma;

      % Cálculo de fint
      %fint = fint+B'*dSigma*m_pesoPG(iPG);
      fint = fint+(B'*dSigmaiEstab+B5'*dSigma5Estab)*m_pesoPG(iPG);

      % Cálculo de matriz elemental
      %kt = kt+B'*ct*B*m_pesoPG(iPG);
      kt = kt+(B'*ctEstab*B+B5'*ct5Estab*B5)*m_pesoPG(iPG);

   end

   % Se ordena las matrices como vectores columnas
   sigma_new = sigma_new(:);
   hvar_new = hvar_new(:);
   aux_var = aux_var(:);
   eps_new = eps_new(:);
   eps_fluct = eps_fluct(:);

end