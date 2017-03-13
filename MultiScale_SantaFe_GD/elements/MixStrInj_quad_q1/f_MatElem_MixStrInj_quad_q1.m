function [kt,fint,sigma_new,eps_new,eps_fluct,hvar_new,m_VarHistElemNew,aux_var,m_TensorTang] = ...
   f_MatElem_MixStrInj_quad_q1(...
   u,hvar_old,m_VarHistElemOld,m_VarAuxPG,m_VarAuxElem,e_DatElemSet,e_DatMatSet,m_Be,m_DetJe,DefMacro,...
   sigma_old,e_VG)

   % Variable globales
   ntens = e_VG.ntens;
   
   dofpe = e_DatElemSet.dofpe;
   nPG = e_DatElemSet.npg; 
   wg = e_DatElemSet.wg;

   % Propiedades materiales
   sihvarpg = e_DatMatSet.sihvarpg;
   siavarpg = e_DatMatSet.siavarpg;
   conshyp  = e_DatMatSet.conshyp;
   esImplex = e_DatMatSet.esImplex;
   
   %Variable histórica del elemento
   m_SigmaEstabOld = reshape(m_VarHistElemOld,ntens,nPG);
   
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
   m_DSigmaEstab = zeros(ntens,nPG);
   
   % Redimensionado de matrices
   hvar_old = reshape(hvar_old,sihvarpg,[]);
   aux_var = reshape(m_VarAuxPG,siavarpg,nPG);
   sigma_old = reshape(sigma_old,ntens,[]);
   
   %En el peso de gauss ya viene multiplicado el espesor.
   %El término constante de este elemento, proveniente del problema mixto, se integra con 1 solo punto de
   %gauss, en el lugar de 4. Por ello al PG 5 se incorpora el peso 4, que corresponde al caso de un solo punto
   %de gauss (wg viene con cero para realizar correctamente otras integraciones dentro del código).
   wg(5) = 4;
   m_pesoPG = m_DetJe.*wg;
   
   %Variable condBif, por ahora se sigue la siguiente convención
   %condBif=1: Pertenece al dominio de elementos mixtos.
   %condBif=0: No pertenece al dominio de elementos mixtos.
   
   % Parámetros de estabilización por PG
   if condBif>0
      %Si condBif es igual a 1, significa que se detectó bifurcación.
      %Ver si no vale la pena activar la subintegración unos pasos después que bifurcó (condBif=2), o eso se
      %deja para la SD.
      factEstab = e_DatElemSet.estabBif;      
   elseif condBif==0
      %Si condBif es cero, significa que todavía no se detectó la bifurcación.
      factEstab = e_DatElemSet.estabNoBif;
   end 
   m_FactEstab = ones(nPG,1)*factEstab;
   %Se asume que el PG central está en la posición 5 de la lista de puntos de gauss.
   m_FactEstab(5) = 1-m_FactEstab(5);

   for iPG = 1:nPG

      e_VG.iPG = iPG;

      B = m_Be(:,:,iPG);

      eps_fluct(:,iPG) = B*u;

      % Deformacion aplicada a cada punto de Gauss
      eps_new(:,iPG) = DefMacro(:,iPG)+eps_fluct(:,iPG);
      %eps_new(:,iPG) = DefMacro+eps_fluct(:,iPG);

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
            error('Matrices Elementales MixStrInj_Quad_q1: Modelo constitutivo no definido.')
      end

      %Se almacena para los tensor tangente constitutivos para realizar homogeneización y análisis de
      %bifurcación
      if esImplex
         % IMPLEX: En el caso de tener activado el implex se almacena las matrices constitutivas tangentes 
         %implícitas para el análisis de bifurcación.
         %Por ahora se asume que el modelo constitutivo siempre devuelve una estructura cuando está
         %seleccionado el implex.
         ctImplexEstab = ct.Implex*m_FactEstab(iPG);
         m_TensorTang(:,:,iPG) = ctImplexEstab;
         %Se almacena los tensores tangentes constitutivo implícitos para análisis de bifurcación como si
         %fuera PG adicionales, tantos como nPG. Se almacena en los índices (:,:,nPG+1:2*nPG).
         %Hay que tener cuidado cuando se hace un loop por la componente de PG sea hasta nPG, y no usando la
         %tercera dimensión de m_TensorTang (size(m_TensorTang,3)).
         %m_TensorTang(:,:,iPG+nPG) = ct.Impli*m_FactEstab(iPG);
         m_TensorTang(:,:,iPG+nPG) = ct.Impli;
         %En los cálculos para el ensamblaje se utiliza el implex.
         ctEstab = ctImplexEstab;
      else
         ctEstab = ct*m_FactEstab(iPG);
         %m_TensorTang(:,:,iPG) = ctEstab;
         %Para el análisis de bifurcación se guarda los tensores constitutivos sin estabilizar (se usa los de
         %los 4 primeros PG, no central, para la homgeneización).
         m_TensorTang(:,:,iPG) = ct;
      end
      
      m_DSigma = sigma_new(:,iPG)-sigma_old(:,iPG);
      m_DSigmaEstab(:,iPG) = m_DSigma*m_FactEstab(iPG);

      % Cálculo de fint
      fint = fint+B'*m_DSigmaEstab(:,iPG)*m_pesoPG(iPG);

      % Cálculo de matriz elemental
      kt = kt+B'*ctEstab*B*m_pesoPG(iPG);

   end
   
   % Actualización de las tensiones estabilizadas.
   %Notar que en m_DSigma viene el valor del último punto de gauss (5), que es el término constante que se
   %suma a las tensiones (ver que m_DSigmaEstab(:,5) = m_DSigma5*(1-factEstab))
   %m_DSigmaEstab = m_DSigmaEstab+m_DSigma*(1-m_FactEstab)';
   m_DSigmaEstab(:,1:4) = bsxfun(@plus,m_DSigmaEstab(:,1:4),m_DSigmaEstab(:,5));
   %En realidad el valor del 5to PG no se utiliza, excepto para imprimir, ya que se integra los valores
   %estabilizados en los 4 primeros PGs.
   m_DSigmaEstab(:,5) = m_DSigma;
   m_SigmaEstabNew = m_SigmaEstabOld+m_DSigmaEstab;
   
   % Tensor tangente estabilizado
   %Notar que en ct viene el correspondiente al último de gauss (5), que es el término constante del campo de
   %tensores tangentes.
   if esImplex
      m_TensorTang(:,:,1:4) = bsxfun(@plus,m_TensorTang(:,:,1:4),ctEstab);
      m_TensorTang(:,:,5) = ct.Implex;
%       m_TensorTang(:,:,6:9) = bsxfun(@plus,m_TensorTang(:,:,6:9),m_TensorTang(:,:,10));
%       m_TensorTang(:,:,10) = ct.Impli;
   else
%       m_TensorTang(:,:,1:4) = bsxfun(@plus,m_TensorTang(:,:,1:4),ctEstab);
%       m_TensorTang(:,:,5) = ct;
   end

   % Se ordena las matrices como vectores columnas
   sigma_new = sigma_new(:);
   hvar_new = hvar_new(:);
   aux_var = aux_var(:);
   eps_new = eps_new(:);
   eps_fluct = eps_fluct(:);
   m_VarHistElemNew = m_SigmaEstabNew(:);

end