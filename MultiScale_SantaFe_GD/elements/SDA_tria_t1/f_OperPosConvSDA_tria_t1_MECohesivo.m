function [sigma_new,eps_new,eps_fluct_new,hvar_new,m_VarHistElemNew,e_VarAuxPG,m_VarAuxElem] = ...
   f_OperPosConvSDA_tria_t1_MECohesivo(...
   sigma_new,eps_new,eps_fluct_new,hvar_new,hvar_old,m_VarHistElemNew,e_VarAuxPG,m_VarAuxElem,m_CT,...
   nElem,m_BT,e_DatMatSet,e_DatElemSet,m_NumElem,m_ElemPGImpr,e_VG)

   %Se considera que este elemento suele puede tener 2 puntos de gauss, uno para la parte
   %regular y otra para la singular.
   ntens = e_VG.ntens;
   ndime = e_VG.ndime;
   nomArchCompl = e_VG.fileCompleto;
   sihvarpg = e_DatMatSet.sihvarpg;
   siavarpg = e_DatMatSet.siavarpg;
   esImplex = e_DatMatSet.esImplex;
   iStep = e_VG.istep;
   %Variables micro
   e_VGMicro     = e_DatMatSet.e_VG;
   e_DatSetMicro = e_DatMatSet.e_DatSet;
   xx            = e_DatMatSet.xx;
   
   %parfor iElem = 1:nElem
   for iElem = 1:nElem
      
      condBif = m_VarAuxElem(1,iElem);     
      
      % Operaciones posterior a la convergencia de los problemas micro.
      %Considerando que previo a la bifurcación solo interesa el PG regular, se realiza los cambios en este
      %punto de gauss únicamente.
      if condBif>1
         nPG = 2;
      else
         nPG = 1;
      end
      for iPG = 1:nPG
         %No se puede realizar este llamado a este nivel ya que no está disponible c_CT de los elementos
         %micro, ya que fueron descartados para no transferir tanta información (Se podría determinar c_CT en
         %este punto del programa).
%          [e_VarEst_new,e_VarAux] = f_OperPosConv(e_VarEst_new,e_VarEst_old,e_VarAux,e_DatSet,c_CT,e_VGMicro);
         %Por ahora se asume que en el llamado de la función a nivel micro no se necesita del c_CT, por
         %ejemplo que el instante de bifurcación fue determinado por el valor de fload.
         c_CTMicro = cell(e_VGMicro.nSet,1);
         %Se asume como generalmente va pasar en los modelos multiescala, al ser estructuras, que se tiene una
         %componente de las variables históricas por PG.
         e_VarEst_newMicro = hvar_new(iPG,iElem).e_VarEst;
         e_VarEst_oldMicro = hvar_old(iPG,iElem).e_VarEst;
         e_VarAuxMicro = hvar_new(iPG,iElem).e_VarAux;
         %
         [e_DatSetMicro,hvar_new(iPG,iElem).e_VarEst,hvar_new(iPG,iElem).e_VarAux,...
             e_VGMicro] =...
          f_OperPosConv(hvar_new(iPG,iElem).u , xx ,...
            e_VarEst_newMicro,e_VarEst_oldMicro,e_VarAuxMicro,e_DatSetMicro,c_CTMicro,e_VGMicro);
      end
      
      %Se guarda los elementos del dominio localizado obtenidos por el criterio de incremento de deformación
      %para impresión (ver si no hacerlo en el postproceso de multiescala, para los PG que se imprimen).
      %Hacerlo fuera de los condicionales permite ver la evolución del área que localiza según este criterio, .
      %Notar que se utiliza como deformación sobre la que se proyecta la macro del PG regular, pero en
      %realidad debería ser BetaxNormal.
%       e_VarAuxPG(1,iElem).m_ElemLocCalc = f_DominioLoc(hvar_new(1,iElem).e_VarEst,...
%          hvar_old(1,iElem).e_VarEst,eps_new(1:ntens,iElem),e_DatSetMicro,e_VGMicro);
      %
      % Operaciones que se realiza hasta que se detecta la bifurcación.       
      if ~condBif
         %Mientras que no bifurca se calcula el dominio localizado proyectando el incremento de deformación
         %fluctuante sobre la deformación macro del punto de gauss regular.
         [e_VarAuxPG(1,iElem).m_ElemLocCalc,e_VarAuxPG(1,iElem).c_ProyIncrEps] = f_DominioLoc(...
            hvar_new(1,iElem).e_VarEst,hvar_old(1,iElem).e_VarEst,eps_new(1:ntens,iElem),e_DatSetMicro,...
            e_VGMicro);
         %
         %Para el análisis de bifurcación se toma el tensor de punto de gauss regular (justamente 
         %es el único que interesa previo a la misma).
         if esImplex
            %Tensores tangentes constitutivos para bifurcación (se toma el tensor implícito).
            m_CTBif = m_CT(:,:,3,iElem);
         else
            m_CTBif = m_CT(:,:,1,iElem);
         end
         if any(m_ElemPGImpr(1,m_ElemPGImpr(2,:)==1)==iElem)
            %Para impresión de det(Q).
            %Se cambia el nombre del archivo para individualizar el elemento y el punto de Gauss 
            %que corresponde a la celda unitaria.
            nombrArchBif = [nomArchCompl,'_E',int2str(m_NumElem(iElem)),'.analisisBif'];
            [condBif,m_angBif] = f_CondBifct(m_CTBif,nombrArchBif,e_VG);
         else
            [condBif,m_angBif] = f_CondBifct(m_CTBif,e_VG);
         end
         %
         if condBif==1
            %Impresión en pantalla
            fprintf('** Se detectó bifurcación en el elemento %d.\n',m_NumElem(iElem))
            %
            if ~isfield(e_DatMatSet,'angNormalImp')
               m_vecBif = [cos(m_angBif);sin(m_angBif)];
               %Se busca la dirección más paralela.
               %m_vecDirPref = [cos(m_angBifPref);sin(m_angBifPref)];
               %[~,ind] = max(abs(m_vecBif'*m_vecDirPref));
               %Se busca la dirección más perpendicular.
               if isfield(e_DatMatSet,'angBanda')
                  angBanda = e_DatMatSet.angBanda;
                  m_vecBand = [cos(angBanda);sin(angBanda)];
                  [~,ind] = min(abs(m_vecBif'*m_vecBand));
               else
                  %Elegir otro criterio para elegir el ángulo del análisis que se debe adoptar.
                  angBanda = pi/2;
                  m_vecBand = [cos(angBanda);sin(angBanda)];
                  [~,ind] = min(abs(m_vecBif'*m_vecBand));
               end
               angBif = m_angBif(ind);
               m_Normal = m_vecBif(:,ind);
               if length(m_angBif)==1
                  angBetaBif = angBif;
               else
                  angBetaBif = m_angBif(1:end~=ind);
               end
            else
               %Esto falla porque falta definir cuál es el angBetaBif. Corregir esto!!!
               angBif = e_DatMatSet.angNormalImp;
               m_Normal = [cos(angBif);sin(angBif)];
               %Implementar que se pueda imponer el ángulo de beta de bifurcación.
               angBetaBif = 0;
               warning(['Como no está definido el ángulo beta al imponer el ángulo de la normal, ',...
                  'se impone como ángulo de beta %f'],angBetaBif)
            end          
            %
            if isfield(e_DatElemSet,'m_NodSolitElem')&&any(e_DatElemSet.m_NodSolitElem(iElem,:))
               %Se toma las matrices de deformación del primer punto de gauss.
               B = m_BT(:,:,1,iElem);
               dN_xy = zeros(e_VG.ndime,e_DatElemSet.npe);
               dN_xy(1,:) = B(1,1:e_VG.ndn:e_DatElemSet.dofpe);
               dN_xy(2,:) = B(2,2:e_VG.ndn:e_DatElemSet.dofpe);
               m_GradPhi = dN_xy(:,e_DatElemSet.m_NodSolitElem(iElem,:));
               %La normal de la fisura se cambia de sentido para apuntar a ese nodo.
               if m_Normal'*m_GradPhi<0
                  m_Normal = -m_Normal;
                  %Como ahora se guarda el ángulo de bifurcación también se lo invierte.
                  angBif = angBif+pi;
               end
            else
               %Se toma el dominio omega_phi+ como el que que contiene el nodo único del elemento triangular.
               [m_Normal,m_GradPhi] = get_solitary_node(m_Normal,m_BT(:,:,1,iElem),e_DatElemSet,e_VG);
            end
            %
            %Para asegurar un modo de apertura de carga en la deformación para determinar el dominio
            %localizado, se invierte el signo del vector de beta de bifurcación para que así se produzca.
            m_BetaBif = [cos(angBetaBif);sin(angBetaBif)];
            if m_Normal'*m_BetaBif<0
               angBetaBif = angBetaBif+pi;
            end
            %
            %Determinación de matriz de proyección tipo beta x normal de la deformación regular.
            %Se determina el vector t perpendicular al vector normal. No interesa la dirección ya que al hacer
            %t x t cualquiera de los dos sentidos da el mismo signo.
            %Está pensado para el problema 2D.
%             m_vecPerNor = [-m_Normal(2);m_Normal(1)];
%             m_TensProy = e_VG.FOAT1*f_Vec2Voigt2D(m_vecPerNor,e_VG)*m_vecPerNor;
%             m_TensProy = e_VG.FOSIT-m_TensProy*m_TensProy';
            %
            %Almacenamiento de variables
            %Se guarda en notación de Voigt el vector normal.
            m_VarAuxElem(1,iElem) = condBif;
            m_VarAuxElem(2:9,iElem) = reshape(f_Vec2Voigt2D(m_Normal,e_VG),[],1);
            m_VarAuxElem(10:17,iElem) = reshape(f_Vec2Voigt2D(m_GradPhi,e_VG),[],1);
            pasoBif = iStep;
            m_VarAuxElem(18,iElem) = pasoBif;
            m_VarAuxElem(22,iElem) = angBif;
            m_VarAuxElem(23,iElem) = angBetaBif;
%             hvar_new(1,iElem).m_TensProy = m_TensProy;
            %
            %Impresión en archivo
            %(se está asumiendo un solo PG regular previo a la bifurcación)
            nomArch = [e_VG.fileCompleto,'.elemBif'];
            fId = fopen(nomArch,'at');
            %#Análisis de bifurcación
            %#| Paso | Elemento | Ángulos de la normal usado (ángulos del análisis de bifurcación)
            if length(m_angBif)==1
               fprintf(fId,['Paso=%d (t=%f), Elemento=%d, Ángulo de la normal ',...
                  'usado=%g° (%g°): Se detecta bifurcación.\n'],iStep,e_VG.Dtime*iStep,...
                  m_NumElem(iElem),angBif/pi*180,m_angBif/pi*180);
            else
               fprintf(fId,['Paso=%d (t=%f), Elemento=%d, Ángulo de la normal ',...
                  'usado=%g° (%g° y %g°): Se detecta bifurcación.\n'],iStep,...
                  e_VG.Dtime*iStep,m_NumElem(iElem),angBif/pi*180,m_angBif/pi*180);
            end            
            fclose(fId);
         end
      end
      %
      %Se separa los condicionales por si se quiere implementar la discontinuidad fuerte desde el paso que se
      %detectó la bifurcación (en el paso siguiente se tiene en cuenta la modificaciones de esta función, ya
      %que se aplican al final del paso).
      %
      condBif = m_VarAuxElem(1,iElem);
      %if condBif==1
      if condBif>0
         if ~isfield(e_DatMatSet,'angNormalImp')
            %Desde que bifurcó hasta que se activa la discontinuidad fuerte se calcula el dominio localizado
            %proyectando el incremento de deformación fluctuante sobre una deformación macro BetaxNormal
            m_Normal = reshape(m_VarAuxElem(2:9,iElem),ntens,ndime);
            angBetaBif = m_VarAuxElem(23,iElem);          
            m_BetaBif = [cos(angBetaBif);sin(angBetaBif)];
            m_DefBNBif = m_Normal*m_BetaBif;
            [e_VarAuxPG(1,iElem).m_ElemLocCalc,e_VarAuxPG(1,iElem).c_ProyIncrEps] = ...
               f_DominioLoc(hvar_new(1,iElem).e_VarEst,hvar_old(1,iElem).e_VarEst,m_DefBNBif,...
               e_DatSetMicro,e_VGMicro);
         else
            %En el caso que se imponga el ángulo de la normal, se sigue calculando el dominio poyectando sobre
            %sobre la deformación (macro) del punto de gauss regular, ya que no se conoce la dirección de
            %beta.
            [e_VarAuxPG(1,iElem).m_ElemLocCalc,e_VarAuxPG(1,iElem).c_ProyIncrEps] = ...
               f_DominioLoc(hvar_new(1,iElem).e_VarEst,hvar_old(1,iElem).e_VarEst,eps_new(1:ntens,iElem),...
               e_DatSetMicro,e_VGMicro);
         end
      end
      %
      % Cambio de condiciones de borde externa después que bifurcó.
      %Se cambia la condición de borde solo en el rango donde condBif=1, ya que como está programado no se
      %conoce el ángulo de la normal de bifurcación antes que bifurque (solo si se la impone manualmente) y
      %esta condición de borde pisa la existente, y por lo tanto se debe imponer antes de la segunda condición
      %de borde (que se agrega a la existente) para no perderla.
      condBif = m_VarAuxElem(1,iElem);
      pasoBif = m_VarAuxElem(18,iElem);
      noActBif = m_VarAuxElem(21,iElem);
      if condBif==1&&isfield(e_DatMatSet,'pasosPosBifCBExt')&&...
            iStep==pasoBif+e_DatMatSet.pasosPosBifCBExt&&~noActBif
         fprintf('** Se cambia la CB externa de la CU de ambos PG del elemento %d.\n',m_NumElem(iElem))
         angBif = m_VarAuxElem(22,iElem);
         angTang = angBif+pi/2;
         m_CondBord = e_DatMatSet.f_CambioCBExt(angTang,e_DatMatSet.xx);
         e_VGMicroAux = e_VGMicro;
         e_VGMicroAux.m_CondBord = m_CondBord;
         [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicroAux,e_DatMatSet.xx,...
            e_DatSetMicro,e_VGMicroAux.m_ConecFront);
         %Se cambia en en el punto de gauss regular ya que el punto de gauss singular se ignora hasta que se
         %activa la bifurcación. Luego los datos son copiados al singular cuando se activa la discontinuidad
         %fuerte. Habría que ver que hacer en la descarga del punto de gauss regular, por ahora se recupera
         %las originales.
         hvar_new(1,iElem).m_LinCond = m_LinCondMicro;
         hvar_new(1,iElem).doff = doffMicro;
         hvar_new(1,iElem).dofl = doflMicro;
      end
      %
      % Operaciones que ocurren una vez cuando se activa el SDA, en el mismo paso o ciertos pasos
      % después que se detecta la bifurcación.
      condBif = m_VarAuxElem(1,iElem);
      pasoBif = m_VarAuxElem(18,iElem);
      noActBif = m_VarAuxElem(21,iElem);
      if condBif==1&&iStep==pasoBif+e_DatMatSet.pasosPosBif&&~noActBif
         %Impresión en pantalla
         fprintf('** Se activa discontinuidad fuerte en el elemento %d.\n',m_NumElem(iElem))
         %Impresión en archivo
         %(se está asumiendo un solo PG, el regular, previo a la bifurcación)
         nomArch = [e_VG.fileCompleto,'.elemBif'];
         fId = fopen(nomArch,'at');
         fprintf(fId,'Fin de Paso=%d (t=%f), Elemento=%d: Se activa discontinuidad fuerte.\n',...
            iStep,e_VG.Dtime*iStep,m_NumElem(iElem));
         fclose(fId);
         %
         %Se cambia el tensor de proyección para que no lo siga aplicando en ninguno de los dos puntos de
         %gauss después de que se activó la discontinuidad fuerte.
%          hvar_new(1,iElem).m_TensProy = e_VG.FOSIT;
         %
         condBif = 2;
         m_VarAuxElem(1,iElem) = condBif;
         %Se copia las variables históricas del PG regular al singular.
         sigma_new(ntens+1:2*ntens,iElem) = sigma_new(1:ntens,iElem);
         eps_new(ntens+1:2*ntens,iElem) = eps_new(1:ntens,iElem);
         eps_fluct_new(ntens+1:2*ntens,iElem) = eps_fluct_new(1:ntens,iElem);
         hvar_new(sihvarpg+1:2*sihvarpg,iElem) = hvar_new(1:sihvarpg,iElem);
         %También se transfiere las variables auxiliares
         e_VarAuxPG(siavarpg+1,iElem) = e_VarAuxPG(1,iElem);
         %
         if isfield(e_DatMatSet,'m_ElemLocImp')
            %Se guarda el punto de gauss regular, ya que después se realiza la copia en el punto de gauss
            %singular cuando se impone la discontinuidad fuerte.
            m_ElemLoc = e_DatMatSet.m_ElemLocImp;
         else
            %Función de determinación del dominio localizado.
%             error(['Operaciones de fin de paso: Triángulo SDA lineal del modelo cohesivo: ',...
%                'Determinación del dominio localizado no implementado.'])
            m_ElemLoc = e_VarAuxPG(1,iElem).m_ElemLocCalc;
         end
         %Área del dominio localizado
         omegaMicroL = sum(e_DatMatSet.m_VolElem(m_ElemLoc));
         %Longitudes y normales micro.
         if isfield(e_DatMatSet,'lMicroImp')
            lMicro = e_DatMatSet.lMicroImp;
            m_Normal = reshape(m_VarAuxElem(2:9,iElem),ntens,ndime);
            if isnumeric(lMicro)&&~isnan(lMicro)
               c_NormalesMicro = cell(e_VGMicro.nSet,2);
               %Normales micro iguales para todos los elementos
               c_NormalesMicro(:,1) = arrayfun(@(x)repmat(m_Normal,x.nElem),e_DatSetMicro,...
                  'UniformOutput',false);
               %Longitudes micro iguales para todos los elementos
               c_NormalesMicro(:,2) = arrayfun(@(x)lMicro*ones(1,x.nElem),e_DatSetMicro,...
                  'UniformOutput',false);
            else               
               % Lectura de un archivo las normales micro y longitudes micro               
               c_NormalesMicro = f_NormalesMicro(m_ElemLoc,hvar_new(2,iElem).e_VarEst,...
                  e_DatSetMicro,m_Normal,e_VGMicro);
               %Por ahora hasta que se programe bien el tensor tangente del elemento del SDA, se utiliza como
               %ancho localizado para el cálculo del mismo, una media volumétrica de las medias indicada en el
               %archivo.
               c_lMicro = c_NormalesMicro(:,2);
               m_VolElemMicro = e_DatMatSet.m_VolElem;
               %Recordar que viene inf en el caso no está definido la longitud micro. Por eso, en el caso de
               %no pertenecer al dominio localizado, se la excluye y si está en el dominio localizado
               %el producto del área por la longitud micro va devolver Inf, indicando que falta ingresar ese
               %valor.
               m_lMicro = zeros(1,e_VGMicro.nElem);
               m_lMicro([e_DatSetMicro.m_IndElemSet]) = [c_lMicro{:}];
               %No se usa f_HomogArea para no tener que repetir los valores en los 4PG del elemento y ya se
               %tiene una matriz con las áreas.
               lMicro = m_VolElemMicro(m_ElemLoc)*m_lMicro(m_ElemLoc)'/omegaMicroL;
            end
         else
            %Acá debe colocarse la función de determinación de la longitud micro.
            error(['Operaciones de fin de paso: Triángulo SDA lineal del modelo cohesivo: ',...
               'Determinación de la longitud micro no implementada.'])

         end
         %
         if isfield(e_DatMatSet,'longFisImp')
            longFis = e_DatMatSet.longFisImp;
         else
            %Acá debe colocarse la función de determinación de la longitud de fisura.
            error(['Operaciones de fin de paso: Triángulo SDA lineal del modelo cohesivo: ',...
               'Determinación de la longitud de fisura no implementada.'])
         end
         if isfield(e_DatMatSet,'facNormMicroImp')
            facNormMicro = e_DatMatSet.facNormMicroImp;
         else
            %Acá debe colocarse la función de determinación de la media de las normales micro a largo de la
            %fisura (eje central del dominio localizado).
            error(['Operaciones de fin de paso: Triángulo SDA lineal del modelo cohesivo: ',...
               'Determinación de la media de las normales micro.'])
         end
         %
         lMacro = lMicro*e_DatMatSet.omegaMicro/omegaMicroL;
         %
         %Como la matriz m_CondBord no se guarda (ver si no guardarlas en las variables históricas), por
         %ahora, cuando se cambia las CB externas se la vuelve a calcular para tenerlas.
         if isfield(e_DatMatSet,'pasosPosBifCBExt')&&...
               e_DatMatSet.pasosPosBifCBExt>=0&&pasoBif+e_DatMatSet.pasosPosBifCBExt<=iStep
            %Primero se recupera la condición del punto de gauss regular a la original que tenía para la
            %descarga, ya que desde que se detecta la bifurcación hasta que se impone la discontinuidad fuerte
            %y se activa el modelo cohesivo se cambió la condición de borde del punto de gauss regular por que
            %como está programado es el único que se resuelve previo a la discontinuidad fuerte.
            %m_CondBordPGR = e_VGMicro.m_CondBord;
            [m_LinCondMicroPGR,~,~,doffMicroPGR,doflMicroPGR] = f_CondBord(e_VGMicro,...
               e_DatMatSet.xx,e_DatSetMicro,e_VGMicro.m_ConecFront);
            hvar_new(1,iElem).m_LinCond = m_LinCondMicroPGR;
            hvar_new(1,iElem).doff = doffMicroPGR;
            hvar_new(1,iElem).dofl = doflMicroPGR;
            %
            angBif = m_VarAuxElem(22,iElem);
            angTang = angBif+pi/2;
            m_CondBord = e_DatMatSet.f_CambioCBExt(angTang,e_DatMatSet.xx);
         else
            %En el caso de no haberse cambiado las condiciones de borde de partida, se utiliza las mismas.
            m_CondBord = e_VGMicro.m_CondBord;
         end
         %
         m_ConecFrontL = f_LadFront(find(m_ElemLoc),e_DatMatSet.m_ElemVec,e_DatSetMicro,e_VGMicro);
         e_VGMicroAux = e_VGMicro;
         switch e_DatMatSet.tipoCBL
            case 0  %No se impone CB en el dominio localizado
               %Se copia las mismas matrices que el punto regular (cuidado acá si se cambia la CB entre la
               %bifurcación y la discontinuidad fuerte, que ambos PG recupera la condiciones de borde
               %originales, es decir las previas a ese cambio).
               m_LinCondMicro = hvar_new(1,iElem).m_LinCond;
               doffMicro = hvar_new(1,iElem).doff;
               doflMicro = hvar_new(1,iElem).dofl;
            case 1  %CB de mínima restricción
               c_ConecFront = {e_VGMicro.m_ConecFront,m_ConecFrontL};
               m_CteMinRestrL = [0;0;0]*omegaMicroL;
               e_VGMicroAux.m_CondBord = m_CondBord; 
               [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicroAux,e_DatMatSet.xx,...
                  e_DatSetMicro,c_ConecFront,m_CteMinRestrL);
            case 3  %CB lineal
               m_CondBord = f_CBLineal(m_CondBord,m_ConecFrontL,e_VGMicro.ndn);
               e_VGMicroAux.m_CondBord = m_CondBord; 
               [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicroAux,e_DatMatSet.xx,...
                  e_DatSetMicro,e_VGMicroAux.m_ConecFront);
            case 4  %CB Taylor
               m_CondBord = f_CBTaylor(m_ElemLoc,e_DatSetMicro,m_CondBord,hvar_new(2,iElem).doff,e_VGMicro);
               e_VGMicroAux.m_CondBord = m_CondBord; 
               [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicroAux,e_DatMatSet.xx,...
                  e_DatSetMicro,e_VGMicroAux.m_ConecFront);
            otherwise
               error(['Operaciones de fin de paso: Triángulo SDA lineal del modelo cohesivo: ',...
                  'Condición de borde del dominio localizado no implementado'])
         end
         %Inicialización de la tracción
         %Se considera como tracción inicial la del PG regular multiplicada por la normal macro.
         m_TraccionInicial = m_Normal'*sigma_new(1:ntens,iElem);
         %
         % Discontinuidad de las tracciones
         [m_difTracc,m_TracHomogL,m_TracHomog,m_NormalMicroMedia] = f_DiscTracc(...
            c_NormalesMicro,m_Normal,hvar_new(2,iElem).e_VarEst,...
            sigma_new(1:ntens,iElem),omegaMicroL,longFis,m_ElemLoc,e_DatSetMicro,e_VGMicro);
         nomArch = [e_VG.fileCompleto,'.discTrac'];
         fId = fopen(nomArch,'wt');
         fprintf(fId,'** Elemento %d.\n',m_NumElem(iElem));
         fprintf(fId,'(en el momento de activación del modelo cohesivo multiescala).\n');
         fprintf(fId,'* Diferencia de tracción macro: DTx = %g y DTy = %g (módulo %g).\n',m_difTracc(1),...
            m_difTracc(2),norm(m_difTracc));
         fprintf(fId,'* Tracción macro homogeneizada en el dominio localizado: TLx = %g y TLy = %g (módulo %g).\n',...
            m_TracHomogL(1),m_TracHomogL(2),norm(m_TracHomogL));
         fprintf(fId,'* Tracción macro homogeneizada en el todo el dominio micro: Tx = %g y Ty = %g (módulo %g).\n',...
            m_TracHomog(1),m_TracHomog(2),norm(m_TracHomog));
         fprintf(fId,'* Normal micro homogeneizada en el dominio localizado: Nhx = %g y Nhy = %g (ángulo %g°).\n',...
            m_NormalMicroMedia(1),m_NormalMicroMedia(2),atan(m_NormalMicroMedia(2)/m_NormalMicroMedia(1))*180/pi);
         fclose(fId);
         %
         hvar_new(2,iElem).m_ElemLoc = m_ElemLoc;
         hvar_new(2,iElem).omegaMicroL = omegaMicroL;
         hvar_new(2,iElem).lMicro = lMicro;
         hvar_new(2,iElem).lMacro = lMacro;            
         hvar_new(2,iElem).m_LinCond = m_LinCondMicro;
         hvar_new(2,iElem).doff = doffMicro;
         hvar_new(2,iElem).dofl = doflMicro;
         hvar_new(2,iElem).c_NormalesMicro = c_NormalesMicro;
         hvar_new(2,iElem).longFis = longFis;
         hvar_new(2,iElem).facNormMicro = facNormMicro;
         %
         m_VarHistElemNew(ndime+1:2*ndime,iElem) = m_TraccionInicial;
      end
      %
      % Operaciones que ocurren desde que se impone el SDA.
      condBif = m_VarAuxElem(1,iElem);
      if condBif>1
         %Se cambia la condición de bifurcación a 3 cuando ya ha pasado un paso desde que se activó la SD y el
         %modelo cohesivo. Esto se utiliza si se quiere realizar una operación sólo en un cierta cantidad de
         %después que se activa.
%          if condBif==2&&iStep==pasoBif+e_DatMatSet.pasosPosBif+1
%             m_VarAuxElem(1,iElem) = condBif+1;
%          end
%          if condBif==2&&iStep==pasoBif+e_DatMatSet.pasosPosBif+10
%             m_VarAuxElem(1,iElem) = condBif+1;
%          end
      end
      %
   end
   
end

function [m_difTracc,m_TracHomogL,m_TracHomog,m_NormalMicroMedia] = f_DiscTracc(...
   c_NormalesMicro,m_NormalMacro,e_VarEst,sigmaMacro,omegaMicroL,longFis,m_ElemLoc,e_DatSet,e_VG)

   %%Discontinuidad de tracciones
   %Determina las tensiones instántaneas, y con ellas las tracciones instantáneas, en el dominio micro total y
   %localizado para ver la discontinuidad de las tracciones. En el dominio localizado se utiliza las normales
   %micro, mientras que en el dominio total se utiliza la normal macro. En este último caso directamente se
   %usa la tensión homogeneizada del punto de gauss regular macro (se asume que en el caso en el elemento
   %tipo, de inyección de deformación, se envía la tensión homogeneizada a partir de las estabilizadas micro).
   %Tracciones micro         
   nSet = e_VG.nSet;
   nTens = e_VG.ntens;
   nDime = e_VG.ndime;
   c_TraccL = cell(nSet,1);
   for iSet = 1:e_VG.nSet
      nElem = e_DatSet(iSet).nElem;
      nPG = e_DatSet(iSet).e_DatElem.npg;
      eltype = e_DatSet(iSet).e_DatElem.eltype;
      m_NormalesMicro = c_NormalesMicro{iSet,1};
      m_lMicro = c_NormalesMicro{iSet,2};
      m_TraccL = zeros(nDime,nPG,nElem);
      switch eltype
         case {2,4,8}
            m_Sigma = e_VarEst(iSet).sigma;
         case 20
            %En el elemento mixto se debe utilizar las tensiones micro estabilizadas.
            m_Sigma = e_VarEst(iSet).VarHistElem;
         otherwise
            error(['Modelo Multiescala: Salto de Tracciones: Tensiones de homogeneización: Modelo ',...
               'constitutivo no definido.'])
      end
      m_Sigma = reshape(m_Sigma,nTens,nPG,nElem);
      for iElem = 1:nElem
         m_NormalMicroElem = m_NormalesMicro(:,:,iElem);
         m_lMicroElem = m_lMicro(iElem);
         for iPG = 1:nPG
            %m_TraccL(:,iPG,iElem) = m_NormalesMicro(:,:,iPG,iElem)'*m_Sigma(:,iPG,iElem);
            m_TraccL(:,iPG,iElem) = m_NormalMicroElem'*m_Sigma(:,iPG,iElem)/m_lMicroElem;
         end
      end
      c_TraccL{iSet} = m_TraccL;
   end
   %No sería necesario realizar esto ya que las normales micro son nulas cuando se está fuera del
   %dominio localizado.
   c_DetJTLoc = arrayfun(@(x)bsxfun(@times,x.m_DetJT,m_ElemLoc(x.m_IndElemSet)),e_DatSet,...
      'UniformOutput',false);
   %m_TracHomogL = f_HomogArea(c_TraccL,nDime,omegaMicroL,c_DetJTLoc,e_DatSet,e_VG);
   m_TracHomogL = f_HomogArea(c_TraccL,nDime,longFis,c_DetJTLoc,e_DatSet,e_VG);
   m_TracHomog = m_NormalMacro'*sigmaMacro;
   %Diferencia de tracciones          
   m_difTracc = m_TracHomog-m_TracHomogL;
   
   %%Normal media
   %Se calcula una normal media de las normales micro seleccionadas en el dominio localizado.
   %Se transforma el vector en notación de voigt en notación normal.
   %c_NormalesMicro = cellfun(@(x)[x(1,1,:,:);x(2,2,:,:)],c_NormalesMicro(:,1),'UniformOutput',false);
   c_NormalesMicro = cellfun(@(x)[x(1,1,:);x(2,2,:)],c_NormalesMicro(:,1),'UniformOutput',false);
   %Cómo ya no se lleva las normales micro por PG, sino por elemento, se repite las normales en los 4 PG para
   %poder integrar.
   c_NormalesMicro = arrayfun(@(x,y)repmat(x{1},[1,y.e_DatElem.npg,1]),c_NormalesMicro,e_DatSet,...
      'UniformOutput',false);
   m_NormalMicroMedia = f_HomogArea(c_NormalesMicro,nDime,omegaMicroL,c_DetJTLoc,e_DatSet,e_VG);
   m_NormalMicroMedia = m_NormalMicroMedia/norm(m_NormalMicroMedia);
   
end