function [hvar_new,m_VarHistElemNew,e_VarAuxPG,m_VarAuxElem] = ...
   f_OperPosConvSDA_QSFe(uElem,xx,conec,EF_sidesCutCPF,...
   hvar_new,hvar_old,m_VarHistElemNew,...
   m_VarHistElemOld, e_VarAuxPG,m_VarAuxElem,m_CT,...
   nElem,e_DatMatSet,e_DatElemSet,e_DatSet,e_VG)

%ntens = e_VG.ntens;
%ndime = e_VG.ndime;
%nomArchCompl = e_VG.fileCompleto;
%
%sihvarpg = e_DatMatSet.sihvarpg;
%siavarpg = e_DatMatSet.siavarpg;
esImplex = e_DatMatSet.esImplex;
%m_ElemPGImpr = e_DatMatSet.m_ElemPGImpr;
iStep = e_VG.istep;
%Variables micro
e_VGMicro       = e_DatMatSet.e_VG;
e_DatSetMicro   = e_DatMatSet.e_DatSet;
m_NumElem       = e_DatSet.m_NumElem;

p_indActSTmacro =  e_DatElemSet.pointersVHE.p_indActSTmacro    ;
p_condBif       =  e_DatElemSet.pointersVAE.p_condBif ;
p_leq_elem      =  e_DatElemSet.pointersVAE.p_leq_elem ;
%p_normal_bif    =  e_DatElemSet.pointersVAE.p_normal_bif ;
pAngBif = e_DatElemSet.pointersVAE.pAngBif;
p_n_tens        =  e_DatElemSet.pointersVAE.p_n_tens    ;
p_ref_vector    =  e_DatElemSet.pointersVAE.p_ref_vector ;
p_CeBif = e_DatElemSet.pointersVAE.pCeBif;
i_vectVHElem = e_DatElemSet.pointersVHE.i_vectVHElem;
i_indST = e_DatElemSet.pointersVHE.i_indST;

dN_xy           = e_DatSet.dN_xy  ;

%Nombre de archivo que se guarda los datos de bifurcación.
nomArchBif = [e_VG.fileCompleto,'.elemBif'];
%Se abre el archivo.
fId = fopen(nomArchBif,'at');
%La mayoría de los cálculos y cambios se realiza en el PG 6, que es el que está llevando de la deformación
%localizada.
iPG = 6;

%parfor iElem = 1:nElem
for iElem = 1:nElem
   
   condBif           =  m_VarAuxElem(p_condBif,iElem);
   %leq_elem          =  m_VarAuxElem(p_leq_elem  ,iElem) ;
   %n_bifurc          =  m_VarAuxElem(p_normal_bif,iElem) ;
   %n_bifurc          =  reshape(n_bifurc,2,2);
   ind_ActState_new  =  m_VarHistElemNew (p_indActSTmacro   ,iElem) ;
   ind_ActState_old  =  m_VarHistElemOld (p_indActSTmacro   ,iElem) ;
   %ind_state_old =  m_VarHistElemOld(i_indST,iElem);
   ind_state_new =  m_VarHistElemNew(i_indST,iElem);
   
   e_VG.iElemSet = iElem;
   e_VG.iElemNum = m_NumElem(iElem);
   
   %VER SI ESTAS VARIABLES NO PASARLAS COMO VARIABLES HISTÓRICAS.
   %Para realizar operaciones una sola vez después que bifurcó sin guardar en la variables históricas.
   %Notar que esto solo funciona en el caso que se permita hacer operaciones inmediatamente que bifurcó.
   condBifPrev = condBif;
   %Se guarda la normal de bifurcación seleccionada previa (old) para ver si vuelve a cambiar durante el
   %periodo intermedio entre la bifurcación y la activación de la SD, y así realizar ciertas operaciones
   %únicamente en esos casos. También para recuperar los datos de bifurcación previos cuando se detecta 
   %bifurcación.
   n_tensPrev = m_VarAuxElem(p_n_tens,iElem);
   n_tensPrev = [n_tensPrev(1);n_tensPrev(6)];
   m_AngBifPrev = m_VarAuxElem(pAngBif,iElem);
   
   % Operaciones que se realiza hasta que se detecta la bifurcación.
   if ~condBif
 
      % evaluar condicion de bifurcacion en el PG macro    
      if esImplex
         m_CTBif = m_CT(:,:,12,iElem);
      else
         m_CTBif = m_CT(:,:,6,iElem);
      end  
%       [condBif,~,leq_elem,n_bifurc,n_tens] = ...
%          bifurcation_condition_MACRO...
%          (m_CTBif,leq_elem,condBif,e_VG,xx,...
%          conec(iElem,:),uElem(:,iElem),e_DatElemSet,dN_xy(:,:,iElem),e_DatSet);
      [condBif,m_AngBif] = f_CondBifct(m_CTBif,e_VG);
      
      %Impresión del tensor tangente homogeneizado
      if 1
         nomArchCt = [e_VG.fileCompleto,'.Ct'];
         format = ['%d %d',repmat(' %.15g',1,e_VG.ntens^2),'\n'];
         if e_VG.istep==1
            %Para inicializar el archivo cada vez se corre de nuevo, se abre distinto.
            fIdCt = fopen(nomArchCt,'wt');
         else            
            fIdCt = fopen(nomArchCt,'at');
         end
         fprintf(fIdCt,format,iStep,e_VG.iElemNum,m_CTBif);
         fclose(fIdCt);
      end
      %Ver si es necesario seleccionar una normal previo a la bifurcación. Pareciera que se usa en la 
      %SDA_Properties para determinar el GradPhi y los nodos libres, pero esto se usa únicamente cuando se
      %inyecta.
      m_NormBif = [cos(m_AngBif);sin(m_AngBif)];
      n_tens = Normal_vector_selection(xx,conec(iElem,:),uElem(:,iElem),m_NormBif,0,...
            e_VG,e_DatElemSet,e_DatMatSet,dN_xy(:,:,iElem),EF_sidesCutCPF(iElem),...
            hvar_new(:,iElem),hvar_old(:,iElem));
      
      m_VarAuxElem(p_condBif,iElem) = condBif;
      %m_VarAuxElem(p_leq_elem  ,iElem) = leq_elem    ;
      %m_VarAuxElem(p_normal_bif,iElem) = reshape(n_bifurc,4,1);
      %Se almacena los ángulos de bifurcación.
      m_VarAuxElem(pAngBif,iElem) = m_AngBif;
      %m_VarAuxElem(p_n_tens    ,iElem) = reshape(n_tens,8,1);
      m_VarAuxElem(p_n_tens,iElem) = n_tens(:);
    
   end
   %      elseif (condBif == 1) && (indActSTmacro==1)
   %Sebastian: Cambiado esta línea: elseif  (condBif == 1 && indSTmacro <= 1) || (ind_ActState_new==2 && ind_ActState_old==0)
   %else
   %Para realizar los cálculos inmediatamente después que se detecta bifurcación.
   if condBif
      
      %Se guarda como tensor de bifurcación de desacrga elástico como es que corresponde al paso que se
      %detecta bifurcación.
      if ~condBifPrev
         if esImplex
            %El implex ya utiliza el tensor elástico (por el daño) como tensor tangente.
            m_CeBif = m_CT(:,:,6,iElem);
         else
            %Para obtener el tensor elástico
            e_VGAux = e_VG;
            e_VGAux.elast = 1;
            e_VGAux.iPG = 5;
            condBifR = 0;
            %
            %Se impone un deformación arbitraria, no interesa para el cálculo de tensor tangente.
            m_IDefR = [1;0;0;0];
            m_IDefS = [0;0;0;0];
            ksd = 0;
            m_vectVHElem_new = m_VarHistElemNew(i_vectVHElem,iElem);
            m_CeBif = f_RMap_MEStafe(m_IDefR,m_IDefS,hvar_new(5,iElem),e_DatMatSet,condBifR,...
               ind_ActState_old,e_VGAux,m_vectVHElem_new,ksd);            
         end
         m_VarAuxElem(p_CeBif,iElem) = reshape(m_CeBif,e_VG.ntens*e_VG.ntens,[]);
      end
      
      %Se ejecuta entre la bifurcación y la activación de la SD, y mientras no se active la SD (una vez
      %activada la bifurcación por primera vez no se entra más).
      %Cómo esta parte está previo a realizar un cambio de estado en ind_ActState_new (más abajo), se obliga
      %entrar al menos una vez (por ejemplo que se pase directamente de ind_ActState de 0 a 2,
      %ind_ActState_new==2 && ind_ActState_old==0, que no verifica ind_state_new<=1, pero de esta forma sí lo
      %hace).
      if ind_state_new<=1
         
         %Se sigue almacenando el tensor tangente durante la primera vez que se está en estado 1.
         if esImplex
            m_Ce = m_CT(:,:,12,iElem);
         else
            m_Ce = m_CT(:,:,6,iElem);
         end
         nomArchCt = [e_VG.fileCompleto,'.Ct'];
         format = ['%d %d',repmat(' %.15g',1,e_VG.ntens^2),'\n'];
         fIdCt = fopen(nomArchCt,'at');
         fprintf(fIdCt,format,iStep,e_VG.iElemNum,m_Ce);
         fclose(fIdCt);
         
         %Como criterio se utiliza los datos del análisis de bifurcación previo a detectar bifurcación
         %(detQ<0), estos han mostrado ser más cercanos a una normal de bifurcación más coherente con el
         %problema.
         %En los problemas monoescala no se vio lo mismo (VERIFICAR!)
         if ~condBifPrev
            m_VarAuxElem(pAngBif,iElem) = m_AngBifPrev;
            %Esta línea siguiente no es necesario ya que en la función Normal_vector_selection se selecciona
            %una nueva normal con los dos ángulos posibles.
            %m_VarAuxElem(p_n_tens,iElem) = n_tensPrev(:);
         end
         %
         m_AngBif = m_VarAuxElem(pAngBif,iElem)';
         m_NormBif = [cos(m_AngBif);sin(m_AngBif)];

         %Selección de normal de bifurcación         
         [n_tens,m_tens,ref_vector] = Normal_vector_selection(...
            xx,conec(iElem,:),uElem(:,iElem),m_NormBif,0,e_VG,e_DatElemSet,e_DatMatSet,dN_xy(:,:,iElem),...
            EF_sidesCutCPF(iElem),hvar_new(:,iElem),hvar_old(:,iElem));
         
         %Definición del signo del producto (Beta xsym Normal) para proyectar las deformaciones fluctuante
         %micro y así determinar el dominio localizado.
         m_DefBNBif = [n_tens(1,1)*m_tens(1,1);n_tens(2,2)*m_tens(2,2);0;
            (n_tens(1,1)*m_tens(2,2)+n_tens(2,2)*m_tens(1,1))];
         %Para asegurar un modo de apertura de carga en la deformación Beta xsym Normal para determinar 
         %el dominio localizado, se invierte el signo del vector de beta de bifurcación para que así 
         %se produzca (Beta * Normal = traza(Beta xsym Normal) = (Beta xsym Normal)_11+(Beta xsym Normal)_22).
         if m_DefBNBif(1)+m_DefBNBif(2)<0
            m_DefBNBif = -m_DefBNBif;
         end
         %
         %m_DefBNBif = [1;0;0;0];
         %
         m_VarAuxElem(p_n_tens,iElem) = n_tens(:);
         m_VarAuxElem(p_ref_vector,iElem) = ref_vector;

         %Cálculo del dominio de localizacion para impresión
         [e_VarAuxPG(iPG,iElem).m_ElemLocCalc,e_VarAuxPG(iPG,iElem).c_ProyIncrEps] = ...
            f_DominioLoc(hvar_new(iPG,iElem).e_VarEst,hvar_old(iPG,iElem).e_VarEst,m_DefBNBif,...
            e_DatSetMicro,e_VGMicro);
        
         if ~condBifPrev
            %Impresión en pantalla
            fprintf('** Se detectó bifurcación en el elemento %d.\n',m_NumElem(iElem))
            %Impresión en archivo
            fprintf(fId,'Fin de Paso=%d (t=%f), Elemento=%d: Se detecta bifurcación.\n',...
               iStep,e_VG.Dtime*iStep,m_NumElem(iElem));
         end
         
         %Dentro de este condicional se entra únicamente en el primer paso que cambió el ind_ActState de 0
         %a 1 ó 2. Ver que en el condicional superior se entra únicamente cuando se detecta por primera vez
         %bifurcación (es decir cuando es ind_state_old <= 1).
         %También se hace una verificación si la normal de bifurcación seleccionada (entre las dos posibles)
         %cambia durante la etapa intermedia.
         n_tensVec = [n_tens(1);n_tens(6)];
         if ~condBifPrev||norm(n_tensVec-n_tensPrev)>0
            
            %Impresión en archivo
            fprintf(fId,['Fin de Paso=%d (t=%f), Elemento=%d: Ángulo de la normal ',...
                  'seleccionado=%g° (%g° y %g°).\n'],iStep,e_VG.Dtime*iStep,...
                  m_NumElem(iElem),atan2(n_tensVec(2),n_tensVec(1))/pi*180,...
                  atan2(m_NormBif(2,:),m_NormBif(1,:))/pi*180);
            
            if isfield(e_DatMatSet,'f_CambioCBExt')
               fprintf('** Se cambia la CB externa de la CU de todos los PG del elemento %d.\n',...
                  m_NumElem(iElem))
               %angTang = atan2(n_tens(1,1),-n_tens(2,2));
               %n_NormalPer = n_tensVec;
               n_NormalPer = m_VarAuxElem(p_ref_vector,iElem);
               %ang = 10/180*pi;
               %n_NormalPer = [cos(ang);sin(ang)];
               m_CondBord = e_DatMatSet.f_CambioCBExt(e_DatMatSet.xx,e_DatSetMicro,n_NormalPer,...
                  e_VGMicro);
               e_VGMicroAux = e_VGMicro;
               e_VGMicroAux.m_CondBord = m_CondBord;
               [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicroAux,e_DatMatSet.xx,
                  e_DatSetMicro,e_VGMicroAux.m_ConecFront);
               %Se cambia las condiciones de borde de todos los puntos de gauss, ya que PG=5 copia del PG=6
               %cuando se activa la fisura.
               %Sobre los PGRs se mantiene las condiciones de borde periódicas después de la activación de
               %la fisura (descarga elástica), como también el PGS (ablandamiento del modelo cohesivo).
               [hvar_new(1:6,iElem).m_LinCond] = deal(m_LinCondMicro);
               [hvar_new(1:6,iElem).doff] = deal(doffMicro);
               [hvar_new(1:6,iElem).dofl] = deal(doflMicro);
            end
            
            % Calculo longitud de Oliver (Oliver X./89)
            m_CoordElem = f_CoordElem(xx,conec(iElem,:));
            npe = e_DatElemSet.npe ;
            new_coord_n_MACRO = zeros(size(m_CoordElem));
            xcg = sum(m_CoordElem(1,:))/npe;
            ycg = sum(m_CoordElem(2,:))/npe;
            bif_angle_rad = atan(n_tensVec(2)/n_tensVec(1));
            %
            cosAngBif = cos(bif_angle_rad);
            senAngBif = sin(bif_angle_rad);
            for inode = 1:npe
               new_coord_n_MACRO(1,inode) = (m_CoordElem(1,inode)-xcg)*cosAngBif+...
                  (m_CoordElem(2,inode)-ycg)*senAngBif;
               new_coord_n_MACRO(2,inode) = -(m_CoordElem(1,inode)-xcg)*senAngBif+...
                  (m_CoordElem(2,inode)-ycg)*cosAngBif;
            end
            %
            % Calculo longitud equivalente para los elementos macro (Oliver/89)
            %switch eltype
               %case {4,10,21,22,23} % Elementos cuadrilaterales - indistinto de si es enriquecido o no.
                  leq_elem = le_quad_q1_epd(m_CoordElem,new_coord_n_MACRO,bif_angle_rad);
                  %[leq_element_MACRO,fii] = le_quad_q1_epd(coord_n_MACRO,new_coord_n_MACRO,bif_angle_rad);
               %otherwise
                  %error('This case is not implemented yet!');
            %end
            %
            m_VarAuxElem(p_leq_elem,iElem) = leq_elem;
            %
         end        

      end
      
      %% Cambio de estado del elemento.
      m_vectVHElemNew = m_VarHistElemNew(i_vectVHElem,iElem);
      %switch conshyp
         %case {4,11,44,53,54}
            ro = m_vectVHElemNew(1);
            ro_inj = m_vectVHElemNew(2);
         %otherwise
            %error('Cuadrángulo SDA bilineal: Factores ro: Modelo Constitutivo no definido.')
      %end   
      SLI_n = m_vectVHElemNew(5);
      CPI_n = e_DatElemSet.NumsidesCutCPF(iElem);

      %Para tener en cuenta que la media del criterio de carga/descarga puede ser muy levemente mayor que cero
      %pero son muy poco elementos los están dañando.
      limCeroCarga = 1e-5;
      % FROM (STATE==0) TO (STATE==1 || STATE==2)
      if ind_ActState_new==0

         % CASE 1
         %if ((RLI_n>0) && (ro<ro_inj)) || (((RLI_n>0) && (ro>ro_inj)) && (CPI_n==0 || CPI_n==4))
         %if ((SLI_n>0) && (ro<ro_inj)) || (((SLI_n>0) && (ro>ro_inj)) && (CPI_n==0 || CPI_n==4))
         if (SLI_n>limCeroCarga&&ro<ro_inj) || (SLI_n>limCeroCarga&&ro>ro_inj&&(CPI_n==0||CPI_n==4))
            ind_state_new = 1;
            ind_ActState_new = 1;
         % CASE 2
         %elseif (((SLI_n>0) && (ro>ro_inj)) || ((RLI_n>0) && (ro>ro_inj))) && (CPI_n==2)
         %elseif (((SLI_n>0) && (ro>ro_inj)) || ((SLI_n>0) && (ro>ro_inj))) && (CPI_n==2)
         elseif SLI_n>limCeroCarga&&ro>ro_inj&&CPI_n==2
            ind_state_new = 2;
            ind_ActState_new = 2;
         end

      % FROM (STATE==1) TO (STATE==0 || STATE==2)
      elseif ind_ActState_new==1

%          fprintf('Criterio de activación del SD ro=%f vs ro_inj=%f en el elemento %d.\n',...
%             ro,ro_inj,m_NumElem(iElem))
%          fprintf('Criterio de carga %f en el elemento %d.\n',...
%             SLI_n,m_NumElem(iElem))
         % CASE 3
         if (SLI_n<=limCeroCarga)
            ind_state_new = 1;
            ind_ActState_new = 0;
         % CASE 4
         elseif (SLI_n>limCeroCarga) && (ro>ro_inj) && (CPI_n==2)
            ind_state_new = 2;
            ind_ActState_new = 2;
         end

      % FROM (STATE==2) TO (STATE==0 || STATE==1)
      elseif ind_ActState_new==2

         % CASE 5
         if SLI_n<=limCeroCarga&&(CPI_n==0||CPI_n==4)
         %if (SLI_n<=limCeroCarga)
           %ind_state_new = 2;
           ind_ActState_new = 0;
         % CASE 6
         elseif (SLI_n>limCeroCarga) && (CPI_n==0 || CPI_n==4)
           %ind_state_new = 2;
           ind_ActState_new = 1;
         end
         %Si no se quiere que descargue del SD comentar la líneas previas.

      end
      %
      %Se realiza las siguientes operaciones si hubo cambio de estado.
      if ind_ActState_new~=ind_ActState_old
         %Impresión en archivo
         fprintf(fId,'Fin de Paso=%d (t=%f), Elemento=%d: Cambió de estado %d al %d.\n',...
            iStep,e_VG.Dtime*iStep,m_NumElem(iElem),ind_ActState_old,ind_ActState_new);
         %
         m_VarHistElemNew(p_indActSTmacro,iElem) = ind_ActState_new;
         m_VarHistElemNew(i_indST,iElem) = ind_state_new;
      end
      %
      %%
      ind_state_old = m_VarHistElemOld(i_indST,iElem);
      %Las operaciones siguientes se realizan únicamente en la primera vez que se activa el SD.
      if ind_ActState_new==2&&ind_state_old<=1

         %Impresión en pantalla
         fprintf('** Se activa la discontinuidad fuerte en el elemento %d.\n',m_NumElem(iElem));
         %Impresión en archivo
         fprintf(fId,'Fin de Paso=%d (t=%f), Elemento=%d: Se activa discontinuidad fuerte.\n',...
            iStep,e_VG.Dtime*iStep,m_NumElem(iElem));

         %Área del dominio localizado
         if isfield(e_DatMatSet,'m_ElemLocImp')
            m_ElemLoc = e_DatMatSet.m_ElemLocImp;
         else
            %Función de determinación del dominio localizado.
            m_ElemLoc = e_VarAuxPG(iPG,iElem).m_ElemLocCalc;
         end
         
         %Definición de las normal usada para orientar las normales micro.
         if ~isfield(e_DatMatSet,'sentNormales')||isempty(e_DatMatSet.sentNormales)
            %Por defecto se utiliza la normal de bifurcación seleccionada.
            sentNormales = 1;
         else
            sentNormales = e_DatMatSet.sentNormales;
         end
         if ~(length(sentNormales)>1)
            switch sentNormales
               case 1
                  %Se utiliza la normal de bifurcación seleccionado.
                  m_NormSent = m_VarAuxElem(p_n_tens,iElem);
                  m_NormSent = [m_NormSent(1);m_NormSent(4)];
               case 2
                  %Se utiliza la normal al path-crack.
                  m_NormSent = m_VarAuxElem(p_ref_vector,iElem);
               otherwise
                  error(['Operaciones de fin de paso: Cuadrángulo SDA bilineal del modelo multiescala ',...
                     'cohesivo: Tipo de vector para dar sentido a las normales micro no definido'])
            end
         end         

         e_DatSetMicro = e_DatMatSet.e_DatSet;
         nSetMicro = e_VGMicro.nSet;
         %Inicialización
         longFis         = 0;
         omegaMicroL     = 0;
         lMicro          = 0;
         m_ElemLoc_Bandas = false(1,length(m_ElemLoc));         
         N_barra = [0;0];
         for iSet = 1:nSetMicro
            e_DatElemMicro = e_DatSetMicro(iSet).e_DatElem;
            eltype = e_DatElemMicro.eltype;
            switch eltype
               case {2,4,8}
               case 32
                  m_IndElemSet = e_DatSetMicro(iSet). m_IndElemSet;
                  normal_micro = e_DatSetMicro(iSet).e_DatElem.normal_micro;
                  le_Elem = e_DatSetMicro(iSet).e_DatElem.le_Elem;
                  ksb = e_DatSetMicro(iSet).e_DatElem.ksb;
                  thickness = e_DatSetMicro(iSet).e_DatElem.thickness;
                  %Solo se cambia en el PG=6 que es el que lleva los datos del SD.
                  m_SentNormMicro = hvar_new(iPG,iElem).e_VarAux(iSet).VarAuxElem(1,:);
                  %
                  VolElem = e_DatMatSet.m_VolElem;
                  nElem_micro  = e_DatSetMicro(iSet).nElem;
                  for iElem_micro = 1:nElem_micro
                     ind_ELE_glob = m_IndElemSet(iElem_micro);
                     %Determinación de la normal media de las normales micro y de la longitud de fisura.
                     %Se calcula considerando solo los elementos que pertenecen al dominio localizado.
                     if m_ElemLoc(ind_ELE_glob)
                        m_NormMicroi = m_SentNormMicro(iElem_micro)*normal_micro(:,iElem_micro);
                        %Orientación de las normales micro.
                        if m_NormSent'*m_NormMicroi<0
                           %normal_micro(:,iElem_micro)= -normal_micro(:,iElem_micro);
                           m_NormMicroi = -m_NormMicroi;
                           %Se guarda en la variable auxiliar del elemento que lleva el signo.
                           m_SentNormMicro(iElem_micro) = -m_SentNormMicro(iElem_micro);
                        end
                        %Para tener en cuenta espesores no unitarios en la microcelda, se calcula la
                        %longitudes como superficie perpendicular al plano x-y.
                        deltaSFis = le_Elem(iElem_micro)*thickness;
                        longFis = longFis+deltaSFis;
                        m_ElemLoc_Bandas(ind_ELE_glob) = 1;
                        omegaMicroL = omegaMicroL+VolElem(ind_ELE_glob);
                        N_barra = N_barra+m_NormMicroi*deltaSFis;
                        %
                        lMicro = lMicro+ksb(iElem_micro);
                     end
                  end
                  %e_DatSetMicro(iSet).e_DatElem.normal_micro = normal_micro;
                  hvar_new(iPG,iElem).e_VarAux(iSet).VarAuxElem(1,:) = m_SentNormMicro;
               otherwise
                  error(['Operaciones de fin de paso: Cuadrángulo SDA bilineal del modelo multiescala ',...
                     'cohesivo: Normal micro no definido para este tipo de elemento.']);
            end
         end
         %e_DatMatSet.e_DatSet = e_DatSetMicro;

         %%%%%%%%%%%%%%% SACAR

%          incr = longFis / 6.5 ;
%          %incr = 0;
%          longFis = longFis+incr;
%          omegaMicroL = omegaMicroL+incr*0.043;
%          %    omegaMicroL=omegaMicroL +  incr*0.021 ;
%          N_barra = N_barra+ (N_barra/norm(N_barra))*incr;

         %%%%%%%%%%%%%%% SACAR
         %Determinación de la normal media.
         N_barra      = N_barra/longFis;
         facNormMicro = norm(N_barra);
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %facNormMicro = 1;
         %facNormMicro = norm([-sqrt(2)/2;sqrt(2)/2]*1+[sqrt(2)/2;sqrt(2)/2]*1)/2;
         %longFis = sqrt(2)*25;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         N_barra      = N_barra / norm(N_barra);
         m_VarAuxElem(p_n_tens    ,iElem)  = [N_barra(1) 0 0 N_barra(2) 0 N_barra(2) 0 N_barra(1)  ] ;
         
         %Porque excluye los elementos que no son elementos 32, esto en forma general está mal, ya que al
         %en general al dominio localizado puede también pertenecer elementos elásticos.
         m_ElemLoc    = m_ElemLoc_Bandas;

         lMicro       = lMicro/sum(m_ElemLoc);
         lMacro       = lMicro*e_DatMatSet.omegaMicro/omegaMicroL;      

         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %  redireccion de n_micro con N_barra promediado
%          for iSet = 1: nSet
%             e_DatElemMicro = e_DatSetMicro(iSet).e_DatElem;
%             eltype = e_DatElemMicro.eltype;
%             switch eltype
%                case 32
%                   %m_IndElemSet = e_DatSetMicro(iSet). m_IndElemSet;
%                   normal_micro = e_DatSetMicro(iSet).e_DatElem.normal_micro;
%                   nElem_micro  = e_DatSetMicro(iSet).nElem;
%                   for iElem_micro = 1:nElem_micro
%                      if (N_barra'* normal_micro(:,iElem_micro)) < 0
%                         fprintf(['Operaciones de fin de paso: Cuadrángulo SDA bilineal del modelo ',...
%                            'multiescala cohesivo: Elemento %d: Fue necesario un cambio en el sentido de ',...
%                            'normales micro por por segunda vez.\n'],m_NumElem(iElem))
%                         normal_micro(:,iElem_micro)= -normal_micro(:,iElem_micro);
%                      end
%                   end
%                   e_DatSetMicro(iSet).e_DatElem.normal_micro = normal_micro;
%             end
%          end
%          e_DatMatSet.e_DatSet= e_DatSetMicro;
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%         
         
         %Impresión en archivo de la normal utilizada
         fprintf(fId,'Fin de Paso=%d (t=%f), Elemento=%d: Ángulo de la normal seleccionada=%g°.\n',...
            iStep,e_VG.Dtime*iStep,m_NumElem(iElem),atan2(N_barra(2),N_barra(1))/pi*180);
         
         %Impresión en archivo del factor de tortuosidad
         fprintf(fId,'Fin de Paso=%d (t=%f), Elemento=%d: Factor de tortuosidad=%g.\n',...
            iStep,e_VG.Dtime*iStep,m_NumElem(iElem),facNormMicro);
         
         %Determinación de la tracción inicial (previo al primer paso que se resuelve con SD), usando el 
         %dominio localizado definitivo.
         m_Tracc = f_HomogTracc(e_DatSetMicro,hvar_new(iPG,iElem).e_VarEst,hvar_new(iPG,iElem).e_VarAux,...
            m_ElemLoc,facNormMicro,longFis,e_VGMicro);
         fprintf(fId,'Fin de Paso=%d (t=%f), Elemento=%d: Tracciones iniciales: Tx=%.15g y Ty=%.15g.\n',...
            iStep,e_VG.Dtime*iStep,m_NumElem(iElem),m_Tracc(1),m_Tracc(2));
         
         %Condiciones de borde original.
         m_CondBord = e_VGMicro.m_CondBord;
         %Se recupera las condiciones de borde impuesta en el estado 1, para mantenerla en el estado 2.
         if isfield(e_DatMatSet,'f_CambioCBExt')
            fprintf('** Se cambia la CB externa de la CU de todos los PG del elemento %d.\n',...
               m_NumElem(iElem))
            %n_NormalPer = [n_tens(1);n_tens(6)];
            %n_NormalPer = m_NormSent;
            n_NormalPer = m_VarAuxElem(p_ref_vector,iElem);
            %ang = 10/180*pi;
            %n_NormalPer = [cos(ang);sin(ang)];
            m_CondBord = e_DatMatSet.f_CambioCBExt(e_DatMatSet.xx,e_DatSetMicro,n_NormalPer,...
               e_VGMicro);
         end

         m_ConecFrontL =  f_LadFront(find(m_ElemLoc),e_DatMatSet.m_ElemVec,e_DatSetMicro,e_VGMicro);
         e_VGMicroAux = e_VGMicro; 
         switch e_DatMatSet.tipoCBL
            case 0  %No se impone CB en el dominio localizado
               %Se copia las mismas matrices que el punto regular (cuidado acá si se cambia la CB entre la
               %bifurcación y la discontinuidad fuerte, que ambos PG recupera la condiciones de borde
               %originales, es decir las previas a ese cambio).
               m_LinCondMicro = hvar_new(5,iElem).m_LinCond;
               doffMicro = hvar_new(5,iElem).doff;
               doflMicro = hvar_new(5,iElem).dofl;
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
               m_CondBord = f_CBTaylor(m_ElemLoc,e_DatSetMicro,m_CondBord,hvar_new(iPG,iElem).doff,e_VGMicro);
               %Alfredo: quito las condiciones de borde en las bandas
               %m_CondBord = m_CondBord(m_CondBord(:,2)~=55,:);
               %e_DatMatSet.e_VG.m_CondBord =  m_CondBord;
               e_VGMicroAux.m_CondBord = m_CondBord;
               [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicroAux,e_DatMatSet.xx,...
                  e_DatSetMicro,e_VGMicroAux.m_ConecFront);
            case 10 %CB Taylor en todo el dominio de la microcelda
               m_CondBord = [(1:e_VGMicro.nnod)',repmat([11,0,0],e_VGMicro.nnod,1)];
               e_VGMicroAux.m_CondBord = m_CondBord;
               [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicroAux,e_DatMatSet.xx,...
                  e_DatSetMicro,e_VGMicroAux.m_ConecFront);
            otherwise
               error(['Operaciones de fin de paso: Cuadrángulo SDA bilineal del modelo multiescala cohesivo: ',...
                  'Condición de borde del dominio localizado no implementado'])
         end
         %
         hvar_new(iPG,iElem).m_ElemLoc    = m_ElemLoc;
         hvar_new(iPG,iElem).omegaMicroL  = omegaMicroL;
         hvar_new(iPG,iElem).lMicro       = lMicro;
         hvar_new(iPG,iElem).lMacro       = lMacro;
         hvar_new(iPG,iElem).m_LinCond    = m_LinCondMicro;
         hvar_new(iPG,iElem).doff         = doffMicro;
         hvar_new(iPG,iElem).dofl         = doflMicro;
         hvar_new(iPG,iElem).longFis      = longFis;
         hvar_new(iPG,iElem).facNormMicro = facNormMicro;
         %        
      end
      
   end
   %
end
%Se cierra archivo de impresión de datos de bifurcación.
fclose(fId);

end

