function [m_VarAuxElem,m_VarHistElemNew] = ...
   f_OperPosConvSDA_QuadQ1(uElem,xx,conec,m_VarHistElemNew,m_VarHistElemOld,EF_sidesCutCPF,...
   m_VarAuxElem,m_CT,e_DatSet_iSet,e_VG,hvar_new,hvar_old)

iStep = e_VG.istep;
%
e_DatMat = e_DatSet_iSet.e_DatMat;
e_DatElem = e_DatSet_iSet.e_DatElem;
nElem = e_DatSet_iSet.nElem;
dN_xy = e_DatSet_iSet.dN_xy;
m_NumElem  = e_DatSet_iSet.m_NumElem;
%
eltype = e_DatElem.eltype;
esImplex = e_DatMat.esImplex;
conshyp = e_DatMat.conshyp;
%
p_condBif     = e_DatElem.pointersVAE.p_condBif ;
p_leq_elem    = e_DatElem.pointersVAE.p_leq_elem ;
%p_normal_bif  = e_DatElem.pointersVAE.p_normal_bif ;
pAngBif = e_DatElem.pointersVAE.pAngBif;
p_n_tens      = e_DatElem.pointersVAE.p_n_tens ;
p_injFactor   = e_DatElem.pointersVAE.p_injFactor ;
p_ref_vector  = e_DatElem.pointersVAE.p_ref_vector ;
p_indActST = e_DatElem.pointersVHE.p_indActSTmacro    ;
i_indST       = e_DatElem.pointersVHE.i_indST    ;
i_vectVHElem = e_DatElem.pointersVHE.i_vectVHElem;

%Nombre de archivo que se guarda los datos de bifurcación.
nomArchBif = [e_VG.fileCompleto,'.elemBif'];
%Se abre el archivo.
fId = fopen(nomArchBif,'at');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%parfor iElem = 1:nElem
for iElem = 1:nElem   
   
   condBif        =  m_VarAuxElem(p_condBif   ,iElem) ;
   %leq_elem       =  m_VarAuxElem(p_leq_elem  ,iElem) ;
   %n_bifurc       =  m_VarAuxElem(p_normal_bif,iElem) ;
   %n_bifurc       =  reshape(n_bifurc,2,2);
   %inject_factor  =  m_VarAuxElem(p_injFactor  ,iElem) ;
   ind_ActState_new = m_VarHistElemNew (p_indActST,iElem);
   ind_ActState_old = m_VarHistElemNew (p_indActST,iElem);
   ind_state_new = m_VarHistElemNew(i_indST,iElem);
   
   e_VG.iElemSet = iElem;
   e_VG.iElemNum = m_NumElem(iElem);
   
   %Para realizar operaciones una sola vez después que bifurcó.
   %Notar que esto solo funciona en el caso que se permita hacer operaciones inmediatamente que bifurcó.
   condBifPrev = condBif;
   %Se guarda la normal de bifurcación seleccionada previa (old) para ver si vuelve a cambiar durante el
   %periodo intermedio entre la bifurcación y la activación de la SD, y así realizar ciertas operaciones
   %únicamente en esos casos. También para recuperar los datos de bifurcación previos cuando se detecta 
   %bifurcación.
   n_tensPrev = m_VarAuxElem(p_n_tens,iElem);
   n_tensPrev = [n_tensPrev(1);n_tensPrev(6)];
   %m_AngBifPrev = m_VarAuxElem(pAngBif,iElem);
   
   % Operaciones que se realiza hasta que se detecta la bifurcación.
   if ~condBif
   
      % Evaluación de la condición de bifurcación.
      if esImplex
         m_CTBif = m_CT(:,:,12,iElem);
      else
         m_CTBif = m_CT(:,:,6,iElem);
      end
%      [condBif,~,leq_elem,n_bifurc,n_tens] = ...
%         bifurcation_condition_MACRO...
%         (m_CTBif,leq_elem,condBif,e_VG,xx,...
%         conec(iElem,:),uElem(:,iElem),e_DatElem,dN_xy(:,:,iElem),e_DatSet_iSet);
      [condBif,m_AngBif] = f_CondBifct(m_CTBif,e_VG);
      
      %Ver si es necesario seleccionar una normal previo a la bifurcación. Pareciera que se usa en la 
      %SDA_Properties para determinar el GradPhi y los nodos libres, pero esto se usa únicamente cuando se
      %inyecta.
      m_NormBif = [cos(m_AngBif);sin(m_AngBif)];
      n_tens = Normal_vector_selection(xx,conec(iElem,:),uElem(:,iElem),m_NormBif,0,...
            e_VG,e_DatElem,e_DatMat,dN_xy(:,:,iElem),EF_sidesCutCPF(iElem),...
            hvar_new(:,iElem),hvar_old(:,iElem));
      
      %Se guarda en notación de Voigt el vector normal.
      m_VarAuxElem(p_condBif   ,iElem) = condBif     ;
      %m_VarAuxElem(p_leq_elem  ,iElem) = leq_elem    ;
      %m_VarAuxElem(p_normal_bif,iElem) = reshape(n_bifurc,4,1);
      %Se almacena los ángulos de bifurcación.
      m_VarAuxElem(pAngBif,iElem) = m_AngBif;
      %m_VarAuxElem(p_n_tens,iElem) = reshape(n_tens,8,1);
      m_VarAuxElem(p_n_tens,iElem) = n_tens(:);
      
   end
   %elseif (condBif == 1) && (indSTmacro==1)
   %else
   if condBif
      
      %Se ejecuta entre la bifurcación y la activación de la SD, y mientras no se active la SD (una vez
      %activada la bifurcación por primera vez no se entra más).
      %Cómo esta parte está previo a realizar un cambio de estado en ind_ActState_new (más abajo), se obliga
      %entrar al menos una vez (por ejemplo que se pase directamente de ind_ActState de 0 a 2,
      %ind_ActState_new==2 && ind_ActState_old==0, que no verifica ind_state_new<=1, pero de esta forma sí lo
      %hace).
      if ind_state_new<=1
         
         %Como criterio se utiliza los datos del análisis de bifurcación previo a detectar bifurcación
         %(detQ<0), estos han mostrado ser más cercanos a una normal de bifurcación más coherente con el
         %problema.
         %En modelos constitutivos simples (como los fenomenológicos) no conviene tomar el ángulo obtenido del
         %análisis previo a la bifurcación, ya que generalmente es el elástico y entonces los ángulos no
         %resultan adecuados. En donde resultó más conveniente fue en los modelos multiescala, que se tiene
         %mucha disipación estable.
         %if ~condBifPrev
            %m_VarAuxElem(pAngBif,iElem) = m_AngBifPrev;
            %Esta línea siguiente no es necesario ya que en la función Normal_vector_selection se selecciona
            %una nueva normal con los dos ángulos posibles.
            %m_VarAuxElem(p_n_tens,iElem) = n_tensPrev(:);
         %end
         %
         m_AngBif = m_VarAuxElem(pAngBif,iElem)';
         m_NormBif = [cos(m_AngBif);sin(m_AngBif)];
         
         %Selección de normal de bifurcación
         %[n_tens,m_tensor(:,:,iElem),ref_vector,inject_factor] ...
         [n_tens,m_tensor(:,:,iElem),ref_vector] = ...
            Normal_vector_selection(xx,conec(iElem,:),uElem(:,iElem),m_NormBif,0,...
            e_VG,e_DatElem,e_DatMat,dN_xy(:,:,iElem),EF_sidesCutCPF(iElem),...
            hvar_new(:,iElem),hvar_old(:,iElem));
         %
         %
         m_VarAuxElem(p_n_tens,iElem) = n_tens(:);
         m_VarAuxElem(p_ref_vector,iElem) = ref_vector;
         %m_VarAuxElem(p_injFactor,iElem) = inject_factor;
                  
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
            
            % Calculo longitud de Oliver (Oliver X./89)
            m_CoordElem = f_CoordElem(xx,conec(iElem,:));
            npe = e_DatElem.npe ;
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
            %   case {4,10,21,22,23} % Elementos cuadrilaterales - indistinto de si es enriquecido o no.
                  leq_elem = le_quad_q1_epd(m_CoordElem,new_coord_n_MACRO,bif_angle_rad);
                  %[leq_element_MACRO,fii] = le_quad_q1_epd(coord_n_MACRO,new_coord_n_MACRO,bif_angle_rad);
            %   otherwise
            %      error('This case is not implemented yet!');
            %end
            %
            m_VarAuxElem(p_leq_elem,iElem) = leq_elem;
            
         end
         
      end
      
      %Cambio de estado del elemento.
      m_vectVHElemNew = m_VarHistElemNew(i_vectVHElem,iElem);
      %switch conshyp
      %   case {4,11,44,53,54}
            ro = m_vectVHElemNew(1);
            ro_inj = m_vectVHElemNew(2);
      %   otherwise
      %      error('Cuadrángulo SDA bilineal: Factores ro: Modelo Constitutivo no definido.')
      %end
      %rTrial = m_vectVHElemNew(9);
      %facInyRo = 0.5;
      %
      SLI_n = m_vectVHElemNew(5);
      CPI_n = e_DatElem.NumsidesCutCPF(iElem);

      %Para tener en cuenta que la media del criterio de carga/descarga puede ser muy levemente mayor que cero
      %pero son muy poco elementos los están dañando.
      limCeroCarga = 1e-5;
      % FROM (STATE==0) TO (STATE==1 || STATE==2)
      if ind_ActState_new==0
         
         % CASE 1
         %if ((RLI_n>0) && (ro<ro_inj)) || (((RLI_n>0) && (ro>ro_inj)) && (CPI_n==0 || CPI_n==4))
         %if ((SLI_n>0) && (ro<ro_inj)) || (((SLI_n>0) && (ro>ro_inj)) && (CPI_n==0 || CPI_n==4))
         if (SLI_n>limCeroCarga&&ro<ro_inj) || (SLI_n>limCeroCarga&&ro>ro_inj&&(CPI_n==0||CPI_n==4))
         %Tendría que ser rActual para delimitar cuando cambia del estado 0 al 1.
         %if (rTrial>facInyRo*e_DatMat.r_0&&ro<ro_inj) || (rTrial>facInyRo*e_DatMat.r_0&&ro>ro_inj&&(CPI_n==0||CPI_n==4))
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
         
         % CASE 3
         if (SLI_n<=limCeroCarga)
         %if rTrial<facInyRo*e_DatMat.r_0
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
         %if rTrial<facInyRo*e_DatMat.r_0&&(CPI_n==0||CPI_n==4)
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
         m_VarHistElemNew(p_indActST,iElem) = ind_ActState_new;
         m_VarHistElemNew(i_indST,iElem) = ind_state_new;
      end
      %
      %
      ind_state_old = m_VarHistElemOld(i_indST,iElem);
      %Las operaciones siguientes se realizan únicamente en la primera vez que se activa el SD.
      if ind_ActState_new==2&&ind_state_old<=1
         
         %Impresión en pantalla
         fprintf('** Se activa la discontinuidad fuerte en el elemento %d.\n',m_NumElem(iElem));
         %Impresión en archivo
         fprintf(fId,'Fin de Paso=%d (t=%f), Elemento=%d: Se activa discontinuidad fuerte.\n',...
            iStep,e_VG.Dtime*iStep,m_NumElem(iElem));
        
      end
     
   end
   
end
%Se cierra archivo de impresión de datos de bifurcación.
fclose(fId);
%
end
