
function [KT,Fint,c_GdlCond,e_VarEst_new,e_VarAux,c_TensorTang] = ...
  f_Par_MatGlobales_MICRO(xx,u,du,EPS_ROMI,eta,Q_ROMI,...
   c_GdlCond,e_DatSet,e_VarEst_new,e_VarEst_old,e_VarAux,DefMacro,e_VG)

   % FOR HROM
   % ********
   switch e_VG.MOD_TYPE
       
       case 2
           HROM_ONLINE ;
           return
           
       case 1
           ntens = e_VG.ntens;
           nModesEPS_TOT   = e_VG.nModesEPS_TOT;
           ROMI.KT_ROM     = zeros(nModesEPS_TOT);
           ROMI.KT_ROM_IMP = zeros(nModesEPS_TOT);
           ROMI.Fint_ROM   = zeros(nModesEPS_TOT,1);
           ROMI.Q_ROM      = zeros(nModesEPS_TOT,ntens) ;
           ROMI.Q_ROM_IMP  = zeros(nModesEPS_TOT,ntens) ;
           
           PHI_GEN = e_VG.PHI_GEN;

       otherwise
           ROMI=[];
           eta =[];
           ModoPHI_EPS = [];
   end

   % Recupera variables globales
   ntens   = e_VG.ntens;
   nSet = e_VG.nSet;
   ndime = e_VG.ndime;
   
   % Inicializaciones
   c_Ke = cell(nSet,1);
   c_Fint = cell(nSet,1);
   c_Fil = cell(nSet,1);
   c_Col = cell(nSet,1);
   c_FilFza = cell(nSet,1);
   c_TensorTang = cell(nSet,1); 
   
   %KT = sparse([]);
   %Fint=sparse([]);
   
   %Loop sobre el set de elementos (tipo y material del elemento)
   for iSet = 1:nSet
       e_DatMatSet = e_DatSet(iSet).e_DatMat;
       e_DatElemSet = e_DatSet(iSet).e_DatElem;
       nElem   = e_DatSet(iSet).nElem;
       npg     = e_DatElemSet.npg;
%        if e_DatMatSet.conshyp==1 %para los modelos elasticos se evita la evaluacion de Fint y K_e
%             ce=e_DatMatSet.ce;
%             CE= reshape(ce,[],1);
%             m_DofElem = e_DatSet(iSet).m_DofElem;
%             dofElemSet = m_DofElem(:);
%             uElemSet  = reshape(u(dofElemSet),[],nElem);
%             dofpe     = e_DatElemSet.dofpe;
%             m_BT      = e_DatSet(iSet).m_BT;
%             sigma_new = zeros(ntens,npg,nElem);
%             eps_new   = zeros(ntens,npg,nElem);
%             
%             emacro_xx= DefMacro{1}(1,1,1); 
%             emacro_yy= DefMacro{1}(2,1,1) ;
%             emacro_xy= DefMacro{1}(4,1,1) ;
%             
%             for jipg=1:npg
%                 for kdof=1:2:dofpe
%                     uelemex= reshape(uElemSet(kdof, :),[],1);
%                     uelemey= reshape(uElemSet(kdof+1, :),[],1);
%                     belemex = reshape(m_BT(1,kdof,jipg,:),[],1);
%                     belemey = reshape(m_BT(2,kdof+1,jipg,:),[],1);
%                     ex= reshape(eps_new(1,jipg,:),[],1);
%                     ey= reshape(eps_new(2,jipg,:),[],1);
%                     exy= reshape(eps_new(4,jipg,:),[],1);
%                     eps_new(1,jipg,:)  = ex + belemex.*uelemex;
%                     eps_new(2,jipg,:)  = ey + belemey.*uelemey;
%                     eps_new(4,jipg,:)  = exy + belemex.*uelemey + belemey.*uelemex;
%                     
%                 end
%                 
%                 sigma_new(1,jipg,:)= sigma_new(1,jipg,:)+ ce(1,1)*(eps_new(1,jipg,:)+emacro_xx)+ce(1,2)*(eps_new(2,jipg,:)+emacro_yy);
%                 sigma_new(2,jipg,:)= sigma_new(2,jipg,:)+ ce(2,1)*(eps_new(1,jipg,:)+emacro_xx)+ce(2,2)*(eps_new(2,jipg,:)+emacro_yy);
%                 sigma_new(3,jipg,:)= sigma_new(3,jipg,:)+ ce(3,1)*(eps_new(1,jipg,:)+emacro_xx)+ce(3,2)*(eps_new(2,jipg,:)+emacro_yy);
%                 sigma_new(4,jipg,:)= sigma_new(4,jipg,:)+ ce(4,4)*(eps_new(4,jipg,:)+emacro_xy);
%                     
%             end
%             
%             sigma_new = reshape(sigma_new, 4*npg,[]);
%             eps_new = reshape(eps_new, 4*npg,[]);
%             e_VarEst_new(iSet).sigma = sigma_new;
%             e_VarEst_new(iSet).eps   = eps_new;
% 
%            if isfield(e_DatMatSet,'esImplex')&&e_DatMatSet.esImplex
%                CE = repmat(CE,2*npg*nElem,1);
%                m_TensorTang = reshape(CE,ntens,ntens,2*npg,nElem);
%            else
%                CE = repmat(CE,npg*nElem,1);
%                m_TensorTang = reshape(CE,ntens,ntens,npg,nElem);
%            end
%        
%            c_TensorTang{iSet} = m_TensorTang;
%        else
           % Recuperación de variables
           %(evita llamar desde las estructura de los sets en el bucle de los elementos, ver que solo
           %recupera punteros)
           eltype  = e_DatElemSet.eltype;
           dofpe   = e_DatElemSet.dofpe;
           %conec = e_DatSet(iSet).conec;
           m_DofElem = e_DatSet(iSet).m_DofElem;
           sigma_old = e_VarEst_old(iSet).sigma;
           eps_old = e_VarEst_old(iSet).eps;
           hvar_old = e_VarEst_old(iSet).hvar;
           m_VarHistElemOld = e_VarEst_old(iSet).VarHistElem;
           aux_var = e_VarAux(iSet).VarAuxGP;
           m_VarAuxElem = e_VarAux(iSet).VarAuxElem;
           m_BT = e_DatSet(iSet).m_BT;
           m_DetJT = e_DatSet(iSet).m_DetJT;
           %siavare = e_DatSET(iSet).siavare;
           %npe = e_DatElemSet.npe;
           %m_IndElemSet = e_DatSet(iSet).m_IndElemSet;
           m_NumElem = e_DatSet(iSet).m_NumElem;
           m_DefMacroSet = DefMacro{iSet};
           
           % Inicializaciones
           %aux_var = zeros(siavare,nElem);
           
           %La verificación si esImplex es field de la estructura no sería necesario si a todos los modelos
           %constitutivos se agrega este campo.
           if isfield(e_DatMatSet,'esImplex')&&e_DatMatSet.esImplex
               %No se puede disminuir la matriz m_TensorTang cuando se calcula
               m_TensorTang = zeros(ntens,ntens,2*npg,nElem);
               %esImplex = 1;
           else
               m_TensorTang = zeros(ntens,ntens,npg,nElem);
               %esImplex = 0;
           end
           m_Ke = zeros(dofpe,dofpe,nElem);
           m_Fint = zeros(dofpe,nElem);
           %Como por convención no debería usarse las news dentro de las funciones, solo cuando converja, no se las
           %envía como argumento, y por lo tanto se evita esa transferencia (tiene más influencia en el caso de
           %cluster) de esas variables new. Por ello se debe generar en cada función la porción que corresponda de
           %las variables new (que tiene su costo, pero en el caso de un cluster es más conveniente), y tener en
           %cuenta que si es multiescala es una estructura y si es material estándar es una matriz.
           sigma_new   = e_VarEst_new(iSet).sigma;
           eps_new     = e_VarEst_new(iSet).eps;
           hvar_new    = e_VarEst_new(iSet).hvar;
           eps_fluct   = e_VarEst_new(iSet).eps_fluct;
           m_VarHistElemNew = e_VarEst_new(iSet).VarHistElem;
           
           % Tensor constitutivo elástico
           %En este momento se está precalculando en el read_data.
           %ce = c_elas(e_DatMatSet.young,e_DatMatSet.poiss,e_VG);
           
           % Coordenadas de los nodos de todos los elementos del set
           %Como se precalcula la matriz de deformación, las coordenadas de los nodos no se necesita.
           %coord_n = reshape(f_CoordElem(xx,conec')',[],npe,nElem);
           coord_n = [];
           
           % Grados de libertad y coordenadas de los nodos de los elementos del set
           %dofElemSet = f_DofElem(reshape(conec',1,[]),ndn);
           dofElemSet = m_DofElem(:);
           %m_IndDof = repmat(dofElemSet,dofpe,1);
           %Ver si hay alguna forma de hacerlo más rápido o mejor.
           %permute parece que hace copia de memoria, así que es probable que sea mejor usar dos repmat.
           %m_Fil = reshape(permute(reshape(m_IndDof,dofpe,dofpe,[]),[2,1,3]),[],1);
           %m_Col = m_IndDof(:);
           m_FilFza = dofElemSet';
           m_Fil = reshape(repmat(reshape(dofElemSet,dofpe,[]),dofpe,1),1,[]);
           m_Col = reshape(repmat(m_FilFza,dofpe,1),1,[]);
           uElemSet  = reshape(u(dofElemSet),[],nElem);
           duElemSet = reshape(du(dofElemSet),[],nElem);
           %m_DefMacroSet = DefMacro(:,m_IndElemSet);
           
           %Solo para debug, se guarda el número de set en e_VG (ya que esta variable se pasa a todas las
           %funciones).
           e_VG.iSet = iSet;
           
           %Se podría poner el parfor de dentro de cada eltype, así evitar hacer la verificación de tipo
           %de elemento para cada set.
           %parfor iElem = 1:nElem
           %for iElem = 1:nElem
           
           %Esta línea no puede estar si está el parfor activado.
           %e_VG.iElem = iElem;
           %No se puede modificar una variable definida previa al loop, si no es sliced y si no es
           %interpretada como de reducción, ya que no sabe como interpretarla el MatLab (puede haber
           %superposición de resultados, al hacer reducción). Por eso cada Lab debe tener su "copia",
           %para modificarla en forma independiente.
           %Lo que puede ser lento es copiar para cada lab esa e_VG_Aux, que puede ser grande.
           %En realidad como cada procesador tiene su copia local del e_VG (ya que MatLab realiza una
           %copia por cada Lab de todas la variables, excepto las sliced, donde solo copia la parte
           %que le corresponde al procesador), y al "copiarse" esta al e_VG_Aux, lo que único que se
           %hace es copiarse el puntero, ya que no se está modificando el e_VG_Aux. Luego lo que se
           %realiza es la modificación de un valor de un campo (field), donde supuestamente
           %MatLab al cambiar un field de una estructura no hace copia de toda la estructura de nuevo,
           %si solo del field, por lo tanto no tendría que ser mucho más lenta.
           %Otra sería pasar la variable iElem como argumento en las funciones que llama dentro del
           %parfor.
           %e_VG_Aux = e_VG;
           %e_VG_Aux.iElem = m_IndElemSet(iElem);
           %e_VG_Aux.iElemSet = iElem;
           %fprintf('Número de elemento: %d\n',iElem)
           
           % Fuerza interna y tensor tangente del elemento
           switch eltype
               case 2
%                   parfor iElem = 1:nElem
                   for iElem = 1:nElem
                       e_VG_Aux2 = e_VG;
                       e_VG_Aux2.iElemSet = iElem;
                       e_VG_Aux2.iElemNum = m_NumElem(iElem);
                       %Triángulo Estándar de 3 nodos
                       [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                           eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                           m_TensorTang(:,:,:,iElem)] = ...
                           f_MatElem_tria_t1(...
                           coord_n,uElemSet(:,iElem),hvar_old(:,iElem),aux_var(:,iElem),...
                           e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                           m_DefMacroSet(:,:,iElem),sigma_old(:,iElem),e_VG_Aux2);
                   end
               case 4     %Cuadrángulo Estándar de 4 nodos
                   %               parfor iElem = 1:nElem
                   m_VolElem    =  e_DatSet(iSet).m_VolElem;
                   for iElem = 1:nElem
                       e_VG_Aux4 = e_VG;
                       e_VG_Aux4.iElemSet = iElem;
                       e_VG_Aux4.iElemNum = m_NumElem(iElem);
                       %
                       [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                           eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                           m_TensorTang(:,:,:,iElem)] =...
                           f_MatElem_quad_q1(...
                           uElemSet(:,iElem),eps_old(:,iElem), hvar_old(:,iElem),aux_var(:,iElem),...
                           e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                           m_DefMacroSet(:,:,iElem),sigma_old(:,iElem),...
                           m_VolElem(iElem),e_VG_Aux4);
                   end
               case 31    % ELEMENTO TIPO BANDA : Cuadrángulo de 4 nodos
                   %               parfor iElem = 1:nElem
                   p_IntDissip  =  e_DatElemSet.pointersVHE.i_IntDissip  ;
                   p_IntEnergy  =  e_DatElemSet.pointersVHE.i_pIntEnergy;
                   m_IntDissip  =  m_VarHistElemOld (p_IntDissip ,:)  ;
                   m_pIntEnergy =  m_VarHistElemOld (p_IntEnergy ,:)  ;
                   m_pIntEnergy =  reshape(m_pIntEnergy, npg,[])  ;
                   ksb          =  e_DatElemSet.ksb;
                   for iElem = 1:nElem
                       e_VG_Aux4 = e_VG;
                       e_VG_Aux4.iElemSet = iElem;
                       e_VG_Aux4.iElemNum = m_NumElem(iElem);
                       %
                       %m_incrDissip =  m_VarAuxElem (p_incrDissip   ,iElem)  ;
                       %
                       [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                           eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                           m_TensorTang(:,:,:,iElem),m_IntDissip(iElem),...
                           m_pIntEnergy(:,iElem)] =...
                           f_MatElem_Banda_31(...
                           uElemSet(:,iElem),eps_old(:,iElem), hvar_old(:,iElem),aux_var(:,iElem),...
                           e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                           m_DefMacroSet(:,:,iElem),sigma_old(:,iElem),...
                           m_pIntEnergy(:,iElem), ksb(iElem),e_VG_Aux4);
                       %
                   end
                   m_VarHistElemNew (p_IntDissip ,: )      =  m_IntDissip(:)   ;
                   %%      m_pIntEnergy =  reshape(m_pIntEnergy,[],1)  ;
                   m_VarHistElemNew (p_IntEnergy ,: ) =  m_pIntEnergy   ;
               case 32
                   p_IntDissip  =  e_DatElemSet.pointersVHE.i_IntDissip  ;
                   p_IntEnergy  =  e_DatElemSet.pointersVHE.i_pIntEnergy;
                   m_IntDissip  =  m_VarHistElemOld (p_IntDissip ,:)  ;
                   m_pIntEnergy =  m_VarHistElemOld (p_IntEnergy ,:)  ;
                   m_pIntEnergy =  reshape(m_pIntEnergy, npg,[])  ;
                   ksb          =  e_DatElemSet.ksb;
                   parfor iElem = 1:nElem
%                    for iElem = 1:nElem
                       e_VG_Aux4 = e_VG;
                       e_VG_Aux4.iElemSet = iElem;
                       e_VG_Aux4.iElemNum = m_NumElem(iElem);
                       %
                       [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                           eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                           m_TensorTang(:,:,:,iElem),m_IntDissip(iElem),...
                           m_pIntEnergy(:,iElem)] =...
                           f_MatElem_Banda_Tria1(...
                           uElemSet(:,iElem),eps_old(:,iElem), hvar_old(:,iElem),aux_var(:,iElem),...
                           e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                           m_DefMacroSet(:,:,iElem),sigma_old(:,iElem),...
                           m_pIntEnergy(:,iElem), ksb(iElem),e_VG_Aux4);
                       %
                   end
                   m_VarHistElemNew (p_IntDissip ,: )      =  m_IntDissip(:)   ;
                   %%      m_pIntEnergy =  reshape(m_pIntEnergy,[],1)  ;
                   m_VarHistElemNew (p_IntEnergy ,: ) =  m_pIntEnergy   ;
               case 5
                   parfor iElem = 1:nElem
                       %                for iElem = 1:nElem
                       %Faltaría juntar las funciones de las fuerzas internas y matriz de rigidez
                       %elemental. También las funciones constitutivas.
                       %Falta precalcular las matrices de deformaciones y determinante de J. Esto está hecho,
                       %falta pasarselo a las funciones.
                       [m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),hvar_new(:,iElem),...
                           aux_var(:,iElem)] = elem_int_force_barra2D(coord_n,uElemSet(:,iElem),...
                           hvar_old(:,iElem),e_DatMatEl,DefMacro_new,e_VG);
                       [m_Ke(:,:,iElem),m_TensorTang(:,:,:,iElem)] = elem_stiff_barra2D(coord_n,...
                           uElemSet(:,iElem),hvar_old(:,iElem),aux_var(:,iElem),e_DatMatEl,ce,...
                           DefMacro_new,e_VG);
                   end
               case 7
                   error('Matrices Globales: Elemento finito no completo.')
               case 8      %Cuadrángulo BBar de 4 nodos
%                    parfor iElem = 1:nElem
%                    %for iElem = 1:nElem
%                        e_VG_Aux8 = e_VG;
%                        e_VG_Aux8.iElemSet = iElem;
%                        e_VG_Aux8.iElemNum = m_NumElem(iElem);
%                        [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
%                            eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
%                            m_TensorTang(:,:,:,iElem)] = ...
%                            f_MatElem_bbar_q1(...
%                            uElemSet(:,iElem),hvar_old(:,iElem),aux_var(:,iElem),...
%                            e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
%                            m_DefMacroSet(:,:,iElem),sigma_old(:,iElem),e_VG_Aux8);
%                    end 

                   
%               parfor iElem = 1:nElem
                %m_VolElem   =  e_DatSet(iSet).m_VolElem;
                p_leq_elem   =  e_DatElemSet.pointersVAE.p_leq_elem ;
                leq          =  m_VarAuxElem(p_leq_elem,:); 
                p_IntEnergy  =  e_DatElemSet.pointersVHE.i_pIntEnergy;
                m_pIntEnergy =  m_VarHistElemOld (p_IntEnergy ,:)  ; 
                m_pIntEnergy =  reshape(m_pIntEnergy, npg,[])  ; 
                
%               parfor iElem = 1:nElem
                for iElem = 1:nElem
                  e_VG_Aux8 = e_VG;
                  e_VG_Aux8.iElemSet = iElem;
                  e_VG_Aux8.iElemNum = m_NumElem(iElem);
                  %
                  e_PGsID = (m_NumElem(iElem)-1)*ntens*npg+1:m_NumElem(iElem)*ntens*npg ;
                  if e_VG.MOD_TYPE==1
                      ModoPHI_EPS = PHI_GEN(e_PGsID,:);
                      kin_Var = ModoPHI_EPS*EPS_ROMI;
                  else
                      kin_Var = uElemSet(:,iElem);
                  end
                  
                  [m_Ke(:,:,iElem),m_Fint(:,iElem),ROMI,sigma_new(:,iElem),eps_new(:,iElem),...
                       eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),...
                       m_TensorTang(:,:,:,iElem),m_pIntEnergy(:,iElem)] = ...
                  f_MatElem_bbar_q1_ROM(...
                       kin_Var,eta,ModoPHI_EPS,ROMI,eps_old(:,iElem),hvar_old(:,iElem),aux_var(:,iElem),...
                       e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                       m_DefMacroSet(:,:,iElem),sigma_old(:,iElem),leq(iElem),m_pIntEnergy(:,iElem),e_VG_Aux8);
               end
               m_VarHistElemNew (p_IntEnergy ,: ) =  m_pIntEnergy   ;                    
                   
                   
                            
              %LARGE DEFORMATIONS
              case 108  %Cuadrángulo de 4 nodos FBar.
                 %m_VolElem = e_DatSet(iSet).m_VolElem;
        %          parfor iElem = 1:nElem %JLM
                 for iElem = 1:nElem
                    e_VG_Aux108 = e_VG;
                    e_VG_Aux108.iElemSet = iElem;
                    e_VG_Aux108.iElemNum = m_NumElem(iElem);

                    if e_VG.MOD_TYPE==1 %&& ~e_VG.esME % Microscale ROM test
                        e_PGsID = (m_NumElem(iElem)-1)*ntens*npg+1:m_NumElem(iElem)*ntens*npg ;
                        ModoPHI_EPS = PHI_GEN(e_PGsID,:);
                        kin_Var = ModoPHI_EPS*EPS_ROMI;

                        [m_Ke(:,:,iElem),m_Fint(:,iElem),ROMI,sigma_new(:,iElem),eps_new(:,iElem),...
                           eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),m_TensorTang(:,:,:,iElem)] =...
                           f_MatElem_FBar_q1_ROM(...
                           uElemSet(:,iElem),kin_Var,eta,ModoPHI_EPS,ROMI,hvar_old(:,iElem),aux_var(:,iElem),...
                           m_VarAuxElem(:,iElem),e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                           m_DefMacroSet(:,iElem),eps_old(:,iElem),sigma_old(:,iElem),...
                           e_VG_Aux108);

                    else
                        [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                           eps_fluct(:,iElem),hvar_new(:,iElem),aux_var(:,iElem),m_TensorTang(:,:,:,iElem)] =...
                           f_MatElem_FBar_q1(...
                           uElemSet(:,iElem),hvar_old(:,iElem),aux_var(:,iElem),m_VarAuxElem(:,iElem),...
                           e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                           m_DefMacroSet(:,iElem),eps_old(:,iElem),sigma_old(:,iElem),...
                           e_VG_Aux108);                
                    end
        %             m_VarHistElemNew (p_IntEnergy ,: ) =  m_pIntEnergy ;
                 end
                   
                   
               case 10
                   m_resbT = zeros(ndime,nElem);
                   m_kbuT = zeros(ndime,dofpe,nElem);
                   m_invkbbT = zeros(ndime,ndime,nElem);
                   %Recuperación de variable de salto (variable interna condesada)
                   m_Beta = c_GdlCond{iSet,1};
                   %                parfor iElem = 1:nElem
                   for iElem = 1:nElem
                       e_VG_Aux10 = e_VG;
                       e_VG_Aux10.iElemSet = iElem;
                       e_VG_Aux10.iElemNum = m_NumElem(iElem);
                       %Triángulo Estándar de 3 nodos con discontinuidad fuerte (SDA)
                       [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),...
                           eps_fluct(:,iElem),hvar_new(:,iElem),m_VarHistElemNew(:,iElem),...
                           aux_var(:,iElem),m_VarAuxElem(:,iElem),m_TensorTang(:,:,:,iElem),...
                           m_resbT(:,iElem),m_kbuT(:,:,iElem),m_invkbbT(:,:,iElem)] = ...
                           f_MatElem_SDA_tria_t1(...
                           uElemSet(:,iElem),m_Beta(:,iElem),hvar_old(:,iElem),m_VarHistElemOld(:,iElem),...
                           aux_var(:,iElem),m_VarAuxElem(:,iElem),e_DatElemSet,e_DatMatSet,...
                           m_BT(:,:,:,iElem),m_DetJT(:,iElem),m_DefMacroSet(:,:,iElem),sigma_old(:,iElem),...
                           eps_old(:,iElem),e_VG_Aux10);
                   end
                   c_GdlCond{iSet,2} = m_resbT;
                   c_GdlCond{iSet,3} = m_kbuT;
                   c_GdlCond{iSet,4} = m_invkbbT;
               case 20    %Cuadrángulo de 4 nodos mixto con inyección de deformación
                   parfor iElem = 1:nElem
                       %                for iElem = 1:nElem
                       e_VG_Aux20 = e_VG;
                       e_VG_Aux20.iElemSet = iElem;
                       e_VG_Aux20.iElemNum = m_NumElem(iElem);
                       [m_Ke(:,:,iElem),m_Fint(:,iElem),sigma_new(:,iElem),eps_new(:,iElem),eps_fluct(:,iElem),...
                           hvar_new(:,iElem),m_VarHistElemNew(:,iElem),aux_var(:,iElem),m_TensorTang(:,:,:,iElem)] =...
                           f_MatElem_MixStrInj_quad_q1(...
                           uElemSet(:,iElem),hvar_old(:,iElem),m_VarHistElemOld(:,iElem),aux_var(:,iElem),...
                           m_VarAuxElem(:,iElem),e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),m_DetJT(:,iElem),...
                           m_DefMacroSet(:,:,iElem),sigma_old(:,iElem),e_VG_Aux20);
                   end
                   %##################################################################
               case 21    %Cuadrángulo de 4 nodos mixto con inyección de deformación
                   %##################################################################
                   %                m_resbT = zeros(ndime,nElem);
                   %                m_kbuT = zeros(ndime,dofpe,nElem);
                   %                m_invkbbT = zeros(ndime,ndime,nElem);
                   %
                   
                   m_Kbu_SDA      = zeros(ndime,dofpe,nElem);
                   m_KbbIN_SDA    = zeros(ndime,ndime,nElem);
                   m_Res_beta_SDA = zeros(ndime,nElem);
                   % MSL ELEMENT
                   m_Kbu_MSL      = zeros(ntens,dofpe,nElem);
                   m_KbbIN_MSL    = zeros(ntens,ntens,nElem);
                   m_Res_beta_MSL = zeros(ntens,nElem);
                   %Recuperación de variable de salto (variable interna condesada)
                   %beta_SDA     = c_GdlCond{iSet,1};
                   Dbeta_SDA    = c_GdlCond{iSet,2};
                   %beta_MSL     = c_GdlCond{iSet,3};
                   Dbeta_MSL    = c_GdlCond{iSet,4};
                   
                   i_indST        =  e_DatElemSet.pointersVHE.i_indST    ;
                   i_vectVHElem   =  e_DatElemSet.pointersVHE.i_vectVHElem ;
                   
                   p_condBif      =  e_DatElemSet.pointersVAE.p_condBif ;
                   p_elem_BifType =  e_DatElemSet.pointersVAE.p_elem_Biftype ;
                   p_leq_elem     =  e_DatElemSet.pointersVAE.p_leq_elem ;
                   p_kSD          =  e_DatElemSet.pointersVAE.p_kSD      ;
                   p_phi_grad     =  e_DatElemSet.pointersVAE.p_phi_grad  ;
                   p_n_tens       =  e_DatElemSet.pointersVAE.p_n_tens    ;
                   p_fii          =  e_DatElemSet.pointersVAE.p_fii       ;
                   p_injFactor    =  e_DatElemSet.pointersVAE.p_injFactor ;
                   
                   
                   m_indSTmacro_old  =  m_VarHistElemOld (i_indST        ,:)  ;
                   m_vectVHElem_old  =  m_VarHistElemOld (i_vectVHElem   ,:)  ;
                   m_indSTmacro      =  m_VarHistElemOld (i_indST        ,:)   ;
                   m_vectVHElem      =  m_VarHistElemOld (i_vectVHElem   ,:)  ;
                   
                   elem_BifType  =    m_VarAuxElem(p_elem_BifType,:)  ;
                   kSD           =    m_VarAuxElem(p_kSD,:);
                   
                   fact_inyect = e_DatElemSet.facIny;

                   %  parfor iElem = 1:nElem
                   for iElem = 1:nElem
                       
                       e_VG_Aux21 = e_VG;
                       e_VG_Aux21.iElemSet = iElem;
                       e_VG_Aux21.iElemNum = m_NumElem(iElem);
                       
                       
                       condBif         =  m_VarAuxElem(p_condBif,iElem) ;
                       leq_elem        =  m_VarAuxElem(p_leq_elem,iElem) ;
                       m_phi_grad      =  m_VarAuxElem(p_phi_grad,iElem) ;
                       m_phi_grad      =  reshape(m_phi_grad,4,2);
                       m_n_tens        =  m_VarAuxElem(p_n_tens,iElem) ;
                       m_n_tens        =  reshape(m_n_tens,4,2);
                       
                       [m_Ke(:,:,iElem),m_Kbu_SDA(:,:,iElem),m_Kbu_MSL(:,:,iElem),m_KbbIN_SDA(:,:,iElem),...
                           m_KbbIN_MSL(:,:,iElem),m_Fint(:,iElem),m_Res_beta_SDA(:,iElem),...
                           m_Res_beta_MSL(:,iElem),sigma_new(:,iElem), ...
                           hvar_new(:,iElem),eps_new(:,iElem),m_TensorTang(:,:,:,iElem),...
                           m_indSTmacro(:,iElem) ,elem_BifType(iElem),kSD(iElem),m_vectVHElem(:,iElem)] = ...
                           f_MatElem_quad_q1_SDA ...
                           (duElemSet(:,iElem),eps_old(:,iElem),...
                           Dbeta_SDA(:,iElem),Dbeta_MSL(:,iElem),...
                           aux_var(:,iElem),condBif, leq_elem ,m_phi_grad,m_n_tens,sigma_old(:,iElem),...
                           hvar_old(:,iElem),e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),...
                           m_DetJT(:,iElem),m_indSTmacro_old(iElem),e_VG_Aux21,...
                           m_vectVHElem_old(:,iElem),fact_inyect);
                       
                   end
                   
                   m_VarHistElemNew (i_indST   ,:)      = m_indSTmacro  ;
                   m_VarHistElemNew (i_vectVHElem   ,:) = m_vectVHElem  ;
                   
                   m_VarAuxElem(p_elem_BifType,:)       =  elem_BifType  ;
                   m_VarAuxElem(p_kSD,:)                =  kSD  ;
                   
                   
                   % SDA ELEMENT
                   c_GdlCond{iSet,5}  =   m_Kbu_SDA ;
                   c_GdlCond{iSet,6}  =   m_KbbIN_SDA ;
                   c_GdlCond{iSet,7}  =   m_Res_beta_SDA;
                   % MSL ELEMENT
                   c_GdlCond{iSet,8}  =   m_Kbu_MSL;
                   c_GdlCond{iSet,9}  =   m_KbbIN_MSL ;
                   c_GdlCond{iSet,10} =   m_Res_beta_MSL;
                   %##################################################################
               case 22    %Cuadrangulo de 4 nodos mixto con inyeccion de deformacion
                   % MODELO SANTA FE
                   %##################################################################
                   m_Kbu_SDA      = zeros(ndime,dofpe,nElem);
                   m_KbbIN_SDA    = zeros(ndime,ndime,nElem);
                   m_Res_beta_SDA = zeros(ndime,nElem);
                   % MSL ELEMENT
                   m_Kbu_MSL      = zeros(ntens,dofpe,nElem);
                   m_KbbIN_MSL    = zeros(ntens,ntens,nElem);
                   m_Res_beta_MSL = zeros(ntens,nElem);
                   m_stressTilde_new = zeros(ntens*npg,nElem);
                   
                   %Recuperaciï¿½n de variable de salto (variable interna condesada)
                   Dbeta_SDA    = c_GdlCond{iSet,2};
                   Dbeta_MSL    = c_GdlCond{iSet,4};
                   
                   i_indST        =  e_DatElemSet.pointersVHE.i_indST    ;
                   p_indActSTmacro=  e_DatElemSet.pointersVHE.p_indActSTmacro    ;
                   i_vectVHElem   =  e_DatElemSet.pointersVHE.i_vectVHElem ;
                   i_stressTilde  =  e_DatElemSet.pointersVHE.i_stressTilde ;
                   
                   p_condBif      =  e_DatElemSet.pointersVAE.p_condBif ;
                   p_elem_BifType =  e_DatElemSet.pointersVAE.p_elem_Biftype ;
                   p_leq_elem     =  e_DatElemSet.pointersVAE.p_leq_elem ;
                   p_kSD          =  e_DatElemSet.pointersVAE.p_kSD      ;
                   p_phi_grad     =  e_DatElemSet.pointersVAE.p_phi_grad  ;
                   p_n_tens       =  e_DatElemSet.pointersVAE.p_n_tens    ;
                   
                   NumsidesCutCPF = e_DatElemSet.NumsidesCutCPF;
                   
                   m_indSTmacro_old  =  m_VarHistElemOld (i_indST        ,:)  ;
                   m_indActSTmacro_old =  m_VarHistElemOld (p_indActSTmacro        ,:)  ;
                   m_vectVHElem_old  =  m_VarHistElemOld (i_vectVHElem   ,:)  ;
                   m_indSTmacro      =  m_VarHistElemOld (i_indST        ,:)   ;
                   m_indActSTmacro     =  m_VarHistElemOld (p_indActSTmacro        ,:)  ;
                   m_vectVHElem      =  m_VarHistElemOld (i_vectVHElem   ,:)  ;
                   m_stressTilde_old =  m_VarHistElemOld (i_stressTilde   ,:)  ;
                   
                   elem_BifType  =    m_VarAuxElem(p_elem_BifType,:)  ;
                   kSD           =    m_VarAuxElem(p_kSD,:);
                   
                   fact_inyect   =    e_VG.fact_inyect;
                   
                   %parfor iElem = 1:nElem
                   for iElem = 1:nElem
                       
                       e_VG_Aux21 = e_VG;
                       e_VG_Aux21.iElemSet = iElem;
                       e_VG_Aux21.iElemNum = m_NumElem(iElem);
                       
                       
                       condBif         =  m_VarAuxElem(p_condBif,iElem) ;
                       leq_elem        =  m_VarAuxElem(p_leq_elem,iElem) ;
                       m_phi_grad      =  m_VarAuxElem(p_phi_grad,iElem) ;
                       m_phi_grad      =  reshape(m_phi_grad,4,2);
                       m_n_tens        =  m_VarAuxElem(p_n_tens,iElem) ;
                       m_n_tens        =  reshape(m_n_tens,4,2);
                       
                       [m_Ke(:,:,iElem),m_Kbu_SDA(:,:,iElem),m_Kbu_MSL(:,:,iElem),m_KbbIN_SDA(:,:,iElem),...
                           m_KbbIN_MSL(:,:,iElem),m_Fint(:,iElem),m_Res_beta_SDA(:,iElem),...
                           m_Res_beta_MSL(:,iElem),sigma_new(:,iElem), ...
                           hvar_new(:,iElem),eps_new(:,iElem),m_TensorTang(:,:,:,iElem),...
                           m_indSTmacro(:,iElem) ,m_indActSTmacro(:,iElem),elem_BifType(iElem),...
                           kSD(iElem),m_vectVHElem(:,iElem),m_stressTilde_new(:,iElem)] = ...
                           f_MatElem_quad_q1_SDA_STAFE ...
                           (duElemSet(:,iElem),eps_old(:,iElem),...
                           Dbeta_SDA(:,iElem),Dbeta_MSL(:,iElem),...
                           aux_var(:,iElem),condBif, leq_elem ,m_phi_grad,m_n_tens,sigma_old(:,iElem),...
                           hvar_old(:,iElem),e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),...
                           m_DetJT(:,iElem),m_indSTmacro_old(iElem),m_indActSTmacro_old(iElem),...
                           e_VG_Aux21,m_vectVHElem_old(:,iElem),fact_inyect,...
                           NumsidesCutCPF(iElem),m_stressTilde_old(:,iElem));
                   end
                   
                   m_VarHistElemNew (i_indST   ,:)      = m_indSTmacro  ;
                   m_VarHistElemNew (p_indActSTmacro   ,:) = m_indActSTmacro  ;
                   m_VarHistElemNew (i_vectVHElem   ,:) = m_vectVHElem  ;
                   m_VarHistElemNew (i_stressTilde  ,:) = m_stressTilde_new  ;
                   
                   m_VarAuxElem(p_elem_BifType,:)       =  elem_BifType  ;
                   m_VarAuxElem(p_kSD,:)                =  kSD  ;
                   
                   
                   % SDA ELEMENT
                   c_GdlCond{iSet,5}  =   m_Kbu_SDA ;
                   c_GdlCond{iSet,6}  =   m_KbbIN_SDA ;
                   c_GdlCond{iSet,7}  =   m_Res_beta_SDA;
                   % MSL ELEMENT
                   c_GdlCond{iSet,8}  =   m_Kbu_MSL;
                   c_GdlCond{iSet,9}  =   m_KbbIN_MSL ;
                   c_GdlCond{iSet,10} =   m_Res_beta_MSL;
                   %##################################################################
               case 23    %Cuadrangulo de 4 nodos mixto SD
                   % MODELO SANTA FE
                   %##################################################################
                   m_Kbu_SDA      = zeros(ndime,dofpe,nElem);
                   m_KbbIN_SDA    = zeros(ndime,ndime,nElem);
                   m_Res_beta_SDA = zeros(ndime,nElem);
                   % MSL ELEMENT
                   m_Kbu_MSL      = zeros(ntens,dofpe,nElem);
                   m_KbbIN_MSL    = zeros(ntens,ntens,nElem);
                   m_Res_beta_MSL = zeros(ntens,nElem);
                   m_stressTilde_new = zeros(ntens*npg,nElem);
                   
                   %Recuperaciï¿½n de variable de salto (variable interna condesada)
                   Dbeta_SDA    = c_GdlCond{iSet,2};
                   Dbeta_MSL    = c_GdlCond{iSet,4};
                   
                   i_indST        =  e_DatElemSet.pointersVHE.i_indST    ;
                   p_indActSTmacro=  e_DatElemSet.pointersVHE.p_indActSTmacro    ;
                   i_vectVHElem   =  e_DatElemSet.pointersVHE.i_vectVHElem ;
                   i_stressTilde  =  e_DatElemSet.pointersVHE.i_stressTilde ;
                   
                   p_condBif      =  e_DatElemSet.pointersVAE.p_condBif ;
                   p_elem_BifType =  e_DatElemSet.pointersVAE.p_elem_Biftype ;
                   p_leq_elem     =  e_DatElemSet.pointersVAE.p_leq_elem ;
                   p_kSD          =  e_DatElemSet.pointersVAE.p_kSD      ;
                   p_phi_grad     =  e_DatElemSet.pointersVAE.p_phi_grad  ;
                   p_n_tens       =  e_DatElemSet.pointersVAE.p_n_tens    ;
                   
                   NumsidesCutCPF = e_DatElemSet.NumsidesCutCPF;
                   
                   m_indSTmacro_old  =  m_VarHistElemOld (i_indST        ,:)  ;
                   m_indActSTmacro_old =  m_VarHistElemOld (p_indActSTmacro        ,:)  ;
                   m_vectVHElem_old  =  m_VarHistElemOld (i_vectVHElem   ,:)  ;
                   m_indSTmacro      =  m_VarHistElemOld (i_indST        ,:)   ;
                   m_indActSTmacro     =  m_VarHistElemOld (p_indActSTmacro        ,:)  ;
                   m_vectVHElem      =  m_VarHistElemOld (i_vectVHElem   ,:)  ;
                   m_stressTilde_old =  m_VarHistElemOld (i_stressTilde   ,:)  ;
                   
                   elem_BifType  =    m_VarAuxElem(p_elem_BifType,:)  ;
                   kSD           =    m_VarAuxElem(p_kSD,:);
                   
                  
                   %parfor iElem = 1:nElem
                   for iElem = 1:nElem
                       
                       e_VG_Aux21 = e_VG;
                       e_VG_Aux21.iElemSet = iElem;
                       e_VG_Aux21.iElemNum = m_NumElem(iElem);
                       
                       
                       condBif         =  m_VarAuxElem(p_condBif,iElem) ;
                       leq_elem        =  m_VarAuxElem(p_leq_elem,iElem) ;
                       m_phi_grad      =  m_VarAuxElem(p_phi_grad,iElem) ;
                       m_phi_grad      =  reshape(m_phi_grad,4,2);
                       m_n_tens        =  m_VarAuxElem(p_n_tens,iElem) ;
                       m_n_tens        =  reshape(m_n_tens,4,2);
                       
                       [m_Ke(:,:,iElem),m_Kbu_SDA(:,:,iElem),m_Kbu_MSL(:,:,iElem),m_KbbIN_SDA(:,:,iElem),...
                           m_KbbIN_MSL(:,:,iElem),m_Fint(:,iElem),m_Res_beta_SDA(:,iElem),...
                           m_Res_beta_MSL(:,iElem),sigma_new(:,iElem), ...
                           hvar_new(:,iElem),eps_new(:,iElem),m_TensorTang(:,:,:,iElem),...
                           m_indSTmacro(:,iElem) ,m_indActSTmacro(:,iElem),elem_BifType(iElem),...
                           kSD(iElem),m_vectVHElem(:,iElem),m_stressTilde_new(:,iElem)] = ...
                           f_MatElem_quad_q23_SDA_STAFE ...
                           (duElemSet(:,iElem),eps_old(:,iElem),...
                           Dbeta_SDA(:,iElem),Dbeta_MSL(:,iElem),...
                           aux_var(:,iElem),condBif, leq_elem ,m_phi_grad,m_n_tens,sigma_old(:,iElem),...
                           hvar_old(:,iElem),e_DatElemSet,e_DatMatSet,m_BT(:,:,:,iElem),...
                           m_DetJT(:,iElem),m_indSTmacro_old(iElem),m_indActSTmacro_old(iElem),...
                           e_VG_Aux21,m_vectVHElem_old(:,iElem),fact_inyect,...
                           NumsidesCutCPF(iElem),m_stressTilde_old(:,iElem));
                   end
                   
                   m_VarHistElemNew (i_indST   ,:)      = m_indSTmacro  ;
                   m_VarHistElemNew (p_indActSTmacro   ,:) = m_indActSTmacro  ;
                   m_VarHistElemNew (i_vectVHElem   ,:) = m_vectVHElem  ;
                   m_VarHistElemNew (i_stressTilde  ,:) = m_stressTilde_new  ;
                   
                   m_VarAuxElem(p_elem_BifType,:)       =  elem_BifType  ;
                   m_VarAuxElem(p_kSD,:)                =  kSD  ;
                   
                   
                   % SDA ELEMENT
                   c_GdlCond{iSet,5}  =   m_Kbu_SDA ;
                   c_GdlCond{iSet,6}  =   m_KbbIN_SDA ;
                   c_GdlCond{iSet,7}  =   m_Res_beta_SDA;
                   % MSL ELEMENT
                   c_GdlCond{iSet,8}  =   m_Kbu_MSL;
                   c_GdlCond{iSet,9}  =   m_KbbIN_MSL ;
                   c_GdlCond{iSet,10} =   m_Res_beta_MSL;
           end
           
           %end
           
           e_VarEst_new(iSet).sigma = sigma_new;
           e_VarEst_new(iSet).eps = eps_new;
           e_VarEst_new(iSet).hvar = hvar_new;
           e_VarEst_new(iSet).eps_fluct = eps_fluct;
           e_VarEst_new(iSet).VarHistElem = m_VarHistElemNew;
           e_VarAux(iSet).VarAuxGP = aux_var;
           e_VarAux(iSet).VarAuxElem = m_VarAuxElem;
           c_Ke{iSet}         = m_Ke(:);
           c_Fint{iSet}       = m_Fint(:);
           c_Fil{iSet}        = m_Fil;
           c_Col{iSet}        = m_Col;
           c_FilFza{iSet}     = m_FilFza;
           c_TensorTang{iSet} = m_TensorTang;
%        end
       
   end
   
   
    if e_VG.MOD_TYPE == 1

        % Ensamble de matriz de rigidez global        
        KT.Implex = [ ROMI.KT_ROM  e_VG.IntPhiTGen;
            e_VG.IntPhiTGen'   zeros(ntens-1) ];
        
        KT.Impli = [ ROMI.KT_ROM_IMP  e_VG.IntPhiTGen;
            e_VG.IntPhiTGen'   zeros(ntens-1) ];
        
        % Ensamble de matriz de fuerzas internas global           
        Fint = [ROMI.Fint_ROM;e_VG.IntPhiTGen'*EPS_ROMI];
        
        KT.Q_ROM = [ROMI.Q_ROM; zeros(ntens-1,ntens)] ;
        KT.Q_ROM_IMP = [ROMI.Q_ROM_IMP; zeros(ntens-1,ntens)];
        
    else
        
        % Ensamble de matriz de fuerzas internas global
        Fint = sparse([c_FilFza{:}],1,cat(1,c_Fint{:}));
        
        % Ensamble de matriz de rigidez global
        KT = sparse([c_Fil{:}],[c_Col{:}],cat(1,c_Ke{:}));
        
    end
   
%    % Ensamble de matriz de fuerzas internas global
%    Fint = sparse([c_FilFza{:}],1,cat(1,c_Fint{:}));
%    
%    % Ensamble de matriz de rigidez global
%    KT = sparse([c_Fil{:}],[c_Col{:}],cat(1,c_Ke{:}));
%    
% %    if ~isempty(e_VG.K_GlobElast)
% %        if length(e_VG.K_GlobElast)> length(KT)
% %           KT(length(e_VG.K_GlobElast),length(e_VG.K_GlobElast))=0;
% %           Fint(length(e_VG.K_GlobElast),1)=0;
% %        end
% %        KT=KT+e_VG.K_GlobElast;
% %        Fint= Fint + e_VG.K_GlobElast*u+ e_VG.BTCe_Glob * DefMacro{1}(:,1,1) ;
% %    end
   
end