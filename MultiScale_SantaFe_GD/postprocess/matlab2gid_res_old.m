function matlab2gid_res(iStep,in,u,c_GdlCond,e_DatSet,e_VarEst,e_VarAux,DefMacro,uTotal,ELOC,e_VG)

%**************************************************
%* RUTINA PARA GENERAR LOS ARCHIVOS DE RESULTADOS *
%* PARA SER POSTPROCESADOS CON GID                *
%* Archivo: NAME.flavia.res                       *
%**************************************************

% Variables globales
ndime        = e_VG.ndime;
%nnod         = e_VG.nnod;
%nElem        = e_VG.nElem;
struhyp      = e_VG.struhyp;
%conshyp      = e_VG.conshyp;
%npg          = e_VG.npg;
ntens        = e_VG.ntens;
%sihvarpg     = e_VG.sihvarpg;
nSet = e_VG.nSet;
%filename_res = e_VG.filename_res;
fileCompleto = e_VG.fileCompleto;

filename_res = [fileCompleto,'.flavia.res'];
fid_res = fopen(filename_res,'at');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fluctuacion de desplazamientos
fprintf(fid_res,'Result "Desplazamientos//Fluctuaciones" "Load Analysis" %d Vector OnNodes\n',iStep);
if ndime==2
    fprintf(fid_res,'ComponentNames "X-DESPL" "Y-DESPL"\n');
elseif ndime==3
    fprintf(fid_res,'ComponentNames "X-DESPL" "Y-DESPL" "Z-DESPL"\n');
end
fprintf(fid_res,'Values\n');
format = ['%d',repmat(' %.15g',1,ndime),'\n'];
fprintf(fid_res,format,[in';[u(1:ndime:length(u)) u(2:ndime:length(u))]']);
fprintf(fid_res,'End Values\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desplazamientos totales
fprintf(fid_res,'Result "Desplazamientos//Totales" "Load Analysis" %d Vector OnNodes\n',iStep);
if ndime==2
    fprintf(fid_res,'ComponentNames "X-DESPL" "Y-DESPL"\n');
elseif ndime==3
    fprintf(fid_res,'ComponentNames "X-DESPL" "Y-DESPL" "Z-DESPL"\n');
end
fprintf(fid_res,'Values\n');
format = ['%d',repmat(' %.15g',1,ndime),'\n'];
fprintf(fid_res,format,[in';uTotal]);
fprintf(fid_res,'End Values\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crack path field
if e_VG.exist_CrackPath
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Crack Path Field
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(fid_res,'Result "Crack Path Field" "Load Analysis" %d Scalar OnNodes\n',iStep);
    fprintf(fid_res,'Values\n');
    format = ['%d', '  %.15g','\n'];
    fprintf(fid_res,format,[in'; e_VG.smooth_dalpha']);
    fprintf(fid_res,'End Values\n');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ALPHA Crack Path
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf(fid_res,'Result "ALPHA Crack Path" "Load Analysis" %d Scalar OnNodes\n',iStep);
    fprintf(fid_res,'Values\n');
    format = ['%d', '  %.15g','\n'];
    fprintf(fid_res,format,[in'; e_VG.smooth_alpha']);
    fprintf(fid_res,'End Values\n');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % vectores normales solo en los elementos bifurcados
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for iSet = 1:nSet
%         eltype = e_DatSet(iSet).e_DatElem.eltype;
%         switch eltype
%             case 53
%                 nElem = e_DatSet(iSet).nElem;
%                 e_DatElemSet = e_DatSet(iSet).e_DatElem;
%                 m_NumElem = e_DatSet(iSet).m_NumElem;
%                 m_VarAuxElem = e_VarAux(iSet).VarAuxElem;
%                 m_VarHistElem = e_VarEst(iSet).VarHistElem;
%                 nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
%                 %*****************************************************************************************************
%                 % VECTOR NORMAL M DE ANALISIS DE BIFURCACION del elemento finito (SEGUNDO VECTOR DEL ANALISIS DE BIFURCACION)
%                 fprintf(fid_res,['Result "Elemento Finito//SDA//Vector Normal" "Load Analysis" %d ',...
%                     'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
%                 if ndime==2
%                     fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP"\n');
%                 elseif ndime==3
%                     fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP" "Z-COMP"\n');
%                 end
%                 fprintf(fid_res,'Values\n');
%                 
%                 i_indST           =  e_DatElemSet.pointersVHE.i_indST    ;
%                 m_indSTmacro       =  m_VarHistElem (i_indST        ,:)  ;
%                 
%                 p_n_tens      = e_DatElemSet.pointersVAE.p_n_tens ;
%                 n_tens        =  m_VarAuxElem(p_n_tens    ,:)  ;
%                 n_tens         =  reshape(n_tens,4,2,[]);
%                 
%                 indi_elem     = find (m_indSTmacro == 2);
%                 nx           =n_tens(1,1,:);
%                 ny           =n_tens(2,2,:);
%                 
%                 for iElem = 1:nElem
%                     if find(iElem==indi_elem)
%                         fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),nx(iElem),ny(iElem),0]');
%                     else
%                         fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),0,0,0]');
%                     end
%                 end
%                 fprintf(fid_res,'End Values  \n');
%         end
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % alpha_sin_suavizar
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iSet = 1:nSet
        eltype = e_DatSet(iSet).e_DatElem.eltype;
        switch eltype
            case 21
                nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
                nomGaussSet     = ['GP_Set_',num2str(iSet)];
                m_NumElem       = e_DatSet(iSet).m_NumElem;
                p_varHistSmooth = 5;
                m_VarHistElem   = e_VarEst(iSet).VarHistElem(p_varHistSmooth,:);
                %
                fprintf(fid_res,['Result "Elemento Finito//SDA//alpha_sin_suavizar" "Load Analysis" %d ',...
                    'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
                fprintf(fid_res,'Values\n');
                %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
                %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
                format = ['%d',repmat(' %.15g',1,1),'\n'];
                %Se imprime la variable del elemento en un punto de gauss central.
                fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
                fprintf(fid_res,'End Values \n');
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % dissipacion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iSet = 1:nSet
        eltype = e_DatSet(iSet).e_DatElem.eltype;
        switch eltype
            case 21
                nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
                nomGaussSet     = ['GP_Set_',num2str(iSet)];
                m_NumElem       = e_DatSet(iSet).m_NumElem;
                p_varHistSmooth = 2;  % disipacion
                m_VarHistElem   = e_VarEst(iSet).VarHistElem(p_varHistSmooth,:);
                %
                fprintf(fid_res,['Result "Elemento Finito//SDA//disipacion" "Load Analysis" %d ',...
                    'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
                fprintf(fid_res,'Values\n');
                %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
                %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
                format = ['%d',repmat(' %.15g',1,1),'\n'];
                %Se imprime la variable del elemento en un punto de gauss central.
                fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
                fprintf(fid_res,'End Values \n');
        end
    end
    for iSet = 1:nSet
        conshyp  = e_DatSet(iSet).e_DatMat.conshyp;
        sihvarpg = e_DatSet(iSet).e_DatMat.sihvarpg;
        switch conshyp
            % case 11
            case 21
                nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
                nomGaussSet     = ['GP_Set_',num2str(iSet)];
                m_NumElem       = e_DatSet(iSet).m_NumElem;
                p_varHistSmooth = sihvarpg*5+1;
                hvar = e_VarEst(iSet).hvar(p_varHistSmooth,:);
                %
                fprintf(fid_res,['Result "Elemento Finito//SDA//danio_elemental" "Load Analysis" %d ',...
                    'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
                fprintf(fid_res,'Values\n');
                %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
                %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
                format = ['%d',repmat(' %.15g',1,1),'\n'];
                %Se imprime la variable del elemento en un punto de gauss central.
                fprintf(fid_res,format,[m_NumElem;hvar]);
                fprintf(fid_res,'End Values \n');
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % variable_inyeccion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iSet = 1:nSet
        eltype = e_DatSet(iSet).e_DatElem.eltype;
        switch eltype
            case 21
                nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
                m_NumElem        = e_DatSet(iSet).m_NumElem;
                p_varHistSmooth  = 2;
                m_VarHistElem    = e_VarEst(iSet).VarHistElem(p_varHistSmooth,:);
                %
                fprintf(fid_res,['Result "Elemento Finito//SDA//variable_inyeccion" "Load Analysis" %d ',...
                    'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
                fprintf(fid_res,'Values\n');
                %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
                %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
                format = ['%d',repmat(' %.15g',1,1),'\n'];
                %Se imprime la variable del elemento en un punto de gauss central.
                fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
                fprintf(fid_res,'End Values \n');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
for iSet = 1:nSet
    
    nElem = e_DatSet(iSet).nElem;
    e_DatElemSet = e_DatSet(iSet).e_DatElem;
    e_DatMatSet  = e_DatSet(iSet).e_DatMat;
    m_NumElem    = e_DatSet(iSet).m_NumElem;
    m_IndElemSet = e_DatSet(iSet).m_IndElemSet;
    npg = e_DatElemSet.npg;
    NPG_ELE = npg;
    conshyp = e_DatMatSet.conshyp;
    eltype = e_DatElemSet.eltype;
    sihvarpg = e_DatMatSet.sihvarpg;
    siavarpg = e_DatMatSet.siavarpg;
    stress = e_VarEst(iSet).sigma;
    strain = e_VarEst(iSet).eps;
    eps_fluct = e_VarEst(iSet).eps_fluct;
    % Variables internas seg�n el modelo constitutivo
    hvar_new = e_VarEst(iSet).hvar;
    
    if (eltype== 21) && (conshyp~=53)
        hvar_new =reshape(hvar_new, sihvarpg,npg,[]);
        hvar_new = hvar_new(:,1:4,:);
        hvar_new =reshape(hvar_new,sihvarpg*4,[]);
    end
    
    aux_var = e_VarAux(iSet).VarAuxGP;
    m_VarAuxElem = e_VarAux(iSet).VarAuxElem;
    m_VarHistElem = e_VarEst(iSet).VarHistElem;
    
    nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
    nomGaussSet = ['GP_Set_',num2str(iSet)];
    
    NPG_ELE = npg;
    switch eltype
       case {4,31}
           switch struhyp
               case {1,2}
                   %Deformaci�n plana y Tensi�n plana
                   tipoDatAlmac = 'PlainDeformationMatrix';
                   nomComponente = '"Tensi�n XX" "Tensi�n YY" "Tensi�n XY" "Tensi�n ZZ"';
               otherwise
                   error(['GidPost: Tensiones en los puntos de Gauss: Componentes: No definidas para ',...
                       'esta hip�tesis de carga']);
           end
           fprintf(fid_res,['Result "Tensiones//Sobre punto de Gauss" "Load Analysis" %d %s ',...
               'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
           fprintf(fid_res,'ComponentNames %s\n',nomComponente);
           fprintf(fid_res,'Values\n');
           format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,npg)];
            fprintf(fid_res,format,[m_NumElem;stress]);
            fprintf(fid_res,'End Values\n');
        case 10    % Tri�ngulo de 3 nodos con discontinuidades fuertes (SDA)
            %          %Se imprime solo el punto de gauss regular, los valores de puntos de gauss singular se
            %          %descarta (ver si no imprimir los valores de este PG de otra forma)
            %          npg = 1;
            %          stress = stress(1:ntens,:);
            %          strain = strain(1:ntens,:);
            %          eps_fluct = eps_fluct(1:ntens,:);
            %          hvar_new = hvar_new(1:sihvarpg,:);
            
            %Como previo a la bifurcaci�n y a la activaci�n de la strong discontinuity (SD) el punto de
            %gauss singular no se actualiza (para ahorrar memoria) se actualiza en la etapa previa solo
            %para imprensi�n, luego tendr�a que seguir distintos caminos.
            m_indNocondBif = m_VarAuxElem(1,:)<2;
            if any(m_indNocondBif)
                stress(ntens+1:2*ntens,m_indNocondBif) = stress(1:ntens,m_indNocondBif);
                strain(ntens+1:2*ntens,m_indNocondBif) = strain(1:ntens,m_indNocondBif);
                eps_fluct(ntens+1:2*ntens,m_indNocondBif) = eps_fluct(1:ntens,m_indNocondBif);
                hvar_new(sihvarpg+1:2*sihvarpg,m_indNocondBif) = hvar_new(1:sihvarpg,m_indNocondBif);
                aux_var(siavarpg+1,m_indNocondBif) = aux_var(1:siavarpg,m_indNocondBif);
            end
            
            %*************************
            % Salto en el elemento
            %*************************
            fprintf(fid_res,['Result "Elemento Finito//SDA//Vector Salto" "Load Analysis" %d ',...
                'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            %Nombre de las componentes
            %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
            fprintf(fid_res,'ComponentNames "BetaX" "BetaY"\n');
            %Valores
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
            format = ['%d',repmat(' %.15g',1,ndime),'\n'];
            m_Beta = c_GdlCond{iSet,1};
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;m_Beta]);
            fprintf(fid_res,'End Values \n');
            
            %********************************************
            % Vector Normal a la fisura en el elemento
            %********************************************
            fprintf(fid_res,['Result "Elemento Finito//SDA//Vector Normal a la fisura" ',...
                '"Load Analysis" %d Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
            fprintf(fid_res,'ComponentNames "nX" "nY"\n');
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
            %caso 3D.
            format = ['%d',repmat(' %.15g',1,ndime),'\n'];
            m_VecNormal = m_VarAuxElem([2;7],:);
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;m_VecNormal]);
            fprintf(fid_res,'End Values \n');
            
            %**********************************************
            % Vector tangente a la fisura en el elemento
            %**********************************************
            %Se asume la convenci�n la regla de la mano derecha para el producto t x n.
            %Se est� asumiendo un solo punto de gauss utilizado en la integraci�n.
            m_VecTan = zeros(2,nElem);
            m_VecTan(1,:) = m_VecNormal(2,:);
            m_VecTan(2,:) = -m_VecNormal(1,:);
            %
            fprintf(fid_res,['Result "Elemento Finito//SDA//Vector Tangente a la fisura" ',...
                '"Load Analysis" %d Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
            fprintf(fid_res,'ComponentNames "tX" "tY"\n');
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
            %caso 3D.
            format = ['%d',repmat(' %.15g',1,ndime),'\n'];
            %
            fprintf(fid_res,format,[m_NumElem;m_VecTan]);
            fprintf(fid_res,'End Values \n');
            
            % Vector Tracci�n
            fprintf(fid_res,['Result "Elemento Finito//SDA//Vector Tracci�n" "Load Analysis" %d ',...
                'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            %Nombre de las componentes
            %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
            fprintf(fid_res,'ComponentNames "TX" "TY"\n');
            %Valores
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
            format = ['%d',repmat(' %.15g',1,ndime),'\n'];
            m_VecTracc = m_VarAuxElem(19:20,:);
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;m_VecTracc]);
            fprintf(fid_res,'End Values \n');
        case 20   %Quadr�ngulo de 4 nodos mixto con inyecci�n de deformaci�n
            % Condici�n de bifurcaci�n o activaci�n de la inyecci�n (dominio mixto)
            m_CondBif = m_VarAuxElem(1,:);
            fprintf(fid_res,['Result "Elemento Finito//Q4MID//Dominio Mixto" "Load Analysis" %d ',...
                'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            %Valores
            fprintf(fid_res,'Values\n');
            format = '%d %d\n';
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;m_CondBif]);
            fprintf(fid_res,'End Values \n');
            % Tensi�n estabilizada en los PGs
            m_TensEstab = m_VarHistElem;
            switch struhyp
                case {1,2}
                    %Deformaci�n plana y Tensi�n plana
                    tipoDatAlmac = 'PlainDeformationMatrix';
                    nomComponente = '"Tensi�n XX" "Tensi�n YY" "Tensi�n XY" "Tensi�n ZZ"';
                otherwise
                    error(['GidPost: Tensiones en los puntos de Gauss: Componentes: No definidas para ',...
                        'esta hip�tesis de carga']);
            end
            fprintf(fid_res,['Result "Elemento Finito//Q4MID//Tensiones Estabilizadas" "Load Analysis" %d %s ',...
                'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
            fprintf(fid_res,'ComponentNames %s\n',nomComponente);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,npg)];
            if (struhyp==1||struhyp==2)
                % Deformaci�n Plana y Tensi�n Plana
                m_Ind = reshape((1:npg*ntens)',ntens,[]);
                m_Ind([ntens,ntens-1],:) = m_Ind([ntens-1,ntens],:);
                fprintf(fid_res,format,[m_NumElem;m_TensEstab(m_Ind(:),:)]);
            else
                fprintf(fid_res,format,[m_NumElem;m_TensEstab]);
            end
            fprintf(fid_res,'End Values\n');
        case 21    % Quadrilatero 4 nodos con discontinuidades fuertes (SDA)
            NPG_ELE = npg-2;
            p_n_tens       =  e_DatElemSet.pointersVAE.p_n_tens    ;
            p_elem_BifType =  e_DatElemSet.pointersVAE.p_elem_Biftype ;
            
            m_elem_BifType  =  m_VarAuxElem (p_elem_BifType , : ) ;
            m_n_tens        =  m_VarAuxElem (p_n_tens  , : ) ;
            %******************************
            % Indice de estado del elemento
            %******************************
            fprintf(fid_res,['Result "Elemento Finito//SDA//Indice Estado" "Load Analysis" %d ',...
                'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            format = ['%d',' %.15g','\n'];
            fprintf(fid_res,format,[m_NumElem; m_elem_BifType]);
            fprintf(fid_res,'End Values \n');
             %*************************
            % Salto en el elemento
            %*************************
            fprintf(fid_res,['Result "Elemento Finito//SDA//Vector Salto" "Load Analysis" %d ',...
                'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            %Nombre de las componentes
            %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
            fprintf(fid_res,'ComponentNames "BetaX" "BetaY"\n');
            %Valores
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
            format = ['%d',repmat(' %.15g',1,ndime),'\n'];
            m_Beta = c_GdlCond{iSet,1};
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;m_Beta]);
            fprintf(fid_res,'End Values \n');
            %********************************************
            % Vector Normal a la fisura en el elemento
            %********************************************
            fprintf(fid_res,['Result "Elemento Finito//SDA//Vector Normal a la fisura" ',...
                '"Load Analysis" %d Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
            fprintf(fid_res,'ComponentNames "nX" "nY"\n');
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
            %caso 3D.
            format = ['%d',repmat(' %.15g',1,ndime),'\n'];
            m_VecNormal = m_n_tens([1;4],:);
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;m_VecNormal]);
            fprintf(fid_res,'End Values \n');

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % vectores normales solo en los elementos bifurcados
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %for iSet = 1:nSet
            %eltype = e_DatSet(iSet).e_DatElem.eltype;
            %switch eltype
            %case 53
            %nElem = e_DatSet(iSet).nElem;
            %e_DatElemSet = e_DatSet(iSet).e_DatElem;
            %m_NumElem = e_DatSet(iSet).m_NumElem;
            %m_VarAuxElem = e_VarAux(iSet).VarAuxElem;
            %m_VarHistElem = e_VarEst(iSet).VarHistElem;
            %nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
            %*****************************************************************************************************
            % VECTOR NORMAL M DE ANALISIS DE BIFURCACION del elemento finito (SEGUNDO VECTOR DEL ANALISIS DE BIFURCACION)
            fprintf(fid_res,['Result "Elemento Finito//SDA//Vector Normal" "Load Analysis" %d ',...
                'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            if ndime==2
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP"\n');
            elseif ndime==3
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP" "Z-COMP"\n');
            end
            fprintf(fid_res,'Values\n');
            
            i_indST           =  e_DatElemSet.pointersVHE.i_indST    ;
            m_indSTmacro       =  m_VarHistElem (i_indST        ,:)  ;
            
            p_n_tens      = e_DatElemSet.pointersVAE.p_n_tens ;
            n_tens        =  m_VarAuxElem(p_n_tens    ,:)  ;
            n_tens         =  reshape(n_tens,4,2,[]);
            
            indi_elem     = find (m_indSTmacro == 2);
            nx           =n_tens(1,1,:);
            ny           =n_tens(2,2,:);
            
            for iElem = 1:nElem
                if find(iElem==indi_elem)
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),nx(iElem),ny(iElem),0]');
                else
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),0,0,0]');
                end
            end
            fprintf(fid_res,'End Values  \n');
            %end
            %end
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            %********************************************
            % Tensiones
            %********************************************
            switch struhyp
                case {1,2} %Deformaci�n plana y Tensi�n plana
                    tipoDatAlmac = 'PlainDeformationMatrix';
                    nomComponente = '"Tensi�n XX" "Tensi�n YY" "Tensi�n XY" "Tensi�n ZZ"';
                case 3 %Tridimensional
                    tipoDatAlmac = 'Matrix';
                    nomComponente = 'Tensi�n XX" "Tensi�n YY" "Tensi�n ZZ" "Tensi�n XY" "Tensi�n XZ" "Tensi�n YZ';
            end
            STRESS= reshape(stress,ntens,npg,[]);
            STRESS=STRESS(:,[1:4],:);
            STRESS= reshape(STRESS,ntens*NPG_ELE,[]);
            
            fprintf(fid_res,['Result "Tensiones//Sobre punto de Gauss" "Load Analysis" %d %s ',...
                'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
            fprintf(fid_res,'ComponentNames %s\n',nomComponente);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,NPG_ELE)];
            fprintf(fid_res,format,[m_NumElem;STRESS]);
            fprintf(fid_res,'End Values\n');
                        %Normal material: normal del an�lisis de bifurcaci�n del material
%             if conshyp==11
%                 if ndime==2
%                     s_NomComp = 'ComponentNames "nX" "nY"\n';
%                 elseif ndime==3
%                     s_NomComp = 'ComponentNames "nX" "nY" "nZ"\n';
%                 end
%                 if e_DatMatSet.esImplex
%                     posNormalCrit1 = 4;
%                     posNormalCrit2 = 5;
%                     posCondBif = 6;
%                 else
%                     posNormalCrit1 = 3;
%                     posNormalCrit2 = 4;
%                     posCondBif = 5;
%                 end
%                 %Normal cr�tica 1
%                 fprintf(fid_res,['Result "Modelo Constitutivo//Normal material 1" ',...
%                     '"Load Analysis" %d Vector OnGaussPoints "',nomGaussSet,'"\n'],iStep);
%                 %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
%                 fprintf(fid_res,s_NomComp);
%                 fprintf(fid_res,'Values\n');
%                 format = ['%d',repmat([repmat(' %.15g',1,ndime),'\n'],1,npg)];
%                 angCrit = hvar_new(posNormalCrit1:sihvarpg:end,:);
%                 m_VecNormal = zeros(ndime*npg,nElem);
%                 m_VecNormal(1:ndime:end) = cos(angCrit);
%                 m_VecNormal(2:ndime:end) = sin(angCrit);
%                 fprintf(fid_res,format,[m_NumElem;m_VecNormal]);
%                 fprintf(fid_res,'End Values \n');
%                 %Normal cr�tica 2
%                 fprintf(fid_res,['Result "Modelo Constitutivo//Normal material 2" ',...
%                     '"Load Analysis" %d Vector OnGaussPoints "',nomGaussSet,'"\n'],iStep);
%                 %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
%                 fprintf(fid_res,s_NomComp);
%                 fprintf(fid_res,'Values\n');
%                 format = ['%d',repmat([repmat(' %.15g',1,ndime),'\n'],1,npg)];
%                 angCrit = hvar_new(posNormalCrit2:sihvarpg:end,:);
%                 m_VecNormal = zeros(ndime*npg,nElem);
%                 m_VecNormal(1:ndime:end) = cos(angCrit);
%                 m_VecNormal(2:ndime:end) = sin(angCrit);
%                 fprintf(fid_res,format,[m_NumElem;m_VecNormal]);
%                 fprintf(fid_res,'End Values \n');
%                 %Condici�n de bifurcaci�n
%                 fprintf(fid_res,['Result "Modelo Constitutivo//Condici�n de bifurcaci�n" ',...
%                     '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
%                 fprintf(fid_res,'Values\n');
%                 format = ['%d',repmat(' %d\n',1,npg)];
%                 fprintf(fid_res,format,[m_NumElem;hvar_new(posCondBif:sihvarpg:end,:)]);
%                 fprintf(fid_res,'End Values\n');
%             end

        otherwise% Tri�ngulo de 3 nodos con discontinuidades fuertes (SDA)
            fprintf(fid_res,['Result "Tensiones//Sobre punto de Gauss" "Load Analysis" %d %s ',...
                'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
            fprintf(fid_res,'ComponentNames %s\n',nomComponente);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,NPG_ELE)];
            if (struhyp==1||struhyp==2)
                % Deformaci�n Plana y Tensi�n Plana
                m_Ind = reshape((1:npg*ntens)',ntens,[]);
                m_Ind([ntens,ntens-1],:) = m_Ind([ntens-1,ntens],:);
                fprintf(fid_res,format,[m_NumElem;stress(m_Ind(:),:)]);
            else
                fprintf(fid_res,format,[m_NumElem;stress]);
            end
            fprintf(fid_res,'End Values\n');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fluctuacion de Deformaciones sobre los puntos de Gauss
            switch struhyp
                case {1,2}
                    %Deformaci�n plana y Tensi�n plana
                    tipoDatAlmac = 'PlainDeformationMatrix';
                    nomComponente = '"Deformaci�n XX" "Deformaci�n YY" "Deformaci�n XY" "Deformaci�n ZZ"';
                case 3
                    %Tridimensional
                    tipoDatAlmac = 'Matrix';
                    nomComponente = ['Deformaci�n XX" "Deformaci�n YY" "Deformaci�n ZZ" "Deformaci�n XY" ',...
                        '"Deformaci�n XZ" "Deformaci�n YZ'];
                    %case 4
                    %Axisimetr�a
                    %case 5
                    %Barras 2D
                otherwise
                    error(['GidPost: Fluctuaciones de Deformaciones en los puntos de Gauss: Componentes: ',...
                        'No definidas para esta hip�tesis de carga']);
            end
            fprintf(fid_res,['Result "Deformaci�n//Fluctuaciones//Sobre punto de Gauss" ',...
                '"Load Analysis" %d %s OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
            fprintf(fid_res,'ComponentNames %s\n',nomComponente);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,npg)];
            if struhyp==1
                % Deformaci�n Plana
                fprintf(fid_res,format,[m_NumElem;eps_fluct(m_Ind(:),:)]);
            else
                fprintf(fid_res,format,[m_NumElem;eps_fluct]);
            end
            fprintf(fid_res,'End Values\n');
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Deformaciones totales sobre los puntos de Gauss
            fprintf(fid_res,['Result "Deformaci�n//Totales//Sobre punto de Gauss" ',...
                '"Load Analysis" %d %s OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
            fprintf(fid_res,'ComponentNames %s\n',nomComponente);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,NPG_ELE)];
            if struhyp==1
                % Deformaci�n Plana
                fprintf(fid_res,format,[m_NumElem;strain(m_Ind(:),:)]);
            else
                fprintf(fid_res,format,[m_NumElem;strain]);
            end
            fprintf(fid_res,'End Values\n');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Deformaci�n macro total aplicada en cada PG
            fprintf(fid_res,['Result "Deformaci�n//Macro aplicada//Sobre elemento" ',...
                '"Load Analysis" %d %s OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
            fprintf(fid_res,'ComponentNames %s\n',nomComponente);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,npg)];
            if struhyp==1||struhyp==2
                m_DefMacro = reshape(DefMacro{iSet}([1,2,4,3],:,:),ntens*npg,nElem);
                % Deformaci�n Plana y Tensi�n Plana
                fprintf(fid_res,format,[m_NumElem;m_DefMacro]);
            else
                m_DefMacro = reshape(DefMacro{iSet},ntens*npg,nElem);
                fprintf(fid_res,format,[m_NumElem;m_DefMacro]);
            end
            fprintf(fid_res,'End Values\n');
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Condici�n de incremento de deformaci�n total que se ha verificado
            fprintf(fid_res,['Result "Deformaci�n//Dominio Localizado//',...
                'Sobre elemento usado" "Load Analysis" %d Scalar OnGaussPoints "',...
                nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            format = '%d %d\n';
            fprintf(fid_res,format,[m_NumElem;ELOC(1,m_IndElemSet)]);
            fprintf(fid_res,'End Values\n');
            
            if size(ELOC,1)>1
                % Condici�n de incremento de deformaci�n total que se ha verificado
                fprintf(fid_res,['Result "Deformaci�n//Dominio Localizado//',...
                    'Sobre elemento seg�n incremento de deformaci�n" "Load Analysis" %d ',...
                    'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
                fprintf(fid_res,'Values\n');
                format = '%d %d\n';
                fprintf(fid_res,format,[m_NumElem;ELOC(2,m_IndElemSet)]);
                fprintf(fid_res,'End Values\n');
            end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    switch conshyp
        %       eltype = e_DatElemSet.eltype;
        % NPG_ELE = npg;
        
        
        case {1,52}
        case {2,9}
            % Plasticidad J2
            % Deformacion plastica equivalente sobre puntos de Gauss
            fprintf(fid_res,['Result  "Modelo Constitutivo//Def. Pl�stica Equiv.//',...
                'Sobre punto de Gauss" "Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
            fprintf(fid_res,format,[m_NumElem;hvar_new(ntens+1:sihvarpg:end,:)]);
            fprintf(fid_res,'End Values\n');
            
            % �ndice de carga fload
            fprintf(fid_res,['Result  "Modelo Constitutivo//�ndice de Carga//Sobre punto de Gauss" ',...
                '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat(' %d\n',1,NPG_ELE)];
            fload = hvar_new(ntens+2:sihvarpg:end,:);
            fprintf(fid_res,format,[m_NumElem;fload]);
            fprintf(fid_res,'End Values\n');
            
            % �ndice de carga fload (por elemento)
            fprintf(fid_res,['Result  "Modelo Constitutivo//�ndice de Carga//',...
                'Sobre elemento (any PG)" "Load Analysis" %d Scalar OnGaussPoints "',...
                nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            format = '%d %d\n';
            %Se adopta si un elemento tiene un PG con fload!=0, el elemento adopta un fload 1.
            fprintf(fid_res,format,[m_NumElem;any(fload)]);
            fprintf(fid_res,'End Values\n');
         case {11}  % elemento banda de la microcelda modelo Barcelona
            % Variable de da�o
            fprintf(fid_res,['Result "Modelo Constitutivo//Damage//Sobre punto de Gauss" ',...
                '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
            fprintf(fid_res,format,[m_NumElem;hvar_new(1:sihvarpg:end,:)]);
            %fprintf(fid_res,format,[m_NumElem;aux_var(1:siavarpg:end,:)]);
            fprintf(fid_res,'End Values\n');
        case {4,5,10}   % antes modelo 11
            % Da�o isotropo est�ndar y regularizado
            % Variable de da�o
            fprintf(fid_res,['Result "Modelo Constitutivo//Damage//Sobre punto de Gauss" ',...
                '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
            %fprintf(fid_res,format,[m_NumElem;hvar_new(3:sihvarpg:end,:)]);
            fprintf(fid_res,format,[m_NumElem;aux_var(1:siavarpg:end,:)]);
            fprintf(fid_res,'End Values\n');
            %*************************************************************************************
            % �ndice de carga FLoad
            %*************************************************************************************
            fprintf(fid_res,['Result "Modelo Constitutivo//�ndice de Carga//Sobre punto de Gauss" ',...
                '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
            %Para utilizar la tabla de resultados previamente creada.
            fprintf(fid_res,['ResultRangesTable ','"Factor de Carga de Da�o Isotr�pico"\n']);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat(' %d\n',1,NPG_ELE)];
            %fprintf(fid_res,format,[m_NumElem;hvar_new(4:sihvarpg:end,:)]);
            fprintf(fid_res,format,[m_NumElem;aux_var(2:siavarpg:end,:)]);
            fprintf(fid_res,'End Values\n');
            
            %Impresi�n de los valores centrales para el elemento Q1MixtoInjDef en forma separada
            if eltype==20
                % Variable de da�o
                fprintf(fid_res,['Result "Modelo Constitutivo//Damage//Punto de gauss Central" ',...
                    '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
                fprintf(fid_res,'Values\n');
                format = '%d %.15g\n';
                fprintf(fid_res,format,[m_NumElem;aux_var((npg-1)*siavarpg+1,:)]);
                fprintf(fid_res,'End Values\n');
                
                % Factor de carga (FLoad)
                fprintf(fid_res,['Result "Modelo Constitutivo//�ndice de Carga//Punto de gauss Central" ',...
                    '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
                %Para utilizar la tabla de resultados previamente creada.
                fprintf(fid_res,['ResultRangesTable ','"Factor de Carga de Da�o Isotr�pico"\n']);
                fprintf(fid_res,'Values\n');
                format = '%d %d\n';
                fprintf(fid_res,format,[m_NumElem;aux_var((npg-1)*siavarpg+2,:)]);
                fprintf(fid_res,'End Values\n');
            end
        case 8
            %Este modelo (conshyp==8) est� realizado para elementos triangulares lineales (donde se
            %utiliza un solo de punto de Gauss). Por ello en forma indiferente a la cantidad de puntos de
            %Gauss, se imprime el salto, normal, tangente en un solo punto de Gauss.
            
            %****************************************************************
            % Salto para el elemento de Fuerzas Centrales con fisuras de da�o
            %****************************************************************
            fprintf(fid_res,['Result "Modelo Constitutivo//Salto//Vector Salto" "Load Analysis" %d ',...
                'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            %Nombre de las componentes
            if ndime==2
                %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
                fprintf(fid_res,'ComponentNames "BetaX" "BetaY"\n');
            elseif ndime==3
                fprintf(fid_res,'ComponentNames "BetaX" "BetaY" "BetaZ"\n');
            end
            %Valores
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
            %caso 3D.
            format = ['%d',repmat(' %.15g',1,ndime),'\n'];
            %La componente X e Y del salto est� en los dos primeros elementos de la variable hvar_new.
            %Se est� asumiendo un solo punto de gauss utilizado en la integraci�n.
            fprintf(fid_res,format,[m_NumElem;hvar_new(1:ndime,:)]);
            fprintf(fid_res,'End Values \n');
            
            %*************************************************************************************
            % Vector Normal a la fisura para el elemento de Fuerzas Centrales con fisuras de da�o
            %*************************************************************************************
            fprintf(fid_res,['Result "Modelo Constitutivo//Salto//Vector Normal a la fisura" ',...
                '"Load Analysis" %d Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            if ndime==2
                %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
                fprintf(fid_res,'ComponentNames "nX" "nY"\n');
            elseif ndime==3
                fprintf(fid_res,'ComponentNames "nX" "nY" "nZ"\n');
            end
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
            %caso 3D.
            format = ['%d',repmat(' %.15g',1,ndime),'\n'];
            %La componente X e Y de la normal est� en la posici�n 7 y 8 de la variable hvar_new.
            %Se est� asumiendo un solo punto de gauss utilizado en la integraci�n.
            fprintf(fid_res,format,[m_NumElem;hvar_new(7:8,:)]);
            fprintf(fid_res,'End Values \n');
            
            %*************************************************************************************
            % Vector tangente a la fisura para el elemento de Fuerzas Centrales con fisuras de da�o
            %*************************************************************************************
            %Se asume la convenci�n la regla de la mano derecha para el producto t x n.
            %La componente X e Y de la normal est� en la posici�n 7 y 8 de la variable hvar_new.
            %Se est� asumiendo un solo punto de gauss utilizado en la integraci�n.
            m_VecTan = hvar_new(8:-1:7,:);
            m_VecTan(2,:) = -m_VecTan(2,:);
            %
            fprintf(fid_res,['Result "Modelo Constitutivo//Salto//Vector Tangente a la fisura" ',...
                '"Load Analysis" %d Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            if ndime==2
                %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
                fprintf(fid_res,'ComponentNames "tX" "tY"\n');
            elseif ndime==3
                fprintf(fid_res,'ComponentNames "tX" "tY" "tZ"\n');
            end
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
            %caso 3D.
            format = ['%d',repmat(' %.15g',1,ndime),'\n'];
            %
            fprintf(fid_res,format,[m_NumElem;m_VecTan]);
            fprintf(fid_res,'End Values \n');
            
            %*************************************************************************************
            % �ndice de carga FLoad para el elemento de Fuerzas Centrales con fisuras de da�o
            %*************************************************************************************
            fprintf(fid_res,['Result "Modelo Constitutivo//�ndice de Carga" "Load Analysis" %d Scalar ',...
                'OnGaussPoints "',nomGaussUnicoSet,'\n'],iStep);
            %Para utilizar la tabla de resultados previamente creada.
            fprintf(fid_res,['ResultRangesTable ','"Factor de Carga de fuerzas centradas"\n']);
            %Valores (al usarse un solo punto de gauss, se aplica un valor por elemento)
            fprintf(fid_res,'Values\n');
            fprintf(fid_res,'%d %d\n',[m_NumElem;hvar_new(9,:)]);
            fprintf(fid_res,'End Values\n');
            
            %*************************************************************************************
            % Tensi�n �ltima (M�dulo de vector de tensi�n) l�mite adoptada en cada elemento
            %*************************************************************************************
            fprintf(fid_res,['Result "Modelo Constitutivo//Tensi�n �ltima adoptada" "Load Analysis" ',...
                '%d Scalar OnGaussPoints "',nomGaussUnicoSet,'\n'],iStep);
            %Valores (al usarse un solo punto de gauss, se aplica un valor por elemento)
            fprintf(fid_res,'Values\n');
            fprintf(fid_res,'%d %d\n',[m_NumElem;hvar_new(4,:)]);
            fprintf(fid_res,'End Values\n');
            
        case 50
            
            %Se recupera las deformaciones de los dos PG, que fue cortada previamente para imprimir un
            %solo PG a nivel macro.
            f_PostProcME(iStep,m_NumElem,reshape(strain,ntens,npg,[]),e_DatMatSet,hvar_new,e_VG)
            
        case 51
            
            %          cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
            %          listeners = cmdWinDoc.getDocumentListeners;
            %          jFxCommandArea = listeners(3);
            %          set(jFxCommandArea,'Background','red');
            %Se recupera las deformaciones de los dos PG, que fue cortada previamente para imprimir un
            %solo PG a nivel macro.
            f_PostProcMECohesivo(iStep,m_NumElem,reshape(strain,ntens,npg,[]),e_DatMatSet,hvar_new,...
                aux_var,e_VG)
            %          set(jFxCommandArea,'Background','yellow');
            
        case 53
            
            %          cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
            %          listeners = cmdWinDoc.getDocumentListeners;
            %          jFxCommandArea = listeners(3);
            %          set(jFxCommandArea,'Background','red');
            %Se recupera las deformaciones de los dos PG, que fue cortada previamente para imprimir un
            %solo PG a nivel macro.
            f_PostProcMEBcna(iStep,m_NumElem,reshape(strain,ntens,npg,[]),e_DatMatSet,hvar_new,...
                aux_var,e_VG)
            %          set(jFxCommandArea,'Background','yellow');
            
        otherwise
            error('PostProceso: Modelo constitutivo no definido.')
    end
end

fclose(fid_res);