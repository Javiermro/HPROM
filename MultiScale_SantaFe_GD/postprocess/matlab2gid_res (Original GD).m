function matlab2gid_res(iStep,in,u,c_GdlCond,e_DatSet,e_VarEst,e_VarAux,DefMacro,uTotal,ELOC,e_VG)

%**************************************************
%* RUTINA PARA GENERAR LOS ARCHIVOS DE RESULTADOS *
%* PARA SER POSTPROCESADOS CON GID                *
%* Archivo: NAME.flavia.res                       *
%**************************************************

% Variables globales
ndime        = e_VG.ndime;
ndn = e_VG.ndn;
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

if e_VG.isPostResMultFile
   [dir,nomb] = fileparts(fileCompleto);
   fileCompleto = fullfile(dir,'GiDRes',nomb);
   %
   %Se guarda en el archivo List el nombre del archivo de los resultados de este paso tiempo.
   fileNameLst = [fileCompleto,'.post.lst'];
   fId = fopen(fileNameLst,'at');
   fprintf(fId,'%s\n',[nomb,'_P',num2str(iStep),'.flavia.res']);
   fclose(fId);
   %
   %Se abre archivo de resultados.
   filename_res = [fileCompleto,'_P',num2str(iStep),'.flavia.res'];
   fid_res = fopen(filename_res,'wt');
   fprintf(fid_res,'GiD Post Results File 1.0 \n');
else
   %Se abre archivo de resultados.
   filename_res = [fileCompleto,'.flavia.res'];
   fid_res = fopen(filename_res,'at');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if e_VG.esME
   % Fluctuacion de desplazamientos
   tipoDesp = '//Fluctuations';
else
   tipoDesp = '';
end
switch struhyp
   case {1,2,20}
fprintf(fid_res,['Result "Displacements',tipoDesp,'" "Load Analysis" %d Vector OnNodes\n'],iStep);
if ndime==2
   fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL"\n');
elseif ndime==3
   fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL" "Z-DISPL"\n');
end
   case 30
      fprintf(fid_res,['Result "Temperature',tipoDesp,'" "Load Analysis" %d Scalar OnNodes\n'],iStep);
   otherwise
      error('GidPost: Nodal Values: Structural hypothesis not defined.');
end
fprintf(fid_res,'Values\n');
format = ['%d',repmat(' %.15g',1,ndn),'\n'];
fprintf(fid_res,format,[in';reshape(u,ndn,[])]);   %[in,u(1:ndime:end),u(2:ndime:end)]'
fprintf(fid_res,'End Values\n');

if e_VG.esME
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Desplazamientos totales
   switch struhyp
      case {1,2,20}
   fprintf(fid_res,'Result "Displacements//Total" "Load Analysis" %d Vector OnNodes\n',iStep);
   if ndime==2
      fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL"\n');
   elseif ndime==3
      fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL" "Z-DISPL"\n');
   end
      case 30
         fprintf(fid_res,['Result "Temperature//Total" "Load Analysis" %d Scalar OnNodes\n'],iStep);
   otherwise
      error('GidPost: Multi-scale: Total Nodal Values: Structural hypothesis not defined.');
   end
   fprintf(fid_res,'Values\n');
   format = ['%d',repmat(' %.15g',1,ndn),'\n'];
   fprintf(fid_res,format,[in';uTotal]);
   fprintf(fid_res,'End Values\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Crack path field
if e_VG.exist_CrackPath
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Crack Path Field
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   fprintf(fid_res,'Result "Finite Element//SDA//Crack Path Field" "Load Analysis" %d Scalar OnNodes\n',iStep);
   fprintf(fid_res,'Values\n');
   format = ['%d', '  %.15g','\n'];
   fprintf(fid_res,format,[in'; e_VG.smooth_dalpha']);
   fprintf(fid_res,'End Values\n');
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % ALPHA Crack Path
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   fprintf(fid_res,'Result "Finite Element//SDA//ALPHA Crack Path" "Load Analysis" %d Scalar OnNodes\n',iStep);
   fprintf(fid_res,'Values\n');
   format = ['%d', '  %.15g','\n'];
   fprintf(fid_res,format,[in'; e_VG.smooth_alpha']);
   fprintf(fid_res,'End Values\n');
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % alpha_sin_suavizar
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for iSet = 1:nSet
      eltype = e_DatSet(iSet).e_DatElem.eltype;
      switch eltype
         case {21,22,23}
            i_vectVHElem   =  e_DatSet(iSet).e_DatElem.pointersVHE.i_vectVHElem ;
            nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
            nomGaussSet     = ['GP_Set_',num2str(iSet)];
            m_NumElem       = e_DatSet(iSet).m_NumElem;
            %p_varHistSmooth = 5;
            p_varHistSmooth= i_vectVHElem(4);
            m_VarHistElem   = e_VarEst(iSet).VarHistElem(p_varHistSmooth,:);
            %
            fprintf(fid_res,['Result "Finite Element//SDA//ALPHA without smoothing" "Load Analysis" %d ',...
               'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
            %Se considera ndime componentes, aunque no estï¿½ implementado para el caso 3D.
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
         case {21,22,23}
            i_vectVHElem   =  e_DatSet(iSet).e_DatElem.pointersVHE.i_vectVHElem;
            nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
            nomGaussSet     = ['GP_Set_',num2str(iSet)];
            m_NumElem       = e_DatSet(iSet).m_NumElem;
            %p_varHistSmooth = 2;  % disipacion
            p_varHistSmooth= i_vectVHElem(1);
            m_VarHistElem   = e_VarEst(iSet).VarHistElem(p_varHistSmooth,:);
            %
            fprintf(fid_res,['Result "Finite Element//SDA//Dissipation" "Load Analysis" %d ',...
               'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
            %Se considera ndime componentes, aunque no estï¿½ implementado para el caso 3D.
            format = ['%d',repmat(' %.15g',1,1),'\n'];
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
            fprintf(fid_res,'End Values \n');
      end
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % carga_descarga elementos
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for iSet = 1:nSet
      eltype = e_DatSet(iSet).e_DatElem.eltype;
      switch eltype
         case {21,22,23}
            i_vectVHElem   =  e_DatSet(iSet).e_DatElem.pointersVHE.i_vectVHElem;
            nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
            m_NumElem        = e_DatSet(iSet).m_NumElem;
            %p_varHistSmooth  = 2;
            p_varHistSmooth= i_vectVHElem(5);
            m_VarHistElem    = e_VarEst(iSet).VarHistElem(p_varHistSmooth,:);
            %
            fprintf(fid_res,['Result "Finite Element//SDA//Loading Elem. (Delta_r)" "Load Analysis" %d ',...
               'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
            %Se considera ndime componentes, aunque no estï¿½ implementado para el caso 3D.
            format = ['%d',repmat(' %.15g',1,1),'\n'];
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
            fprintf(fid_res,'End Values \n');
      end
   end
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % estado historico del elemento finito
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for iSet = 1:nSet
      eltype = e_DatSet(iSet).e_DatElem.eltype;
      switch eltype
         case {21,22,23}
            i_indST = e_DatSet(iSet).e_DatElem.pointersVHE.i_indST;
            nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
            m_NumElem        = e_DatSet(iSet).m_NumElem;
            m_VarHistElem    = e_VarEst(iSet).VarHistElem(i_indST,:);
            %
            fprintf(fid_res,['Result "Finite Element//SDA//Historical State" "Load Analysis" %d ',...
               'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
            %Se considera ndime componentes, aunque no estï¿½ implementado para el caso 3D.
            format = ['%d',repmat(' %.15g',1,1),'\n'];
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
            fprintf(fid_res,'End Values \n');
      end
   end
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % estado actual del elemento finito
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for iSet = 1:nSet
      eltype = e_DatSet(iSet).e_DatElem.eltype;
      switch eltype
         case {21,22,23}
            p_indActSTmacro = e_DatSet(iSet).e_DatElem.pointersVHE.p_indActSTmacro;
            nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
            m_NumElem        = e_DatSet(iSet).m_NumElem;
            m_VarHistElem    = e_VarEst(iSet).VarHistElem(p_indActSTmacro,:);
            %
            fprintf(fid_res,['Result "Finite Element//SDA//Current State" "Load Analysis" %d ',...
               'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
            %Se considera ndime componentes, aunque no estï¿½ implementado para el caso 3D.
            format = ['%d',repmat(' %.15g',1,1),'\n'];
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
            fprintf(fid_res,'End Values \n');
      end
   end
   
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % CPI in cutted elements
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   for iSet = 1:nSet
      eltype = e_DatSet(iSet).e_DatElem.eltype;
      %i_indST = e_DatSet(iSet).e_DatElem.pointersVHE.i_indST;
      switch eltype
         case {21,22,23}
            nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
            m_NumElem        = e_DatSet(iSet).m_NumElem;
            NumsidesCutCPF    = e_DatSet(iSet).e_DatElem.NumsidesCutCPF;
            %
            fprintf(fid_res,['Result "Finite Element//SDA//Current CPI" "Load Analysis" %d ',...
               'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
            %Se considera ndime componentes, aunque no estï¿½ implementado para el caso 3D.
            format = ['%d',repmat(' %.15g',1,1),'\n'];
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;NumsidesCutCPF']);
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
   % Variables internas segun el modelo constitutivo
   hvar_new = e_VarEst(iSet).hvar;
   
   %esImplex = e_DatMatSet.esImplex;   
   
   if((eltype== 21&&conshyp~=53)|| (eltype==22&&conshyp~=54)||(eltype== 23&&conshyp~=54) )
      hvar_new = reshape(hvar_new, sihvarpg,npg,[]);
      hvar_new = hvar_new(:,1:4,:);
      hvar_new = reshape(hvar_new,sihvarpg*4,[]);
   end
   
   aux_var = e_VarAux(iSet).VarAuxGP;
   m_VarAuxElem = e_VarAux(iSet).VarAuxElem;
   m_VarHistElem = e_VarEst(iSet).VarHistElem;
   
   nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
   nomGaussSet = ['GP_Set_',num2str(iSet)];
   
   switch eltype
      case {2,4,8,31,108}
         %
      case 32
         % Elementos localizados de banda
         if ~isempty(ELOC)
            m_ElemLocSet = ELOC(m_IndElemSet);
         else
            m_ElemLocSet = zeros(size(m_NumElem));
         end
         fprintf(fid_res,['Result "Strain//Localized Domain//',...
            'Sobre elemento usado" "Load Analysis" %d Scalar OnGaussPoints "',...
            nomGaussUnicoSet,'"\n'],iStep);
         fprintf(fid_res,'Values\n');
         format = '%d %d\n';
         fprintf(fid_res,format,[m_NumElem;m_ElemLocSet]);
         fprintf(fid_res,'End Values\n');
         %
         % VECTOR NORMAL del elemento finito de banda         
         fprintf(fid_res,['Result "Strain//Normal vector" "Load Analysis" %d ',...
            'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP"\n');
         fprintf(fid_res,'Values\n');
         %
         %No se imprime las normales micro que no pertenecen al dominio localizado.
         if ~isempty(ELOC)
            m_ElemLocSet = ELOC(m_IndElemSet);
            %Normal micro multiplicado por el sentido.
            m_NormMicro = bsxfun(@times,m_VarAuxElem(1,m_ElemLocSet),...
               e_DatElemSet.normal_micro(:,m_ElemLocSet));
            format = '%d %.15g %.15g\n';
            fprintf(fid_res,format,[m_NumElem(m_ElemLocSet);m_NormMicro]);
         end
         fprintf(fid_res,'End Values\n');
      case 10    % Triï¿½ngulo de 3 nodos con discontinuidades fuertes (SDA)
         %Como previo a la bifurcaciï¿½n y a la activaciï¿½n de la strong discontinuity (SD) el punto de
         %gauss singular no se actualiza (para ahorrar memoria) se actualiza en la etapa previa solo
         %para imprensiï¿½n, luego tendrï¿½a que seguir distintos caminos.
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
         fprintf(fid_res,['Result "Finite Element//SDA//Beta (Jump) vector" "Load Analysis" %d ',...
            'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         %Nombre de las componentes
         %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
         fprintf(fid_res,'ComponentNames "BetaX" "BetaY"\n');
         %Valores
         fprintf(fid_res,'Values\n');
         %Se considera ndime componentes, aunque no estï¿½ implementado para el caso 3D.
         format = ['%d',repmat(' %.15g',1,ndime),'\n'];
         m_Beta = c_GdlCond{iSet,1};
         %Se imprime la variable del elemento en un punto de gauss central.
         fprintf(fid_res,format,[m_NumElem;m_Beta]);
         fprintf(fid_res,'End Values \n');
         
         %********************************************
         % Vector Normal a la fisura en el elemento
         %********************************************
         fprintf(fid_res,['Result "Finite Element//SDA//Normal Vector to the crack" ',...
            '"Load Analysis" %d Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
         fprintf(fid_res,'ComponentNames "nX" "nY"\n');
         fprintf(fid_res,'Values\n');
         %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
         %caso 3D.
         format = ['%d',repmat(' %.15g',1,ndime),'\n'];
         m_VecNormal = m_VarAuxElem([2;7],:);
         %Se imprime la variable del elemento en un punto de gauss central.
         fprintf(fid_res,format,[m_NumElem;m_VecNormal]);
         fprintf(fid_res,'End Values \n');
         
         %**********************************************
         % Vector tangente a la fisura en el elemento
         %**********************************************
         %Se asume la convenciï¿½n la regla de la mano derecha para el producto t x n.
         %Se estï¿½ asumiendo un solo punto de gauss utilizado en la integraciï¿½n.
         m_VecTan = zeros(2,nElem);
         m_VecTan(1,:) = m_VecNormal(2,:);
         m_VecTan(2,:) = -m_VecNormal(1,:);
         %
         fprintf(fid_res,['Result "Finite Element//SDA//Tangent vector to the crack" ',...
            '"Load Analysis" %d Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
         fprintf(fid_res,'ComponentNames "tX" "tY"\n');
         fprintf(fid_res,'Values\n');
         %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
         %caso 3D.
         format = ['%d',repmat(' %.15g',1,ndime),'\n'];
         %
         fprintf(fid_res,format,[m_NumElem;m_VecTan]);
         fprintf(fid_res,'End Values \n');
         
         % Vector Tracciï¿½n
         fprintf(fid_res,['Result "Finite Element//SDA//Traction Vector" "Load Analysis" %d ',...
            'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         %Nombre de las componentes
         %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
         fprintf(fid_res,'ComponentNames "TX" "TY"\n');
         %Valores
         fprintf(fid_res,'Values\n');
         %Se considera ndime componentes, aunque no estï¿½ implementado para el caso 3D.
         format = ['%d',repmat(' %.15g',1,ndime),'\n'];
         m_VecTracc = m_VarAuxElem(19:20,:);
         %Se imprime la variable del elemento en un punto de gauss central.
         fprintf(fid_res,format,[m_NumElem;m_VecTracc]);
         fprintf(fid_res,'End Values \n');
      case 20   %Quadrï¿½ngulo de 4 nodos mixto con inyecciï¿½n de deformaciï¿½n
         % Condiciï¿½n de bifurcaciï¿½n o activaciï¿½n de la inyecciï¿½n (dominio mixto)
         m_CondBif = m_VarAuxElem(1,:);
         fprintf(fid_res,['Result "Finite Element//Q4MID//Mixed Domain" "Load Analysis" %d ',...
            'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         %Valores
         fprintf(fid_res,'Values\n');
         format = '%d %d\n';
         %Se imprime la variable del elemento en un punto de gauss central.
         fprintf(fid_res,format,[m_NumElem;m_CondBif]);
         fprintf(fid_res,'End Values \n');
         % Tensiï¿½n estabilizada en los PGs
         m_TensEstab = m_VarHistElem;
         switch struhyp
            case {1,2}
               %Deformaciï¿½n plana y Tensiï¿½n plana
               tipoDatAlmac = 'PlainDeformationMatrix';
               nomComponente = '"Stress XX" "Stress YY" "Stress XY" "Stress ZZ"';
            otherwise
               error(['GidPost: Tensiones en los puntos de Gauss: Componentes: No definidas para ',...
                  'esta hipï¿½tesis de carga']);
         end
         fprintf(fid_res,['Result "Finite Element//Q4MID//Stabilized Stresses" "Load Analysis" %d %s ',...
            'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
         fprintf(fid_res,'ComponentNames %s\n',nomComponente);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,npg)];
         if (struhyp==1||struhyp==2)
            % Deformaciï¿½n Plana y Tensiï¿½n Plana
            m_Ind = reshape((1:npg*ntens)',ntens,[]);
            m_Ind([ntens,ntens-1],:) = m_Ind([ntens-1,ntens],:);
            fprintf(fid_res,format,[m_NumElem;m_TensEstab(m_Ind(:),:)]);
         else
            fprintf(fid_res,format,[m_NumElem;m_TensEstab]);
         end
         fprintf(fid_res,'End Values\n');
      case {21,22,23}    % Quadrilatero 4 nodos con discontinuidades fuertes (SDA)
         NPG_ELE = npg-2;
         %*************************
         % Salto en el elemento
         %*************************
         fprintf(fid_res,['Result "Finite Element//SDA//Beta (Jump) vector" "Load Analysis" %d ',...
            'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         %Nombre de las componentes
         %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
         fprintf(fid_res,'ComponentNames "BetaX" "BetaY"\n');
         %Valores
         fprintf(fid_res,'Values\n');
         %Se considera ndime componentes, aunque no estï¿½ implementado para el caso 3D.
         format = ['%d',repmat(' %.15g',1,ndime),'\n'];
         m_Beta = c_GdlCond{iSet,1};
         %Se imprime la variable del elemento en un punto de gauss central.
         fprintf(fid_res,format,[m_NumElem;m_Beta]);
         fprintf(fid_res,'End Values \n');
         %********************************************
         % Vector Normal a la fisura en el elemento
         %********************************************
         p_n_tens = e_DatElemSet.pointersVAE.p_n_tens;
         m_n_tens = m_VarAuxElem(p_n_tens,:);
         fprintf(fid_res,['Result "Finite Element//SDA//Normal vector to the crack" ',...
            '"Load Analysis" %d Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
         fprintf(fid_res,'ComponentNames "nX" "nY"\n');
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat(' %.15g',1,ndime),'\n'];
         m_VecNormal = m_n_tens([1;4],:);
         fprintf(fid_res,format,[m_NumElem;m_VecNormal]);
         fprintf(fid_res,'End Values \n');
         %********************************************
         % Gradiente_Phi
         %********************************************
         p_phi_grad     =  e_DatElemSet.pointersVAE.p_phi_grad  ;
         m_phi_grad     =  m_VarAuxElem(p_phi_grad,:) ;
         m_VecNormal    = m_phi_grad([1;4],:);
         fprintf(fid_res,['Result "Finite Element//SDA//Gradiente_phi" ',...
            '"Load Analysis" %d Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
         fprintf(fid_res,'ComponentNames "nX" "nY"\n');
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat(' %.15g',1,ndime),'\n'];
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
         fprintf(fid_res,['Result "Finite Element//SDA//Normal vector" "Load Analysis" %d ',...
            'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         if ndime==2
            fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP"\n');
         elseif ndime==3
            fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP" "Z-COMP"\n');
         end
         fprintf(fid_res,'Values\n');
         
         %i_indST           =  e_DatElemSet.pointersVHE.i_indST    ;
         %m_indSTmacro       =  m_VarHistElem (i_indST        ,:)  ;
         %p_elem_BifType =  e_DatElemSet.pointersVAE.p_elem_Biftype ;
         %elem_BifType  =    m_VarAuxElem(p_elem_BifType,:)  ;
         p_indActSTmacro = e_DatElemSet.pointersVHE.p_indActSTmacro ;
         m_indActSTmacro = m_VarHistElem(p_indActSTmacro,:)  ;
         
         p_n_tens      = e_DatElemSet.pointersVAE.p_n_tens ;
         n_tens        =  m_VarAuxElem(p_n_tens    ,:)  ;
         n_tens         =  reshape(n_tens,4,2,[]);
         
         %indi_elem     = find (m_indSTmacro == 2);
         %indi_elem     = find (elem_BifType == 2);
         indi_elem     = find (m_indActSTmacro == 2);
         
         nx = n_tens(1,1,:);
         ny = n_tens(2,2,:);
         
         for iElem = 1:nElem
            if any(iElem==indi_elem)
               fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n',[m_NumElem(iElem),nx(iElem),ny(iElem),0]');
            else
               fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n',[m_NumElem(iElem),0,0,0]');
            end
         end
         fprintf(fid_res,'End Values  \n');      
         %
         %********************************************
         % Tensiones
         %********************************************
         % Se obtiene las tensiones de los primeros 4 PGs para imprimirlas.
         m_Ind = reshape((1:NPG_ELE*ntens)',ntens,[]);
         switch struhyp
            case {1,2} % Deformación Plana y Tensión Plana
               tipoDatAlmac = 'PlainDeformationMatrix';
               nomComponente = '"Stress XX" "Stress YY" "Stress XY" "Stress ZZ"';
               m_Ind([ntens,ntens-1],:) = m_Ind([ntens-1,ntens],:);
            case 3 %Tridimensional
               tipoDatAlmac = 'Matrix';
               nomComponente = '"Stress XX" "Stress YY" "Stress ZZ" "Stress XY" "Stress XZ" "Stress YZ"';
         end
         %
         %fprintf(fid_res,format,[m_NumElem;stress(m_Ind(:),:)]);
%       
%          STRESS = reshape(stress,ntens,npg,[]);
%          STRESS = STRESS(:,1:4,:);
%          STRESS = reshape(STRESS,ntens*NPG_ELE,[]);
         
         fprintf(fid_res,['Result "Stresses (Implex) Q1SDA//On Gauss Points" "Load Analysis" %d %s ',...
            'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
         fprintf(fid_res,'ComponentNames %s\n',nomComponente);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,NPG_ELE)];
         fprintf(fid_res,format,[m_NumElem;stress(m_Ind(:),:)]);
         fprintf(fid_res,'End Values\n');         
         
         i_stressTilde  =  e_DatElemSet.pointersVHE.i_stressTilde ;
         m_stressTilde  =   m_VarHistElem(i_stressTilde   ,:)  ;
         
%          STRESS= reshape(m_stressTilde,ntens,npg,[]);
%          STRESS=STRESS(:,1:4,:);
%          STRESS= reshape(STRESS,ntens*NPG_ELE,[]);
         
         fprintf(fid_res,['Result "Stresses (Tilde) Q1SDA//On Gauss Points" "Load Analysis" %d %s ',...
            'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
         fprintf(fid_res,'ComponentNames %s\n',nomComponente);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,NPG_ELE)];
         fprintf(fid_res,format,[m_NumElem;m_stressTilde(m_Ind(:),:)]);
         fprintf(fid_res,'End Values\n');
         
      otherwise         
         error('GidPost: PostProceso: Elemento Finito no definido.')
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if eltype~=[21,22,23]
      %
      switch struhyp
         %%SMALL DEFORMATIONS (Infinitesimal Strain)
         case {1,2} %Deformación plana y Tensión plana
            tipoDatAlmac = 'PlainDeformationMatrix';
            nomComponente = '"Stress XX" "Stress YY" "Stress XY" "Stress ZZ"';
            fprintf(fid_res,['Result "Stresses//On Gauss Points" "Load Analysis" %d %s ',...
               'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
            fprintf(fid_res,'ComponentNames %s\n',nomComponente);
         case 3 %Tridimensional
            tipoDatAlmac = 'Matrix';
            nomComponente = 'Stress XX" "Stress YY" "Stress ZZ" "Stress XY" "Stress XZ" "Stress YZ';
            fprintf(fid_res,['Result "Stresses//On Gauss Points" "Load Analysis" %d %s ',...
               'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
            fprintf(fid_res,'ComponentNames %s\n',nomComponente);
         %%LARGE DEFORMATIONS (Finite Strain)
         case 20  %Plane Deformations State (2D).            
            %tipoDatAlmac = 'PlainDeformationMatrix';
            %nomComponente = '"Stress XX" "Stress YY" "Stress XY" "Stress ZZ"';
            %Parece que no hay forma de imprimir tensores no simétricos en el GiD.
            %Para imprimir el tensor de tensiones Primero de Piola-Kirchhoff se utiliza la impresión
            %como escalar.
            fprintf(fid_res,['ResultGroup "Load Analysis" %d OnGaussPoints "',nomGaussSet,'"\n'],...
               iStep);
            fprintf(fid_res,'ResultDescription "Stresses//On Gauss Points//First Piola-Kirchhoff XX" scalar\n');
            fprintf(fid_res,'ResultDescription "Stresses//On Gauss Points//First Piola-Kirchhoff YY" scalar\n');
            fprintf(fid_res,'ResultDescription "Stresses//On Gauss Points//First Piola-Kirchhoff ZZ" scalar\n');
            %Recordar el se está guardando el tensor P transpuesto.
            fprintf(fid_res,'ResultDescription "Stresses//On Gauss Points//First Piola-Kirchhoff XY" scalar\n');
            fprintf(fid_res,'ResultDescription "Stresses//On Gauss Points//First Piola-Kirchhoff YX" scalar\n'); 
         otherwise
            error(['GidPost: Deformaciones en los puntos de Gauss: Componentes: ',...
               'Hipótesis estructural no implementada'])
      end
      fprintf(fid_res,'Values\n');
      format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,NPG_ELE)];
      switch struhyp
         case {1,2}
            m_Ind = reshape((1:npg*ntens)',ntens,[]);
            %Se intercambia las últimas filas, que en caso de Deformación Plana y Tensión Plana significa
            %cambiar componentes XY por ZZ.
            m_Ind([ntens,ntens-1],:) = m_Ind([ntens-1,ntens],:);
%             fprintf(fid_res,format,[m_NumElem;stress(m_Ind(:),:)]);
            fprintf(fid_res,format,[m_NumElem;stress(m_Ind(:,:))]); %JLM
         otherwise
            fprintf(fid_res,format,[m_NumElem;stress]);
      end
      fprintf(fid_res,'End Values\n');
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Fluctuacion de Deformaciones sobre los puntos de Gauss
      if e_VG.esME
         tipoDef = 'Total//';
      else
         tipoDef = '';
      end
      switch struhyp
         case {1,2}
            %Deformación plana y Tensión plana
            tipoDatAlmac = 'PlainDeformationMatrix';
            nomComponente = '"Strain XX" "Strain YY" "Strain XY" "Strain ZZ"';
            fprintf(fid_res,['Result "Strain//',tipoDef,'On Gauss Points" ',...
               '"Load Analysis" %d %s OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
            fprintf(fid_res,'ComponentNames %s\n',nomComponente);
         case 3
            %Tridimensional
            tipoDatAlmac = 'Matrix';
            nomComponente = ['Strain XX" "Strain YY" "Strain ZZ" "Strain XY" ',...
               '"Strain XZ" "Strain YZ'];
            %case 4
            %Axisimétrica
            %case 5
            %Barras 2D
            fprintf(fid_res,['Result "Strain//',tipoDef,'On Gauss Points" ',...
               '"Load Analysis" %d %s OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
            fprintf(fid_res,'ComponentNames %s\n',nomComponente);
         case 20
            fprintf(fid_res,['ResultGroup "Load Analysis" %d OnGaussPoints "',nomGaussSet,'"\n'],...
               iStep);
            tipoDefLD = ['Strain//',tipoDef,'On Gauss Points//'];
            fprintf(fid_res,['ResultDescription "',tipoDefLD,'Deformation Gradient XX" scalar\n']);
            fprintf(fid_res,['ResultDescription "',tipoDefLD,'Deformation Gradient YY" scalar\n']);
            fprintf(fid_res,['ResultDescription "',tipoDefLD,'Deformation Gradient ZZ" scalar\n']);
            fprintf(fid_res,['ResultDescription "',tipoDefLD,'Deformation Gradient XY" scalar\n']);
            fprintf(fid_res,['ResultDescription "',tipoDefLD,'Deformation Gradient YX" scalar\n']); 
         otherwise
            error(['GidPost: Deformaciones en los puntos de Gauss: Componentes: ',...
               'No definidas para esta hipótesis de carga']);
      end
      fprintf(fid_res,'Values\n');
      format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,npg)];
      switch struhyp
         case {1,2}
            m_Ind = reshape((1:npg*ntens)',ntens,[]);
            m_Ind([ntens,ntens-1],:) = m_Ind([ntens-1,ntens],:);
            fprintf(fid_res,format,[m_NumElem;strain(m_Ind(:),:)]);
         otherwise
            fprintf(fid_res,format,[m_NumElem;strain]);
      end
      fprintf(fid_res,'End Values\n');
      %
      %Por ahora no se imprime estos campos para el caso multiescala de grandes deformaciones
      if e_VG.esME&&struhyp~=20
         %
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Deformaciones totales sobre los puntos de Gauss
         fprintf(fid_res,['Result "Strain//Fluctuations//On Gauss Points" ',...
            '"Load Analysis" %d %s OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
         fprintf(fid_res,'ComponentNames %s\n',nomComponente);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,NPG_ELE)];
         if struhyp==1
            % Deformaciï¿½n Plana
            fprintf(fid_res,format,[m_NumElem;eps_fluct(m_Ind(:),:)]);
         else
            fprintf(fid_res,format,[m_NumElem;eps_fluct]);
         end
         fprintf(fid_res,'End Values\n');
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         % Deformación macro total aplicada en cada PG
         fprintf(fid_res,['Result "Strain//Macro applied//Over Element" ',...
            '"Load Analysis" %d %s OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep,tipoDatAlmac);
         fprintf(fid_res,'ComponentNames %s\n',nomComponente);
         fprintf(fid_res,'Values\n');
         %format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,npg)];
         format = ['%d',repmat(' %.15g',1,ntens),'\n'];
         if struhyp==1||struhyp==2
            m_DefMacro =reshape(DefMacro{iSet}([1,2,4,3],1,:),ntens,nElem); % JLM cambie : por 1
%             m_DefMacro = DefMacro{iSet}([1,2,4,3],1,:); %JLM
            % Deformación Plana y Tensiï¿½n Plana
            fprintf(fid_res,format,[m_NumElem;m_DefMacro]); %JLM
%                 m_DefMacro = reshape(DefMacro{iSet}([1,2,4,3],:,:),ntens*npg,nElem);
%                 % Deformaciï¿½n Plana y Tensiï¿½n Plana
%                 fprintf(fid_res,format,[m_NumElem;m_DefMacro]);
         else
            %m_DefMacro = reshape(DefMacro{iSet},ntens,nElem);
            m_DefMacro = DefMacro{iSet};
            fprintf(fid_res,format,[m_NumElem;m_DefMacro]);
         end
         fprintf(fid_res,'End Values\n');
      end
      %
      %             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %             % Condiciï¿½n de incremento de deformaciï¿½n total que se ha verificado
      %             fprintf(fid_res,['Result "Strain//Localized Domain//',...
      %                 'Sobre elemento usado" "Load Analysis" %d Scalar OnGaussPoints "',...
      %                 nomGaussUnicoSet,'"\n'],iStep);
      %             fprintf(fid_res,'Values\n');
      %             format = '%d %d\n';
      %             fprintf(fid_res,format,[m_NumElem;ELOC(1,m_IndElemSet)]);
      %             fprintf(fid_res,'End Values\n');
      %
      %             if size(ELOC,1)>1
      %                 % Condiciï¿½n de incremento de deformaciï¿½n total que se ha verificado
      %                 fprintf(fid_res,['Result "Strain//Localized Domain//',...
      %                     'Over element depending the strain increment" "Load Analysis" %d ',...
      %                     'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
      %                 fprintf(fid_res,'Values\n');
      %                 format = '%d %d\n';
      %                 fprintf(fid_res,format,[m_NumElem;ELOC(2,m_IndElemSet)]);
      %                 fprintf(fid_res,'End Values\n');
      %             end
   end
   
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   switch conshyp 
      case {1,52,100}
         %
      case {2,9}
         % Plasticidad J2
         % Deformacion plastica equivalente sobre puntos de Gauss
         fprintf(fid_res,['Result  "Constitutive Model//Plastic Equiv. Strain//',...
            'Sobre punto de Gauss" "Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
         fprintf(fid_res,format,[m_NumElem;hvar_new(ntens+1:sihvarpg:end,:)]);
         fprintf(fid_res,'End Values\n');
         
         % Índice de carga fload
         fprintf(fid_res,['Result  "Constitutive Model//Load Index//On Gauss Points" ',...
            '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat(' %d\n',1,NPG_ELE)];
         fload = hvar_new(ntens+2:sihvarpg:end,:);
         fprintf(fid_res,format,[m_NumElem;fload]);
         fprintf(fid_res,'End Values\n');
         
         % Índice de carga fload (por elemento)
         fprintf(fid_res,['Result  "Constitutive Model//Load Index//',...
            'On Element (any Fload PG)" "Load Analysis" %d Scalar OnGaussPoints "',...
            nomGaussUnicoSet,'"\n'],iStep);
         fprintf(fid_res,'Values\n');
         format = '%d %d\n';
         %Se adopta si un elemento tiene un PG con fload!=0, el elemento adopta un fload 1.
         fprintf(fid_res,format,[m_NumElem;any(fload)]);
         fprintf(fid_res,'End Values\n');
      case {11,13}  
         % Variable de daï¿½o
         fprintf(fid_res,['Result "Constitutive Model//Damage//On (1st-4th) Gauss Points" ',...
            '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
         fprintf(fid_res,format,[m_NumElem;hvar_new(1:sihvarpg:end,:)]);
         %fprintf(fid_res,format,[m_NumElem;aux_var(1:siavarpg:end,:)]);
         fprintf(fid_res,'End Values\n');
         
         
         if (eltype == 21 || eltype == 22|| eltype == 23 )
            nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
            %nomGaussSet     = ['GP_Set_',num2str(iSet)];
            m_NumElem       = e_DatSet(iSet).m_NumElem;
            p_varHistSmooth = sihvarpg*5+1;
            hvar = e_VarEst(iSet).hvar(p_varHistSmooth,:);
            %
            fprintf(fid_res,['Result "Finite Element//SDA//Damage PG(6)" "Load Analysis" %d ',...
               'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
            %Se considera ndime componentes, aunque no estï¿½ implementado para el caso 3D.
            format = ['%d',repmat(' %.15g',1,1),'\n'];
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;hvar]);
            fprintf(fid_res,'End Values \n');
            %
            %*************************************************************************************
            % Índice de carga Fload para el PG central 6
            %*************************************************************************************
            fprintf(fid_res,['Result "Constitutive Model//Load Index//PG6" ',...
               '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            %Para utilizar la tabla de resultados previamente creada.
            fprintf(fid_res,['ResultRangesTable ','"Factor de Carga de Daño Isotrópico"\n']);
            fprintf(fid_res,'Values\n');
            format = '%d %d\n';
            fprintf(fid_res,format,[m_NumElem;aux_var((6-1)*siavarpg+1,:)]);
            fprintf(fid_res,'End Values\n');
         end
         %end
         
         %*************************************************************************************
         % Índice de carga Fload
         %*************************************************************************************
         fprintf(fid_res,['Result "Constitutive Model//Load Index//On Gauss Points" ',...
            '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
         %Para utilizar la tabla de resultados previamente creada.
         fprintf(fid_res,['ResultRangesTable ','"Factor de Carga de Daño Isotrópico"\n']);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat(' %d\n',1,NPG_ELE)];
         fprintf(fid_res,format,[m_NumElem;aux_var(1:siavarpg:NPG_ELE*siavarpg,:)]);
         fprintf(fid_res,'End Values\n');
         
         % -------------------------------------
         
         % Incremento de variable interna (strain-like) del modelo de danio isotropo
         fprintf(fid_res,['Result "Constitutive Model//implicit delta_r//On (1st-4th) Gauss Points" ',...
            '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
         fprintf(fid_res,format,[m_NumElem;hvar_new(3:sihvarpg:end,:)]);
         %fprintf(fid_res,format,[m_NumElem;aux_var(1:siavarpg:end,:)]);
         fprintf(fid_res,'End Values\n');
         
         %if (eltype == 21)
         %    nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
         %    %nomGaussSet     = ['GP_Set_',num2str(iSet)];
         %    m_NumElem       = e_DatSet(iSet).m_NumElem;
         %    p_varHistSmooth = sihvarpg*5+3;
         %    hvar = e_VarEst(iSet).hvar(p_varHistSmooth,:);
         %    %
         %    fprintf(fid_res,['Result "Finite Element//SDA//implicit delta_r PG(6)" "Load Analysis" %d ',...
         %        'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         %    fprintf(fid_res,'Values\n');
         %    %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
         %    %Se considera ndime componentes, aunque no estï¿½ implementado para el caso 3D.
         %    format = ['%d',repmat(' %.15g',1,1),'\n'];
         %    %Se imprime la variable del elemento en un punto de gauss central.
         %    fprintf(fid_res,format,[m_NumElem;hvar]);
         %    fprintf(fid_res,'End Values \n');
         %end
         
         % Incremento de variable interna (strain-like) del modelo de danio isotropo
         fprintf(fid_res,['Result "Constitutive Model//implicit rn1//On (1st-4th) Gauss Points" ',...
            '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
         fprintf(fid_res,format,[m_NumElem;hvar_new(5:sihvarpg:end,:)]);
         %fprintf(fid_res,format,[m_NumElem;aux_var(1:siavarpg:end,:)]);
         fprintf(fid_res,'End Values\n');
         
         %if (eltype == 21)
         %    nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
         %    %nomGaussSet     = ['GP_Set_',num2str(iSet)];
         %    m_NumElem       = e_DatSet(iSet).m_NumElem;
         %    p_varHistSmooth = sihvarpg*5+5;
         %    hvar = e_VarEst(iSet).hvar(p_varHistSmooth,:);
         %    %
         %    fprintf(fid_res,['Result "Finite Element//SDA//implicit rn1" "Load Analysis" %d ',...
         %        'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         %    fprintf(fid_res,'Values\n');
         %    %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
         %    %Se considera ndime componentes, aunque no estï¿½ implementado para el caso 3D.
         %    format = ['%d',repmat(' %.15g',1,1),'\n'];
         %    %Se imprime la variable del elemento en un punto de gauss central.
         %    fprintf(fid_res,format,[m_NumElem;hvar]);
         %    fprintf(fid_res,'End Values \n');
         %end
         
         % -----------------------------------------------
         
 
      case {4,5,10,12} 
         % Daño isotropo estándar y regularizado
         % Variable de daï¿½o
         fprintf(fid_res,['Result "Constitutive Model//Damage//On Gauss Points" ',...
            '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
         %fprintf(fid_res,format,[m_NumElem;hvar_new(3:sihvarpg:end,:)]);
         fprintf(fid_res,format,[m_NumElem;aux_var(1:siavarpg:end,:)]);
         fprintf(fid_res,'End Values\n');
         %*************************************************************************************
         % Índice de carga FLoad
         %*************************************************************************************
         fprintf(fid_res,['Result "Constitutive Model//Load Index//On Gauss Points" ',...
            '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
         %Para utilizar la tabla de resultados previamente creada.
         fprintf(fid_res,['ResultRangesTable ','"Factor de Carga de Daño Isotrópico"\n']);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat(' %d\n',1,NPG_ELE)];
         fprintf(fid_res,format,[m_NumElem;aux_var(2:siavarpg:end,:)]);
         fprintf(fid_res,'End Values\n');
         
         %Impresiï¿½n de los valores centrales para el elemento Q1MixtoInjDef en forma separada
         if eltype==20
            % Variable de daï¿½o
            fprintf(fid_res,['Result "Constitutive Model//Damage//Punto de gauss Central" ',...
               '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            format = '%d %.15g\n';
            fprintf(fid_res,format,[m_NumElem;aux_var((npg-1)*siavarpg+1,:)]);
            fprintf(fid_res,'End Values\n');
            
            % Factor de carga (FLoad)
            fprintf(fid_res,['Result "Constitutive Model//Load Index//Punto de gauss Central" ',...
               '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            %Para utilizar la tabla de resultados previamente creada.
            fprintf(fid_res,['ResultRangesTable ','"Factor de Carga de Danio Isotropico"\n']);
            fprintf(fid_res,'Values\n');
            format = '%d %d\n';
            fprintf(fid_res,format,[m_NumElem;aux_var((npg-1)*siavarpg+2,:)]);
            fprintf(fid_res,'End Values\n');
         end
      case 8
         %Este modelo (conshyp==8) estï¿½ realizado para elementos triangulares lineales (donde se
         %utiliza un solo de punto de Gauss). Por ello en forma indiferente a la cantidad de puntos de
         %Gauss, se imprime el salto, normal, tangente en un solo punto de Gauss.
         
         %****************************************************************
         % Salto para el elemento de Fuerzas Centrales con fisuras de daï¿½o
         %****************************************************************
         fprintf(fid_res,['Result "Constitutive Model//Jump (Beta)//Jump (Beta) Vector" "Load Analysis" %d ',...
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
         %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
         %caso 3D.
         format = ['%d',repmat(' %.15g',1,ndime),'\n'];
         %La componente X e Y del salto estï¿½ en los dos primeros elementos de la variable hvar_new.
         %Se estï¿½ asumiendo un solo punto de gauss utilizado en la integraciï¿½n.
         fprintf(fid_res,format,[m_NumElem;hvar_new(1:ndime,:)]);
         fprintf(fid_res,'End Values \n');
         
         %*************************************************************************************
         % Vector Normal a la fisura para el elemento de Fuerzas Centrales con fisuras de daï¿½o
         %*************************************************************************************
         fprintf(fid_res,['Result "Constitutive Model//Jump (Beta)//Normal vector to the crack" ',...
            '"Load Analysis" %d Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         if ndime==2
            %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
            fprintf(fid_res,'ComponentNames "nX" "nY"\n');
         elseif ndime==3
            fprintf(fid_res,'ComponentNames "nX" "nY" "nZ"\n');
         end
         fprintf(fid_res,'Values\n');
         %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
         %caso 3D.
         format = ['%d',repmat(' %.15g',1,ndime),'\n'];
         %La componente X e Y de la normal estï¿½ en la posiciï¿½n 7 y 8 de la variable hvar_new.
         %Se estï¿½ asumiendo un solo punto de gauss utilizado en la integraciï¿½n.
         fprintf(fid_res,format,[m_NumElem;hvar_new(7:8,:)]);
         fprintf(fid_res,'End Values \n');
         
         %*************************************************************************************
         % Vector tangente a la fisura para el elemento de Fuerzas Centrales con fisuras de daï¿½o
         %*************************************************************************************
         %Se asume la convenciï¿½n la regla de la mano derecha para el producto t x n.
         %La componente X e Y de la normal estï¿½ en la posiciï¿½n 7 y 8 de la variable hvar_new.
         %Se estï¿½ asumiendo un solo punto de gauss utilizado en la integraciï¿½n.
         m_VecTan = hvar_new(8:-1:7,:);
         m_VecTan(2,:) = -m_VecTan(2,:);
         %
         fprintf(fid_res,['Result "Constitutive Model//Jump (Beta)//Tangent vector to the crack" ',...
            '"Load Analysis" %d Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
         if ndime==2
            %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
            fprintf(fid_res,'ComponentNames "tX" "tY"\n');
         elseif ndime==3
            fprintf(fid_res,'ComponentNames "tX" "tY" "tZ"\n');
         end
         fprintf(fid_res,'Values\n');
         %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
         %caso 3D.
         format = ['%d',repmat(' %.15g',1,ndime),'\n'];
         %
         fprintf(fid_res,format,[m_NumElem;m_VecTan]);
         fprintf(fid_res,'End Values \n');
         
         %*************************************************************************************
         % ï¿½ndice de carga FLoad para el elemento de Fuerzas Centrales con fisuras de daï¿½o
         %*************************************************************************************
         fprintf(fid_res,['Result "Constitutive Model//Load Index" "Load Analysis" %d Scalar ',...
            'OnGaussPoints "',nomGaussUnicoSet,'\n'],iStep);
         %Para utilizar la tabla de resultados previamente creada.
         fprintf(fid_res,['ResultRangesTable ','"Factor de Carga de fuerzas centradas"\n']);
         %Valores (al usarse un solo punto de gauss, se aplica un valor por elemento)
         fprintf(fid_res,'Values\n');
         fprintf(fid_res,'%d %d\n',[m_NumElem;hvar_new(9,:)]);
         fprintf(fid_res,'End Values\n');
         
         %*************************************************************************************
         % Tensiï¿½n ï¿½ltima (Mï¿½dulo de vector de tensiï¿½n) lï¿½mite adoptada en cada elemento
         %*************************************************************************************
         fprintf(fid_res,['Result "Constitutive Model//Adopted ultimate stress" "Load Analysis" ',...
            '%d Scalar OnGaussPoints "',nomGaussUnicoSet,'\n'],iStep);
         %Valores (al usarse un solo punto de gauss, se aplica un valor por elemento)
         fprintf(fid_res,'Values\n');
         fprintf(fid_res,'%d %d\n',[m_NumElem;hvar_new(4,:)]);
         fprintf(fid_res,'End Values\n');
         
      case {50,55}
         
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
         
      case {53,54}
         if (eltype == 21||eltype == 22||eltype == 23 )
            nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
            i_vectVHElem   =  e_DatSet(iSet).e_DatElem.pointersVHE.i_vectVHElem;
            m_NumElem       = e_DatSet(iSet).m_NumElem;
            p_varHistDamage = i_vectVHElem(8);
            VarHistElem = e_VarEst(iSet).VarHistElem(p_varHistDamage,:);
            %
            fprintf(fid_res,['Result "Finite Element//SDA//Damage PG(6)" "Load Analysis" %d ',...
               'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no estï¿½ implementado este modelo constitutivo para el
            %Se considera ndime componentes, aunque no estï¿½ implementado para el caso 3D.
            format = ['%d',repmat(' %.15g',1,1),'\n'];
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;VarHistElem]);
            fprintf(fid_res,'End Values \n');
         end
         f_PostProcMEBcna(iStep,m_NumElem,reshape(strain,ntens,npg,[]),e_DatMatSet,hvar_new,...
            aux_var,e_VG)
      %
      %Large Deformations
      case 110 % Plasticidad J2 con grandes deformaciones         
         % Deformacion plástica equivalente sobre puntos de Gauss
         fprintf(fid_res,['Result  "Constitutive Model//Plastic Equiv. Strain//',...
            'Sobre punto de Gauss" "Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
         fprintf(fid_res,format,[m_NumElem;hvar_new(1:sihvarpg:end,:)]);
         fprintf(fid_res,'End Values\n');
         
         % Deformacion plástica equivalente sobre puntos de Gauss
         fprintf(fid_res,['Result  "Constitutive Model//Plastic Equiv. Stress//',...
            'Sobre punto de Gauss" "Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
         fprintf(fid_res,format,[m_NumElem;hvar_new(2:sihvarpg:end,:)]);
         fprintf(fid_res,'End Values\n');
         
         % Tensor de Cauchy de deformación derecho plástico (simétrico)
         tipoDatAlmac = 'PlainDeformationMatrix';
         nomComponente = '"Cp XX" "Cp YY" "Cp XY" "Cp ZZ"';
         fprintf(fid_res,['Result "Constitutive Model//Plastic Right Cauchy-Green deform. tensor//',...
            'On Gauss Points" ',...
            '"Load Analysis" %d %s OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
         fprintf(fid_res,'ComponentNames %s\n',nomComponente);
         fprintf(fid_res,'Values\n');
         m_Ind = bsxfun(@plus,(3:3+3)',0:sihvarpg:(NPG_ELE-1)*sihvarpg);
         %Como el Cp viene con la convención de [Cpxx,Cpyy,Cpzz,Cpxy] se invierte para coincidir con el orden
         %del GiD.
         m_Ind([4,3],:) = m_Ind([3,4],:);
         format = ['%d',repmat([repmat(' %.15g',1,4),'\n'],1,NPG_ELE)];
         fprintf(fid_res,format,[m_NumElem;hvar_new(m_Ind(:),:)]);
         fprintf(fid_res,'End Values\n');         
         
         % Índice de carga fload
         fprintf(fid_res,['Result  "Constitutive Model//Load Index//On Gauss Points" ',...
            '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
         fprintf(fid_res,'Values\n');
         format = ['%d',repmat(' %d\n',1,NPG_ELE)];
         fload = aux_var(1:siavarpg:end,:);
         fprintf(fid_res,format,[m_NumElem;fload]);
         fprintf(fid_res,'End Values\n');
         
         if NPG_ELE>1
            % Índice de carga fload (por elemento)
            fprintf(fid_res,['Result  "Constitutive Model//Load Index//',...
               'On Element (any Fload PG)" "Load Analysis" %d Scalar OnGaussPoints "',...
               nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            format = '%d %d\n';
            %Se adopta si un elemento tiene un PG con fload!=0, el elemento adopta un fload 1.
            fprintf(fid_res,format,[m_NumElem;any(fload)]);
            fprintf(fid_res,'End Values\n');
         end
      otherwise
         error('PostProceso: Modelo constitutivo no definido.')
   end
end

fclose(fid_res);