function matlab2gid_PG_res(iStep,iStage,in,u,u_old,c_GdlCond,e_DatSet,e_VarEst,e_VarEst_old,e_VarAux,DefMacro,uTotal,ELOC,e_VG)

%**************************************************
%* RUTINA PARA GENERAR LOS ARCHIVOS DE RESULTADOS *
%* PARA SER POSTPROCESADOS CON GID                *
%* Archivo: NAME.flavia.res                       *
%* SOLO IMPRIME LOS PTOS. DE GAUSS SELECCIONANDOS *
%**************************************************

if ~isfield(e_VG,'MOD_TYPE');     return ; else    MOD_TYPE=e_VG.MOD_TYPE; end

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

% if (isfield(e_VG, 'nModesEPS_TOT')) && (isfield(e_VG, 'nModesSTR_REG')) && (isfield(e_VG, 'nModesSTR_DIS'))
%     % SECOND REDUCTION
%     filename_res = [fileCompleto,'_STAGE_',num2str(iStage),'_nMEPS_',num2str(e_VG.nModesEPS_TOT),...
%         '_nMStrREG_',num2str(e_VG.nModesSTR_REG),'_nMStrDIS_',num2str(e_VG.nModesSTR_DIS),'.flavia.res'];
% elseif (isfield(e_VG, 'nModesEPS_TOT')) && (e_VG.MOD_TYPE==1)
%     % FIRST REDUCTION
%     filename_res = [fileCompleto,'_STAGE_',num2str(iStage),'_nMEPS_',num2str(e_VG.nModesEPS_TOT),'.flavia.res'];
% else % HF & TRAINING
    filename_res = [fileCompleto,'_Selected_GP.flavia.res'];
% end

%filename_res = [fileCompleto,'_STAGE_',num2str(iStage),'.flavia.res'];
fid_res = fopen(filename_res,'at');

% To retrieve the original position (zero angle of rotation) for the small scale
if isfield(e_VG,'RVEAngle') && isfield(e_VG,'ScaleFACT')
    RVEAngle = e_VG.RVEAngle; %ScaleFACT = e_DatMatSet.e_VG.ScaleFACT;
    RotMatrix = [cos(RVEAngle) -sin(RVEAngle) 0 ; sin(RVEAngle) cos(RVEAngle) 0; 0 0 1];
    uROT=reshape(u,ndime,[])'*RotMatrix(1:ndime,1:ndime);
    uROT=uROT'; u=uROT(:); clear uROT
    
    u_oldROT=reshape(u_old,ndime,[])'*RotMatrix(1:ndime,1:ndime);
    u_oldROT=u_oldROT'; u_old=u_oldROT(:); clear u_oldROT
    
    uTotalROT = uTotal'*RotMatrix(1:ndime,1:ndime);
    %uTotalROT=uTotalROT'; uTotal=uTotalROT(:); clear uTotalROT
    uTotal=uTotalROT'; clear uTotalROT
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fluctuacion de desplazamientos
% fprintf(fid_res,'Result "Displacements//Fluctuations" "Load Analysis" %d Vector OnNodes\n',iStep);
% if ndime==2
%     fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL"\n');
% elseif ndime==3
%     fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL" "Z-DISPL"\n');
% end
% fprintf(fid_res,'Values\n');
% format = ['%d',repmat(' %.15g',1,ndime),'\n'];
% fprintf(fid_res,format,[in';[u(1:ndime:length(u)) u(2:ndime:length(u))]']);
% fprintf(fid_res,'End Values\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Desplazamientos totales
% fprintf(fid_res,'Result "Displacements//Total" "Load Analysis" %d Vector OnNodes\n',iStep);
% if ndime==2
%     fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL"\n');
% elseif ndime==3
%     fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL" "Z-DISPL"\n');
% end
% fprintf(fid_res,'Values\n');
% format = ['%d',repmat(' %.15g',1,ndime),'\n'];
% fprintf(fid_res,format,[in';uTotal]);
% fprintf(fid_res,'End Values\n');
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Incremento de Desplazamientos totales
% delta_u = u-u_old;
% fprintf(fid_res,'Result "Displacements//Incremental Fluct." "Load Analysis" %d Vector OnNodes\n',iStep);
% if ndime==2
%     fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL"\n');
% elseif ndime==3
%     fprintf(fid_res,'ComponentNames "X-DISPL" "Y-DISPL" "Z-DISPL"\n');
% end
% fprintf(fid_res,'Values\n');
% format = ['%d',repmat(' %.15g',1,ndime),'\n'];
% fprintf(fid_res,format,[in';[delta_u(1:ndime:length(u)) delta_u(2:ndime:length(u))]']);
% fprintf(fid_res,'End Values\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Crack path field
% if e_VG.exist_CrackPath
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Crack Path Field
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fprintf(fid_res,'Result "Finite Element//SDA//Crack Path Field" "Load Analysis" %d Scalar OnNodes\n',iStep);
%     fprintf(fid_res,'Values\n');
%     format = ['%d', '  %.15g','\n'];
%     fprintf(fid_res,format,[in'; e_VG.smooth_dalpha']);
%     fprintf(fid_res,'End Values\n');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % ALPHA Crack Path
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     fprintf(fid_res,'Result "Finite Element//SDA//ALPHA Crack Path" "Load Analysis" %d Scalar OnNodes\n',iStep);
%     fprintf(fid_res,'Values\n');
%     format = ['%d', '  %.15g','\n'];
%     fprintf(fid_res,format,[in'; e_VG.smooth_alpha']);
%     fprintf(fid_res,'End Values\n');
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % alpha_sin_suavizar
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for iSet = 1:nSet
%         eltype = e_DatSet(iSet).e_DatElem.eltype;
%         switch eltype
%             case {21,22}
%                 i_vectVHElem   =  e_DatSet(iSet).e_DatElem.pointersVHE.i_vectVHElem ;
%                 nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
%                 nomGaussSet     = ['GP_Set_',num2str(iSet)];
%                 m_NumElem       = e_DatSet(iSet).m_NumElem;
%                 %p_varHistSmooth = 5;
%                 p_varHistSmooth= i_vectVHElem(4);
%                 m_VarHistElem   = e_VarEst(iSet).VarHistElem(p_varHistSmooth,:);
%                 %
%                 fprintf(fid_res,['Result "Finite Element//SDA//ALPHA without smoothing" "Load Analysis" %d ',...
%                     'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
%                 fprintf(fid_res,'Values\n');
%                 %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
%                 %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
%                 format = ['%d',repmat(' %.15g',1,1),'\n'];
%                 %Se imprime la variable del elemento en un punto de gauss central.
%                 fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
%                 fprintf(fid_res,'End Values \n');
%         end
%     end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % dissipacion
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for iSet = 1:nSet
%         eltype = e_DatSet(iSet).e_DatElem.eltype;
%         switch eltype
%             case {21,22}
%                 i_vectVHElem   =  e_DatSet(iSet).e_DatElem.pointersVHE.i_vectVHElem;                       
%                 nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
%                 nomGaussSet     = ['GP_Set_',num2str(iSet)];
%                 m_NumElem       = e_DatSet(iSet).m_NumElem;
%                 %p_varHistSmooth = 2;  % disipacion
%                 p_varHistSmooth= i_vectVHElem(1);
%                 m_VarHistElem   = e_VarEst(iSet).VarHistElem(p_varHistSmooth,:);
%                 %
%                 fprintf(fid_res,['Result "Finite Element//SDA//Dissipation" "Load Analysis" %d ',...
%                     'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
%                 fprintf(fid_res,'Values\n');
%                 %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
%                 %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
%                 format = ['%d',repmat(' %.15g',1,1),'\n'];
%                 %Se imprime la variable del elemento en un punto de gauss central.
%                 fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
%                 fprintf(fid_res,'End Values \n');
%         end
%     end
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % rate de dissipacion
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for iSet = 1:nSet
%         eltype = e_DatSet(iSet).e_DatElem.eltype;
%         switch eltype
%             case {21,22}
%                 i_vectVHElem   =  e_DatSet(iSet).e_DatElem.pointersVHE.i_vectVHElem;                       
%                 nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
%                 nomGaussSet     = ['GP_Set_',num2str(iSet)];
%                 m_NumElem       = e_DatSet(iSet).m_NumElem;
%                 %p_varHistSmooth = 2;  % disipacion
%                 p_varHistSmooth= i_vectVHElem(10);
%                 m_VarHistElem   = e_VarEst(iSet).VarHistElem(p_varHistSmooth,:);
%                 %
%                 fprintf(fid_res,['Result "Finite Element//SDA//Rate of Dissipation" "Load Analysis" %d ',...
%                     'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
%                 fprintf(fid_res,'Values\n');
%                 %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
%                 %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
%                 format = ['%d',repmat(' %.15g',1,1),'\n'];
%                 %Se imprime la variable del elemento en un punto de gauss central.
%                 fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
%                 fprintf(fid_res,'End Values \n');
%         end
%     end    
%     
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % carga_descarga elementos
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for iSet = 1:nSet
%         eltype = e_DatSet(iSet).e_DatElem.eltype;
%         switch eltype
%             case {21,22}
%                 i_vectVHElem   =  e_DatSet(iSet).e_DatElem.pointersVHE.i_vectVHElem;                 
%                 nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
%                 m_NumElem        = e_DatSet(iSet).m_NumElem;
%                 %p_varHistSmooth  = 2;
%                 p_varHistSmooth= i_vectVHElem(5);
%                 m_VarHistElem    = e_VarEst(iSet).VarHistElem(p_varHistSmooth,:);
%                 %
%                 fprintf(fid_res,['Result "Finite Element//SDA//Loading Elem. (Delta_r)" "Load Analysis" %d ',...
%                     'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
%                 fprintf(fid_res,'Values\n');
%                 %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
%                 %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
%                 format = ['%d',repmat(' %.15g',1,1),'\n'];
%                 %Se imprime la variable del elemento en un punto de gauss central.
%                 fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
%                 fprintf(fid_res,'End Values \n');
%         end
%     end
% 
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % estado historico del elemento finito
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for iSet = 1:nSet
%         eltype = e_DatSet(iSet).e_DatElem.eltype;        
%         switch eltype
%             case {21,22}
%                 i_indST = e_DatSet(iSet).e_DatElem.pointersVHE.i_indST;                
%                 nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
%                 m_NumElem        = e_DatSet(iSet).m_NumElem;
%                 m_VarHistElem    = e_VarEst(iSet).VarHistElem(i_indST,:);
%                 %
%                 fprintf(fid_res,['Result "Finite Element//SDA//Historical State" "Load Analysis" %d ',...
%                     'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
%                 fprintf(fid_res,'Values\n');
%                 %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
%                 %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
%                 format = ['%d',repmat(' %.15g',1,1),'\n'];
%                 %Se imprime la variable del elemento en un punto de gauss central.
%                 fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
%                 fprintf(fid_res,'End Values \n');
%         end
%     end    
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % estado actual del elemento finito
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for iSet = 1:nSet
%         eltype = e_DatSet(iSet).e_DatElem.eltype;    
%         switch eltype
%             case {21,22}
%                 p_indActSTmacro = e_DatSet(iSet).e_DatElem.pointersVHE.p_indActSTmacro;                    
%                 nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
%                 m_NumElem        = e_DatSet(iSet).m_NumElem;
%                 m_VarHistElem    = e_VarEst(iSet).VarHistElem(p_indActSTmacro,:);
%                 %
%                 fprintf(fid_res,['Result "Finite Element//SDA//Current State" "Load Analysis" %d ',...
%                     'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
%                 fprintf(fid_res,'Values\n');
%                 %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
%                 %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
%                 format = ['%d',repmat(' %.15g',1,1),'\n'];
%                 %Se imprime la variable del elemento en un punto de gauss central.
%                 fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
%                 fprintf(fid_res,'End Values \n');
%         end
%     end
%     
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % CPI in cutted elements
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for iSet = 1:nSet
%         eltype = e_DatSet(iSet).e_DatElem.eltype;
%         %i_indST = e_DatSet(iSet).e_DatElem.pointersVHE.i_indST;        
%         switch eltype
%             case {21,22}
%                 nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
%                 m_NumElem        = e_DatSet(iSet).m_NumElem;                
%                 NumsidesCutCPF    = e_DatSet(iSet).e_DatElem.NumsidesCutCPF;
%                 %
%                 fprintf(fid_res,['Result "Finite Element//SDA//Current CPI" "Load Analysis" %d ',...
%                     'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
%                 fprintf(fid_res,'Values\n');
%                 %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
%                 %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
%                 format = ['%d',repmat(' %.15g',1,1),'\n'];
%                 %Se imprime la variable del elemento en un punto de gauss central.
%                 fprintf(fid_res,format,[m_NumElem;NumsidesCutCPF']);
%                 fprintf(fid_res,'End Values \n');
%         end
%     end 
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % longitud de fractura calculada
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for iSet = 1:nSet
%         eltype = e_DatSet(iSet).e_DatElem.eltype;
%         switch eltype
%             case 21
%                 i_vectVHElem   =  e_DatSet(iSet).e_DatElem.pointersVHE.i_vectVHElem;                       
%                 nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
%                 nomGaussSet     = ['GP_Set_',num2str(iSet)];
%                 m_NumElem       = e_DatSet(iSet).m_NumElem;
%                 %p_varHistSmooth = 2;  % disipacion
%                 %p_varHistSmooth= i_vectVHElem(6);
%                 p_varHistSmooth= i_vectVHElem(7);
%                 m_VarHistElem   = e_VarEst(iSet).VarHistElem(p_varHistSmooth,:);
%                 %
%                 fprintf(fid_res,['Result "Finite Element//SDA//Elem. Fract. Length Smu" "Load Analysis" %d ',...
%                     'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
%                 fprintf(fid_res,'Values\n');
%                 %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
%                 %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
%                 format = ['%d',repmat(' %.15g',1,1),'\n'];
%                 %Se imprime la variable del elemento en un punto de gauss central.
%                 fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
%                 fprintf(fid_res,'End Values \n');
%         end
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % variable de regularizacion de la cinematica macro usando los
%     % parametros automaticos en la escala baja
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for iSet = 1:nSet
%         eltype = e_DatSet(iSet).e_DatElem.eltype;
%         switch eltype
%             case 21
%                 i_vectVHElem   =  e_DatSet(iSet).e_DatElem.pointersVHE.i_vectVHElem;                       
%                 nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
%                 nomGaussSet     = ['GP_Set_',num2str(iSet)];
%                 m_NumElem       = e_DatSet(iSet).m_NumElem;
%                 %p_varHistSmooth = 2;  % disipacion
%                 p_varHistlength= i_vectVHElem(6);
%                 m_VarHistElem   = e_VarEst(iSet).VarHistElem(p_varHistlength,:);
%                 %
%                 fprintf(fid_res,['Result "Finite Element//SDA//Regulariz. length Lmu" "Load Analysis" %d ',...
%                     'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
%                 fprintf(fid_res,'Values\n');
%                 %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
%                 %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
%                 format = ['%d',repmat(' %.15g',1,1),'\n'];
%                 %Se imprime la variable del elemento en un punto de gauss central.
%                 fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
%                 fprintf(fid_res,'End Values \n');
%         end
%     end      
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % VARIABLE DE SUAVIZADO EN DEFS Y TENSIONES
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     for iSet = 1:nSet
%         eltype = e_DatSet(iSet).e_DatElem.eltype;
%         switch eltype
%             case {21,22}
%                 i_vectVHElem   =  e_DatSet(iSet).e_DatElem.pointersVHE.i_vectVHElem;                       
%                 nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
%                 nomGaussSet     = ['GP_Set_',num2str(iSet)];
%                 m_NumElem       = e_DatSet(iSet).m_NumElem;
%                 %p_varHistSmooth = 2;  % disipacion
%                 p_varHistSmoothNEW= i_vectVHElem(13);
%                 m_VarHistElem   = e_VarEst(iSet).VarHistElem(p_varHistSmoothNEW,:);
%                 %
%                 fprintf(fid_res,['Result "Finite Element//SDA//norm_eps_INVc_sigma" "Load Analysis" %d ',...
%                     'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
%                 fprintf(fid_res,'Values\n');
%                 %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
%                 %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
%                 format = ['%d',repmat(' %.15g',1,1),'\n'];
%                 %Se imprime la variable del elemento en un punto de gauss central.
%                 fprintf(fid_res,format,[m_NumElem;m_VarHistElem]);
%                 fprintf(fid_res,'End Values \n');
%         end
%     end    
%     
% end

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
    
%     hvar_new = e_VarEst(iSet).hvar; %JLM agregue yo :P
        
%     %if ~(MOD_TYPE==2 && (e_VG.esME==1))
%     if ~(MOD_TYPE==2) % && (e_VG.esME==1))
%         stress = e_VarEst(iSet).sigma;
%         strain = e_VarEst(iSet).eps;
%         %strain_fluct = e_VarEst(iSet).eps_fluct;
% 
%         eps_fluct = e_VarEst(iSet).eps_fluct;
% 
%         % Variables internas seg�n el modelo constitutivo
%         hvar_new = e_VarEst(iSet).hvar;
%         aux_var = e_VarAux(iSet).VarAuxGP;   
%         
%         % for computing incremental displacement at the small
%         hvar_old = e_VarEst_old.hvar;        
%         
%         if (eltype== 21) && (conshyp~=53)
%             hvar_new =reshape(hvar_new, sihvarpg,npg,[]);
%             hvar_new = hvar_new(:,1:4,:);
%             hvar_new =reshape(hvar_new,sihvarpg*4,[]);
%             
%             aux_var = reshape(aux_var, siavarpg,npg,[]);
%             aux_var = aux_var(:,1:4,:);
%             aux_var = reshape(aux_var,siavarpg*4,[]);
%             
%         end
%         
%         m_VarAuxElem = e_VarAux(iSet).VarAuxElem;
%         m_VarHistElem = e_VarEst(iSet).VarHistElem;
%          
%     end
    
    nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
    nomGaussSet = ['GP_Set_',num2str(iSet)];    

    NPG_ELE = npg;
    switch eltype
       case {4,8,31,32,108} %JLM agrego el 108
                   
               switch MOD_TYPE
                   case {0,1} % (HIGH FIDELITY & FIRST REDUCTION)
                       switch struhyp
                           case {1,2}
                               %Deformacion plana y Tension plana
                               tipoDatAlmac = 'PlainDeformationMatrix';
                               nomComponente = '"Stress XX" "Stress YY" "Stress ZZ" "Stress XY"';
                           otherwise
                               error(['GidPost: Tensiones en los puntos de Gauss: Componentes: No definidas para ',...
                                   'esta hip�tesis de carga']);
                       end
                       
                       if isfield(e_VG,'RVEAngle') && isfield(e_VG,'ScaleFACT')
                           stressROT = zeros(size(stress,1),size(stress,2));
                           for iElem = 1:nElem
                               stressiPG = reshape(stress(:,iElem),ntens,npg);
                               %stressiPG(ntens,:) = 0.5*stressiPG(ntens,:) ;
                               CCC    = arrayfun(@(x)RotMatrix(1:ndime,1:ndime)'*[stressiPG(1,x) stressiPG(4,x);stressiPG(4,x) stressiPG(2,x)]*RotMatrix(1:ndime,1:ndime),1:npg,'UniformOutput',false) ;
                               %CCCaux = arrayfun(@(x)CCC{x}(:),1:npg,'UniformOutput',false)' ;
                               CCCaux = arrayfun(@(x)[CCC{x}(1) CCC{x}(4) stressiPG(3,x) CCC{x}(2)]',1:npg,'UniformOutput',false)' ;
                               stressROT(:,iElem) = cell2mat(CCCaux);
                           end
                           stress=stressROT;
                       end
                       
                       fprintf(fid_res,['Result "Stresses//On Gauss Points" "Load Analysis" %d %s ',...
                           'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
                       fprintf(fid_res,'ComponentNames %s\n',nomComponente);
                       fprintf(fid_res,'Values\n');
                       format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,npg)];
                       fprintf(fid_res,format,[m_NumElem;stress]);
                       fprintf(fid_res,'End Values\n');
                       
                       switch struhyp
                           case {1,2}
                               %Deformacion plana y Tension plana
                               tipoDatAlmac = 'PlainDeformationMatrix';
                               %nomComponente = '"Strain XX" "Strain YY" "Strain XY" "Strain ZZ"';
                               nomComponente = '"Strain XX" "Strain YY" "Strain ZZ" "Strain XY"';
                           otherwise
                               error(['GidPost: Tensiones en los puntos de Gauss: Componentes: No definidas para ',...
                                   'esta hip�tesis de carga']);
                       end
                       
                       %Previous treatment for rotating strains (small scale)
                       if isfield(e_VG,'RVEAngle') && isfield(e_VG,'ScaleFACT')
                           eps_fluctROT = zeros(size(eps_fluct,1),size(eps_fluct,2));
                           for iElem = 1:nElem
                               epsiPG = reshape(eps_fluct(:,iElem),ntens,npg);
                               epsiPG(ntens,:) = 0.5*epsiPG(ntens,:) ;
                               CCC    = arrayfun(@(x)RotMatrix(1:ndime,1:ndime)'*[epsiPG(1,x) epsiPG(4,x);epsiPG(4,x) epsiPG(2,x)]*RotMatrix(1:ndime,1:ndime),1:npg,'UniformOutput',false) ;
                               %CCCaux = arrayfun(@(x)CCC{x}(:),1:npg,'UniformOutput',false)' ;
                               %CCCaux = arrayfun(@(x)[CCC{x}(1) CCC{x}(4) epsiPG(3,x) 2*CCC{x}(2)]',1:npg,'UniformOutput',false)' ;
                               CCCaux = arrayfun(@(x)[CCC{x}(1) CCC{x}(4) epsiPG(3,x) 2*CCC{x}(2)]',1:npg,'UniformOutput',false)' ;
                               eps_fluctROT(:,iElem) = cell2mat(CCCaux);
                           end
                           %eps_fluct=eps_fluctROT;
                           %%% TMP=reshape(eps_fluctROT,ntens,npg,[]);
                           %%% TMP(:,:,:)=TMP([1 2 4 3]',:,:);
                           %%% eps_fluctROT = reshape(TMP,ntens*npg,[]);
                           eps_fluct=eps_fluctROT;
                       end
                       
                       fprintf(fid_res,['Result "StrainsFluct//On Gauss Points" "Load Analysis" %d %s ',...
                           'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
                       fprintf(fid_res,'ComponentNames %s\n',nomComponente);
                       fprintf(fid_res,'Values\n');
                       format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,npg)];
                       fprintf(fid_res,format,[m_NumElem;eps_fluct]);
                       fprintf(fid_res,'End Values\n');
                       
               end           
            
            if (eltype ==32)
                %***************** VECTOR ELEMENTAL *****************
                fprintf(fid_res,['Result "Finite Element//SDA//Normal_Elem" "Load Analysis" %d ',...
                    'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
                if ndime==2
                    fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP"\n');
                elseif ndime==3
                    fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP" "Z-COMP"\n');
                end
                fprintf(fid_res,'Values\n');
                
                N_ELEM_Set = e_DatElemSet.N_ELEM ;
                N_ELEM_Set_x = N_ELEM_Set(1,:);
                N_ELEM_Set_y = N_ELEM_Set(2,:);
                
                for iElem = 1:nElem
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),N_ELEM_Set_x(iElem),N_ELEM_Set_y(iElem),0]');
                end
                fprintf(fid_res,'End Values  \n');
            end
            
            % IMPRESION DE PUNTOS DE GAUSS SELECCIONADOS POR LA SEGUNDA REDUCCION
            %if (e_VG.esME==1) && isfield(e_VG,'MOD_TYPE')
            if isfield(e_VG,'MOD_TYPE') % In order to include problems with and without multiscale
                switch (e_VG.MOD_TYPE)
                    case {0.5,1} % (HIGH FIDELITY & FIRST REDUCTION)
                        % Energy in each gauss point at the small (ROM_MODEL)
                        i_pIntEnergy = e_DatElemSet.pointersVHE.i_pIntEnergy;
                        IntEnergy = e_VarEst(iSet).VarHistElem(i_pIntEnergy ,:) ;   
                        fprintf(fid_res,['Result "Energy//GPs//On (1st-4th) Gauss Points" ',...
                            '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
                        fprintf(fid_res,'Values\n');
                        format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
                        fprintf(fid_res,format,[m_NumElem;IntEnergy]);
                        fprintf(fid_res,'End Values\n');
                    case 2 % (SECOND REDUCTION)
                        %if e_VG.MOD_TYPE==2
                        GaussFLAG = e_VG.GaussFLAG(m_NumElem,:);
                        fprintf(fid_res,['Result "Selected_GPs//GPs//On (1st-4th) Gauss Points" ',...
                            '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
                        fprintf(fid_res,'Values\n');
                        format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
                        %fprintf(fid_res,format,[m_NumElem;hvar_new(1:sihvarpg:end,:)]);
                        fprintf(fid_res,format,[m_NumElem;GaussFLAG(:,2:end)']);
                        fprintf(fid_res,'End Values\n');
                        %end
                end
            end
 
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
            fprintf(fid_res,['Result "Finite Element//SDA//Beta (Jump) vector" "Load Analysis" %d ',...
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
            fprintf(fid_res,['Result "Finite Element//SDA//Normal Vector to the crack" ',...
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
            fprintf(fid_res,['Result "Finite Element//SDA//Tangent vector to the crack" ',...
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
            fprintf(fid_res,['Result "Finite Element//SDA//Traction Vector" "Load Analysis" %d ',...
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
            fprintf(fid_res,['Result "Finite Element//Q4MID//Mixed Domain" "Load Analysis" %d ',...
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
                    nomComponente = '"Stress XX" "Stress YY" "Stress XY" "Stress ZZ"';
                otherwise
                    error(['GidPost: Tensiones en los puntos de Gauss: Componentes: No definidas para ',...
                        'esta hip�tesis de carga']);
            end
            fprintf(fid_res,['Result "Finite Element//Q4MID//Stabilized Stresses" "Load Analysis" %d %s ',...
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
            NPG_ELE = npg;
            %p_n_tens       =  e_DatElemSet.pointersVAE.p_n_tens    ;
            %p_elem_BifType =  e_DatElemSet.pointersVAE.p_elem_Biftype ;
            
            %m_elem_BifType  =  m_VarAuxElem (p_elem_BifType , : ) ;
            %m_n_tens        =  m_VarAuxElem (p_n_tens  , : ) ;
            %******************************
            % Indice de estado del elemento
            %******************************
            %fprintf(fid_res,['Result "Finite Element//SDA//State Index" "Load Analysis" %d ',...
            %    'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            %fprintf(fid_res,'Values\n');
            %format = ['%d',' %.15g','\n'];
            %fprintf(fid_res,format,[m_NumElem; m_elem_BifType]);
            %fprintf(fid_res,'End Values \n');
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
            %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
            format = ['%d',repmat(' %.15g',1,ndime),'\n'];
            m_Beta = c_GdlCond{iSet,1};
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;m_Beta]);
            fprintf(fid_res,'End Values \n');
            %*******************************************************
            % Valor absoluto del salto total para el elemento finito
            %*******************************************************
            fprintf(fid_res,['Result "Finite Element//SDA//Norm (Beta) " "Load Analysis" %d ',...
                'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            %Nombre de las componentes
            %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
            %fprintf(fid_res,'ComponentNames "NORM_BetaX" "NORM_BetaY"\n');
            %Valores
            fprintf(fid_res,'Values\n');
            %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
            format = ['%d',repmat(' %.15g',1,1),'\n'];
            m_NORMBeta = c_GdlCond{iSet,12};
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;m_NORMBeta]);
            fprintf(fid_res,'End Values \n');            
            
            %*******************************
            % SUAVIZADO - Derivada de alpha
            %*******************************
            fprintf(fid_res,['Result "Finite Element//SDA//d_alpha" "Load Analysis" %d ',...
                'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            p_dalpha  = e_DatElemSet.pointersVAE.p_dalpha ;
            dalpha = m_VarAuxElem(p_dalpha , :);  
            format = ['%d',repmat(' %.15g',1,1),'\n'];
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;dalpha]);
            fprintf(fid_res,'End Values \n');     
            
            %*******************
            % BIFURCATION FLAG -   
            %*******************
            fprintf(fid_res,['Result "Finite Element//SDA//BIF_FLAG" "Load Analysis" %d ',...
                'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            p_condBif     = e_DatElemSet.pointersVAE.p_condBif ;            
            %p_dalpha  = e_DatElemSet.pointersVAE.p_dalpha ;
            condBif = m_VarAuxElem(p_condBif , :);  
            format = ['%d',repmat(' %.15g',1,1),'\n'];
            %Se imprime la variable del elemento en un punto de gauss central.
            fprintf(fid_res,format,[m_NumElem;condBif]);
            fprintf(fid_res,'End Values \n');                 
            
            %********************************************
            % Vector Normal a la fisura en el elemento
            %********************************************
            % fprintf(fid_res,['Result "Finite Element//SDA//Normal vector to the crack" ',...
            %     '"Load Analysis" %d Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            % %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
            % fprintf(fid_res,'ComponentNames "nX" "nY"\n');
            % fprintf(fid_res,'Values\n');
            % %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
            % %caso 3D.
            % format = ['%d',repmat(' %.15g',1,ndime),'\n'];
            % m_VecNormal = m_n_tens([1;4],:);
            % %Se imprime la variable del elemento en un punto de gauss central.
            % fprintf(fid_res,format,[m_NumElem;m_VecNormal]);
            % fprintf(fid_res,'End Values \n');

            
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
            
            
            %********************************************
            % VECTOR NORMAL N DE ANALISIS DE BIFURCACION 
            %********************************************
            fprintf(fid_res,['Result "Finite Element//SDA//Normal vector (Bif. Analysis)" "Load Analysis" %d ',...
                'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            if ndime==2
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP"\n');
            elseif ndime==3
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP" "Z-COMP"\n');
            end
            fprintf(fid_res,'Values\n');
            
            p_indActSTmacro = e_DatElemSet.pointersVHE.p_indActSTmacro ;
            m_indActSTmacro = m_VarHistElem(p_indActSTmacro,:)  ; 
            
            p_normal_bif = e_DatElemSet.pointersVAE.p_normal_bif ;
            nBIF_tens    =  m_VarAuxElem(p_normal_bif(1:2),:)  ;

            indi_elem    = find (m_indActSTmacro==1|m_indActSTmacro==2);
            
            nBIFx = nBIF_tens(1,:);
            nBIFy = nBIF_tens(2,:);            
            for iElem = 1:nElem
                if find(iElem==indi_elem)
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),nBIFx(iElem),nBIFy(iElem),0]');
                else
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),0,0,0]');
                end
            end
            fprintf(fid_res,'End Values  \n');
            % *********************************************
            

            %********************************************
            % VECTOR NORMAL M DE ANALISIS DE BIFURCACION 
            %********************************************
            fprintf(fid_res,['Result "Finite Element//SDA//M vector (Bif. Analysis)" "Load Analysis" %d ',...
                'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            if ndime==2
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP"\n');
            elseif ndime==3
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP" "Z-COMP"\n');
            end
            fprintf(fid_res,'Values\n');
            
            p_indActSTmacro = e_DatElemSet.pointersVHE.p_indActSTmacro ;
            m_indActSTmacro = m_VarHistElem(p_indActSTmacro,:)  ; 
            
            p_normal_bif = e_DatElemSet.pointersVAE.p_normal_bif ;
            mBIF_tens    =  m_VarAuxElem(p_normal_bif(3:4),:)  ;

            indi_elem    = find (m_indActSTmacro==1|m_indActSTmacro==2);
            
            mBIFx = mBIF_tens(1,:);
            mBIFy = mBIF_tens(2,:);            
            for iElem = 1:nElem
                if find(iElem==indi_elem)
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),mBIFx(iElem),mBIFy(iElem),0]');
                else
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),0,0,0]');
                end
            end
            fprintf(fid_res,'End Values  \n');
            % *********************************************            
            
            
            %*****************************************************************************************************
            % VECTOR NORMAL N DE ANALISIS DE BIFURCACION del elemento finito (VECTOR DE INYECCION DE LA DISCONTINUIDAD)
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
            %indi_elem     = find (m_indActSTmacro == 2);
            
            % La normal la postprocesara desde el instante en el que ha bifurcado
            indi_elem     = find (m_indActSTmacro==1|m_indActSTmacro==2);
            
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

            %*****************************************************************************************************
            % VECTOR NORMAL M DE ANALISIS DE BIFURCACION del elemento finito (SEGUNDO VECTOR DEL ANALISIS DE BIFURCACION)
            fprintf(fid_res,['Result "Finite Element//SDA//M vector" "Load Analysis" %d ',...
                'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            if ndime==2
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP"\n');
            elseif ndime==3
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP" "Z-COMP"\n');
            end
            fprintf(fid_res,'Values\n');
               
            p_indActSTmacro = e_DatElemSet.pointersVHE.p_indActSTmacro ;
            m_indActSTmacro = m_VarHistElem(p_indActSTmacro,:)  ; 
            
            p_m_tens        = e_DatElemSet.pointersVAE.p_m_tens ;
            m_tens          =  m_VarAuxElem(p_m_tens    ,:)  ;
            m_tens          =  reshape(m_tens,4,2,[]);
            %indi_elem       = find (m_indActSTmacro == 2);
            
            % La normal la postprocesara desde el instante en el que ha bifurcado
            indi_elem       = find (m_indActSTmacro==1|m_indActSTmacro==2);
            
            mx           = m_tens(1,1,:);
            my           = m_tens(2,2,:);
            
            for iElem = 1:nElem
                if find(iElem==indi_elem)
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),mx(iElem),my(iElem),0]');
                else
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),0,0,0]');
                end
            end
            fprintf(fid_res,'End Values  \n');

            
            %*****************************************************************
            % VECTOR DE REFERENCIA PARA LA SELECION DE LA NORMAL DE INYECCION
            %*****************************************************************
            fprintf(fid_res,['Result "Finite Element//SDA//Ref. Vector for Injection" "Load Analysis" %d ',...
                'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            if ndime==2
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP"\n');
            elseif ndime==3
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP" "Z-COMP"\n');
            end
            fprintf(fid_res,'Values\n');
               
            p_indActSTmacro = e_DatElemSet.pointersVHE.p_indActSTmacro ;
            m_indActSTmacro = m_VarHistElem(p_indActSTmacro,:)  ; 
            
            p_ref_vector    = e_DatElemSet.pointersVAE.p_ref_vector ;
            ref_vector      = m_VarAuxElem(p_ref_vector    ,:)  ;
            ref_vector      = reshape(ref_vector,4,2,[]);
            %indi_elem       = find (m_indActSTmacro == 2);
            indi_elem     = find (m_indActSTmacro==1|m_indActSTmacro==2);
            
            refVec_x        = ref_vector(1,1,:);
            refVec_y        = ref_vector(2,2,:);
            
            for iElem = 1:nElem
                if find(iElem==indi_elem)
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),refVec_x(iElem),refVec_y(iElem),0]');
                else
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),0,0,0]');
                end
            end
            fprintf(fid_res,'End Values  \n');
            
            
            %***************************************
            % VECTOR DE SUAVIZADO - DOUBLE SMOOTHING
            %***************************************
            fprintf(fid_res,['Result "Finite Element//SDA//Smoothing Vector" "Load Analysis" %d ',...
                'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            if ndime==2
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP"\n');
            elseif ndime==3
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP" "Z-COMP"\n');
            end
            fprintf(fid_res,'Values\n');
               
            %p_indActSTmacro = e_DatElemSet.pointersVHE.p_indActSTmacro ;
            %m_indActSTmacro = m_VarHistElem(p_indActSTmacro,:)  ; 
            
            p_nSmoothing  = e_DatElemSet.pointersVAE.p_nSmoothing ;
            %p_ref_vector    = e_DatElemSet.pointersVAE.p_ref_vector ;
            smooth_vector      = m_VarAuxElem(p_nSmoothing , :);
            smooth_vector      = reshape(smooth_vector,4,2,[]);
            %indi_elem       = find (m_indActSTmacro == 2);
            
            smooth_vector_x        = smooth_vector(1,1,:);
            smooth_vector_y        = smooth_vector(2,2,:);
            
            for iElem = 1:nElem
                %if find(iElem==indi_elem)
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),smooth_vector_x(iElem),smooth_vector_y(iElem),0]');
                %else
                %    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),0,0,0]');
                %end
            end
            fprintf(fid_res,'End Values  \n');            
            
            %*****************************************************************
            % VECTOR DE REFERENCIA PARA EL DOBLE SUAVIZADO
            %*****************************************************************
            fprintf(fid_res,['Result "Finite Element//SDA//Ref. Vector for Smooth." "Load Analysis" %d ',...
                'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            if ndime==2
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP"\n');
            elseif ndime==3
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP" "Z-COMP"\n');
            end
            fprintf(fid_res,'Values\n');
               
            p_nRefSmoothing = e_DatElemSet.pointersVAE.p_nRefSmoothing ;
            m_nRefSmoothing = m_VarAuxElem(p_nRefSmoothing,:)  ; 
            
            nRefSmoothing_x = m_nRefSmoothing(1,:);
            nRefSmoothing_y = m_nRefSmoothing(2,:);
            
            for iElem = 1:nElem
                %if find(iElem==indi_elem)
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),nRefSmoothing_x(iElem),nRefSmoothing_y(iElem),0]');
                %else
                %    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),0,0,0]');
                %end
            end
            fprintf(fid_res,'End Values  \n');
               
            %*****************************************************************
            % VECTOR DE REFERENCIA PARA EL DOBLE SUAVIZADO
            %*****************************************************************
            fprintf(fid_res,['Result "Finite Element//SDA//GradNormU" "Load Analysis" %d ',...
                'Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
            if ndime==2
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP"\n');
            elseif ndime==3
                fprintf(fid_res,'ComponentNames "X-COMP" "Y-COMP" "Z-COMP"\n');
            end
            fprintf(fid_res,'Values\n');
               
            p_GradNormU = e_DatElemSet.pointersVAE.p_GradNormU ;
            m_GradNormU = m_VarAuxElem(p_GradNormU,:)  ; 
            
            GradNormU_x = m_GradNormU(1,:);
            GradNormU_y = m_GradNormU(2,:);
            
            for iElem = 1:nElem
                %if find(iElem==indi_elem)
                    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),GradNormU_x(iElem),GradNormU_y(iElem),0]');
                %else
                %    fprintf(fid_res,'%7i  %15.5e  %15.5e  %15.5e \n ',[m_NumElem(iElem),0,0,0]');
                %end
            end
            fprintf(fid_res,'End Values  \n');
            
            %********************************************
            % Tensiones
            %********************************************
            switch struhyp
                case {1,2} %Deformaci�n plana y Tensi�n plana
                    tipoDatAlmac = 'PlainDeformationMatrix';
                    nomComponente = '"Stress XX" "Stress YY" "Stress XY" "Stress ZZ"';
                case 3 %Tridimensional
                    tipoDatAlmac = 'Matrix';
                    nomComponente = 'Stress XX" "Stress YY" "Stress ZZ" "Stress XY" "Stress XZ" "Stress YZ';
            end
            
            % TENSIONES IMPLEX
            % ****************
            STRESS= reshape(stress,ntens,npg,[]);
            STRESS_STD_IMPLEX=STRESS(:,1:4,:);
            STRESS_STD_IMPLEX = reshape(STRESS_STD_IMPLEX,ntens*NPG_ELE,[]);
            
            fprintf(fid_res,['Result "Stresses (Implex)//On Gauss Points" "Load Analysis" %d %s ',...
                'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
            fprintf(fid_res,'ComponentNames %s\n',nomComponente);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,NPG_ELE)];
            fprintf(fid_res,format,[m_NumElem;STRESS_STD_IMPLEX]);
            fprintf(fid_res,'End Values\n');
            
            
             % TENSION PUNTO DE GAUSS SINGULAR
             % *******************************            
             STRESS_SING=STRESS(:,6,:); 
             fprintf(fid_res,['Result "Stresses (Tilde) Sing. GP. //On Gauss Points" "Load Analysis" %d %s ',...
                 'OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep,tipoDatAlmac);
             fprintf(fid_res,'ComponentNames %s\n',nomComponente);
             fprintf(fid_res,'Values\n');
             format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,1)];
             fprintf(fid_res,format,[m_NumElem;STRESS_SING]);
             fprintf(fid_res,'End Values\n');             
            
            
            % TENSIONES TILDE (AUTOEQUILIBRADAS)
            % **********************************
            i_stressTilde  =  e_DatElemSet.pointersVHE.i_stressTilde ;
            m_stressTilde  =   m_VarHistElem(i_stressTilde   ,:)  ; 
            
            STRESS= reshape(m_stressTilde,ntens,npg,[]);
            STRESS_STD=STRESS(:,1:4,:);
            STRESS_STD= reshape(STRESS_STD,ntens*NPG_ELE,[]);
            
            fprintf(fid_res,['Result "Stresses (Tilde) //On Gauss Points" "Load Analysis" %d %s ',...
                'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
            fprintf(fid_res,'ComponentNames %s\n',nomComponente);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,NPG_ELE)];
            fprintf(fid_res,format,[m_NumElem;STRESS_STD]);
            fprintf(fid_res,'End Values\n');   
             
            % ************************************************
            % STRAINS
            % ************************************************
            switch struhyp
                case {1,2} %Deformaci�n plana y Tensi�n plana
                    tipoDatAlmac = 'PlainDeformationMatrix';
                    %nomComponente = '"Strain XX" "Strain YY" "Strain XY" "Strain ZZ"';
                    nomComponente = '"Strain XX" "Strain YY" "Strain ZZ" "Strain XY"';
                case 3 %Tridimensional
                    tipoDatAlmac = 'Matrix';
                    nomComponente = 'Strain XX" "Strain YY" "Strain ZZ" "Strain XY" "Strain XZ" "Strain YZ';
            end
            fprintf(fid_res,['Result "Strains//On Gauss Points" "Load Analysis" %d %s ',...
                'OnGaussPoints "',nomGaussSet,'"\n'],iStep,tipoDatAlmac);
            fprintf(fid_res,'ComponentNames %s\n',nomComponente);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat([repmat(' %.15g',1,ntens),'\n'],1,npg)];
            fprintf(fid_res,format,[m_NumElem;strain]);
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

        otherwise% Triangulo de 3 nodos con discontinuidades fuertes (SDA)
            switch struhyp
                case {1,2} %Deformaci�n plana y Tensi�n plana
                    tipoDatAlmac = 'PlainDeformationMatrix';
                    nomComponente = '"Stress XX" "Stress YY" "Stress XY" "Stress ZZ"';
                case 3 %Tridimensional
                    tipoDatAlmac = 'Matrix';
                    nomComponente = 'Stress XX" "Stress YY" "Stress ZZ" "Stress XY" "Stress XZ" "Stress YZ';
            end
            fprintf(fid_res,['Result "Stresses//On Gauss Points" "Load Analysis" %d %s ',...
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
                    nomComponente = '"Strain XX" "Strain YY" "Strain XY" "Strain ZZ"';
                case 3
                    %Tridimensional
                    tipoDatAlmac = 'Matrix';
                    nomComponente = ['Strain XX" "Strain YY" "Strain ZZ" "Strain XY" ',...
                        '"Strain XZ" "Strain YZ'];
                    %case 4
                    %Axisimetr�a
                    %case 5
                    %Barras 2D
                otherwise
                    error(['GidPost: Fluctuaciones de Deformaciones en los puntos de Gauss: Componentes: ',...
                        'No definidas para esta hip�tesis de carga']);
            end
            fprintf(fid_res,['Result "Strain//Fluctuations//On Gauss Points" ',...
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
            fprintf(fid_res,['Result "Strain//Total//On Gauss Points" ',...
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
            fprintf(fid_res,['Result "Strain//Macro applied//Over Element" ',...
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
            fprintf(fid_res,['Result "Strain//Localized Domain//',...
                'Sobre elemento usado" "Load Analysis" %d Scalar OnGaussPoints "',...
                nomGaussUnicoSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            format = '%d %d\n';
            fprintf(fid_res,format,[m_NumElem;ELOC(1,m_IndElemSet)]);
            fprintf(fid_res,'End Values\n');
            
            if size(ELOC,1)>1
                % Condici�n de incremento de deformaci�n total que se ha verificado
                fprintf(fid_res,['Result "Strain//Localized Domain//',...
                    'Over element depending the strain increment" "Load Analysis" %d ',...
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
            
            if MOD_TYPE~=2
                % Plasticidad J2
                % Deformacion plastica equivalente sobre puntos de Gauss
                fprintf(fid_res,['Result  "Constitutive Model//Plastic Equiv. Strain//',...
                    'Sobre punto de Gauss" "Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
                fprintf(fid_res,'Values\n');
                format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
                fprintf(fid_res,format,[m_NumElem;hvar_new(ntens+1:sihvarpg:end,:)]);
                fprintf(fid_res,'End Values\n');
                
                % Deviatoric part S
                fprintf(fid_res,['Result  "Constitutive Model//Norm of deviatoric Stress//',...
                    'Sobre punto de Gauss" "Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
                fprintf(fid_res,'Values\n');
                format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
                fprintf(fid_res,format,[m_NumElem;hvar_new(ntens*2+3:sihvarpg:end,:)]);
                fprintf(fid_res,'End Values\n');
                
                % �ndice de carga fload
                fprintf(fid_res,['Result  "Constitutive Model//Load Index//On Gauss Points" ',...
                    '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
                fprintf(fid_res,'Values\n');
                format = ['%d',repmat(' %d\n',1,NPG_ELE)];
                fload = hvar_new(ntens+2:sihvarpg:end,:);
                fprintf(fid_res,format,[m_NumElem;fload]);
                fprintf(fid_res,'End Values\n');
                
                % �ndice de carga fload (por elemento)
                fprintf(fid_res,['Result  "Constitutive Model//Load Index//',...
                    'Sobre elemento (any PG)" "Load Analysis" %d Scalar OnGaussPoints "',...
                    nomGaussUnicoSet,'"\n'],iStep);
                fprintf(fid_res,'Values\n');
                format = '%d %d\n';
                %Se adopta si un elemento tiene un PG con fload!=0, el elemento adopta un fload 1.
                fprintf(fid_res,format,[m_NumElem;any(fload)]);
                fprintf(fid_res,'End Values\n');
            end
            
         case {10,11,12,13}  % elemento banda de la microcelda modelo Barcelona
            % Variable de da�o
            
            if MOD_TYPE~=2
                
                fprintf(fid_res,['Result "Constitutive Model//Damage//On (1st-4th) Gauss Points" ',...
                    '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
                fprintf(fid_res,'Values\n');
                format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
                %fprintf(fid_res,format,[m_NumElem;hvar_new(1:sihvarpg:end,:)]);
                fprintf(fid_res,format,[m_NumElem;aux_var(1:siavarpg:end,:)]);
                fprintf(fid_res,'End Values\n');
                
                if (eltype == 21)
                    nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
                    %nomGaussSet     = ['GP_Set_',num2str(iSet)];
                    m_NumElem       = e_DatSet(iSet).m_NumElem;
                    %p_varHistSmooth = sihvarpg*5+1;
                    p_varHistSmooth = siavarpg*5+1;
                    %hvar = e_VarEst(iSet).hvar(p_varHistSmooth,:);
                    %avar = e_VarAux(iSet).VarAuxGP;
                    avar = e_VarAux(iSet).VarAuxGP(p_varHistSmooth,:);
                    %dvar = avar(p_varHistSmooth,:);
                    %
                    fprintf(fid_res,['Result "Finite Element//SDA//Damage PG(6)" "Load Analysis" %d ',...
                        'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
                    fprintf(fid_res,'Values\n');
                    %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
                    %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
                    format = ['%d',repmat(' %.15g',1,1),'\n'];
                    %Se imprime la variable del elemento en un punto de gauss central.
                    %fprintf(fid_res,format,[m_NumElem;hvar]);
                    fprintf(fid_res,format,[m_NumElem;avar]);
                    fprintf(fid_res,'End Values \n');
                end
                
                % -------------------------------------
                
                % Incremento de variable interna (strain-like) del modelo de danio isotropo
                fprintf(fid_res,['Result "Constitutive Model//implicit delta_r//On (1st-4th) Gauss Points" ',...
                    '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
                fprintf(fid_res,'Values\n');
                format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
                fprintf(fid_res,format,[m_NumElem;hvar_new(5:sihvarpg:end,:)]);
                %fprintf(fid_res,format,[m_NumElem;hvar_new(3:sihvarpg:end,:)]);
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
                %    %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
                %    %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
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
                fprintf(fid_res,format,[m_NumElem;hvar_new(1:sihvarpg:end,:)]);
                %fprintf(fid_res,format,[m_NumElem;hvar_new(5:sihvarpg:end,:)]);
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
                %    %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
                %    %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
                %    format = ['%d',repmat(' %.15g',1,1),'\n'];
                %    %Se imprime la variable del elemento en un punto de gauss central.
                %    fprintf(fid_res,format,[m_NumElem;hvar]);
                %    fprintf(fid_res,'End Values \n');
                %end
                
                % -----------------------------------------------
                
            end
            
        %case {4,5,10}   % antes modelo 11
        case {4,5}   % antes modelo 11
            % Da�o isotropo est�ndar y regularizado
            % Variable de da�o
            fprintf(fid_res,['Result "Constitutive Model//Damage//On Gauss Points" ',...
                '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat(' %.15g\n',1,NPG_ELE)];
            %fprintf(fid_res,format,[m_NumElem;hvar_new(3:sihvarpg:end,:)]);
            fprintf(fid_res,format,[m_NumElem;aux_var(1:siavarpg:end,:)]);
            fprintf(fid_res,'End Values\n');
            %*************************************************************************************
            % �ndice de carga FLoad
            %*************************************************************************************
            fprintf(fid_res,['Result "Constitutive Model//Load Index//On Gauss Points" ',...
                '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussSet,'"\n'],iStep);
            %Para utilizar la tabla de resultados previamente creada.
            %fprintf(fid_res,['ResultRangesTable ','"Factor de Carga de Danio Isotropico"\n']);
            fprintf(fid_res,'Values\n');
            format = ['%d',repmat(' %d\n',1,NPG_ELE)];
            %fprintf(fid_res,format,[m_NumElem;hvar_new(4:sihvarpg:end,:)]);
            fprintf(fid_res,format,[m_NumElem;aux_var(2:siavarpg:end,:)]);
            fprintf(fid_res,'End Values\n');
            
            %Impresi�n de los valores centrales para el elemento Q1MixtoInjDef en forma separada
            if eltype==20
                % Variable de da�o
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
            %Este modelo (conshyp==8) est� realizado para elementos triangulares lineales (donde se
            %utiliza un solo de punto de Gauss). Por ello en forma indiferente a la cantidad de puntos de
            %Gauss, se imprime el salto, normal, tangente en un solo punto de Gauss.
            
            %****************************************************************
            % Salto para el elemento de Fuerzas Centrales con fisuras de da�o
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
            fprintf(fid_res,['Result "Constitutive Model//Jump (Beta)//Normal vector to the crack" ',...
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
            fprintf(fid_res,['Result "Constitutive Model//Jump (Beta)//Tangent vector to the crack" ',...
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
            fprintf(fid_res,['Result "Constitutive Model//Load Index" "Load Analysis" %d Scalar ',...
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
            fprintf(fid_res,['Result "Constitutive Model//Adopted ultimate stress" "Load Analysis" ',...
                '%d Scalar OnGaussPoints "',nomGaussUnicoSet,'\n'],iStep);
            %Valores (al usarse un solo punto de gauss, se aplica un valor por elemento)
            fprintf(fid_res,'Values\n');
            fprintf(fid_res,'%d %d\n',[m_NumElem;hvar_new(4,:)]);
            fprintf(fid_res,'End Values\n');
            
        case 50
            
            %Se recupera las deformaciones de los dos PG, que fue cortada previamente para imprimir un
            %solo PG a nivel macro.
            f_PostProcME(iStep,iStage,m_NumElem,reshape(strain,ntens,npg,[]),e_DatMatSet,hvar_new,hvar_old,e_VG)
            
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
            
            if (eltype == 21)
                nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];                
                i_vectVHElem   =  e_DatSet(iSet).e_DatElem.pointersVHE.i_vectVHElem;                   
                m_NumElem       = e_DatSet(iSet).m_NumElem;
                p_varHistDamage= i_vectVHElem(8);
                VarHistElem = e_VarEst(iSet).VarHistElem(p_varHistDamage,:);
                %
                fprintf(fid_res,['Result "Finite Element//SDA//Damage PG(6)" "Load Analysis" %d ',...
                    'Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
                fprintf(fid_res,'Values\n');
                %Se considera ndime componentes, aunque no est� implementado este modelo constitutivo para el
                %Se considera ndime componentes, aunque no est� implementado para el caso 3D.
                format = ['%d',repmat(' %.15g',1,1),'\n'];
                %Se imprime la variable del elemento en un punto de gauss central.
                fprintf(fid_res,format,[m_NumElem;VarHistElem]);
                fprintf(fid_res,'End Values \n');
            end  
                        
            %          cmdWinDoc = com.mathworks.mde.cmdwin.CmdWinDocument.getInstance;
            %          listeners = cmdWinDoc.getDocumentListeners;
            %          jFxCommandArea = listeners(3);
            %          set(jFxCommandArea,'Background','red');
            %Se recupera las deformaciones de los dos PG, que fue cortada previamente para imprimir un
            %solo PG a nivel macro.
            f_PostProcMEBcna(iStep,iStage,m_NumElem,reshape(strain,ntens,npg,[]),e_DatMatSet,hvar_new,hvar_old,...
                aux_var,e_VG)
            %          set(jFxCommandArea,'Background','yellow');
        
           
            
        otherwise
%             error('PostProceso: Modelo constitutivo no definido.')
    end
end

fclose(fid_res);
