function f_PostProcMECohesivo(iStep,m_NumElemMacro,m_DefMacro,e_DatMatSetMacro,hvar_newMacro,...
   e_VarAuxPGMacro,e_VGMacro)   
   
   %Se recupera variables micro
   xx = e_DatMatSetMacro.xx;
   in = e_DatMatSetMacro.in;
   e_VG = e_DatMatSetMacro.e_VG;
   m_ElemPGImpr = e_DatMatSetMacro.m_ElemPGImpr;
   e_DatSet = e_DatMatSetMacro.e_DatSet;
   m_SetElem = e_DatMatSetMacro.m_SetElem;
   %Se guarda el paso tiempo a nivel micro.
   e_VG.istep = e_VGMacro.istep;
   %
   postpro_impre_step = e_VG.postpro_impre_step;
   IRES = e_VG.IRES;
   ndime = e_VG.ndime;
   nElem = e_VG.nElem;
   %Nombre completo del archivo micro
   nomArchCompl = e_VG.fileCompleto;
   %   
   for iPGGraf = 1:size(m_ElemPGImpr,2)
      
      iElemSetMacro = m_ElemPGImpr(1,iPGGraf);
      iPGMacro = m_ElemPGImpr(2,iPGGraf);
      %
      m_iDefMacro = m_DefMacro(:,iPGMacro,iElemSetMacro);
      e_VarEst_new = hvar_newMacro(iPGMacro,iElemSetMacro).e_VarEst;
      e_VarAux = hvar_newMacro(iPGMacro,iElemSetMacro).e_VarAux;
      u = hvar_newMacro(iPGMacro,iElemSetMacro).u;
      c_GdlCond = hvar_newMacro(iPGMacro,iElemSetMacro).c_GdlCond;
      Fint = hvar_newMacro(iPGMacro,iElemSetMacro).Fint;
      c_DefMacroCU = hvar_newMacro(iPGMacro,iElemSetMacro).c_DefMacro;
      m_ElemLoc = hvar_newMacro(iPGMacro,iElemSetMacro).m_ElemLoc;
      %     
      %Se asume que cuando ya está activado el modelo multiescala cohesivo cuando la matriz de m_ElemLoc es no
      %vacía.
      if isempty(m_ElemLoc)
         m_ElemLoc = false(1,nElem);
         %Previo a la activación del modelo multiescala cohesivo solo está definido el primer PG.
         c_ProyIncrEps = e_VarAuxPGMacro(1,iElemSetMacro).c_ProyIncrEps;
      else
         c_ProyIncrEps = e_VarAuxPGMacro(iPGMacro,iElemSetMacro).c_ProyIncrEps;
      end
      %Se guarda los elementos localizados, en la primera fila va los elementos de la zona
      %localizada que se utilizó en el cálculo, mientras que en la segunda fila se coloca los
      %elementos localizados obtenidos por el criterio de incremento de deformación (puede
      %diferenciarse en el caso de que se imponga un grupo de elementos como dominio de
      %localización)
      m_ElemLoc = [m_ElemLoc;e_VarAuxPGMacro(iPGMacro,iElemSetMacro).m_ElemLocCalc]; %#ok<AGROW>
      
      %Se cambia el nombre del archivo de postproceso de GiD de la celda unitaria, para
      %individualizar el elemento y el punto de Gauss que corresponde.
      e_VG.fileCompleto = [nomArchCompl,'_EM',num2str(m_NumElemMacro(iElemSetMacro)),'PGM',num2str(iPGMacro)];

      % IMPRESIÓN DE RESULTADOS
      %Notar que se imprime este paso si a nivel macro este paso se imprime, y a la vez, que si a
      %nivel micro debe imprimirse.
      if IRES   
         if ~mod(iStep,postpro_impre_step)
            % DESPLAZAMIENTO TOTAL DE LA MICRO-CELDA
            uTotal = [m_iDefMacro(1),m_iDefMacro(4)/2;m_iDefMacro(4)/2,m_iDefMacro(2)]*xx(:,1:2)'...
               +reshape(u,ndime,[]);
            % POSTPROCESO GiD
            matlab2gid_res(iStep,in,u,c_GdlCond,e_DatSet,e_VarEst_new,e_VarAux,c_DefMacroCU,...
               uTotal,m_ElemLoc,e_VG)
            % Impresión adicional del modelo constitutivo micro para el problema multiescala.
            if e_VG.isPostResMultFile
               [dir,nomb] = fileparts(e_VG.fileCompleto);
               filename_res = fullfile(dir,'GiDRes',[nomb,'_P',num2str(iStep),'.flavia.res']);
            else
               filename_res = [e_VG.fileCompleto,'.flavia.res'];
            end
            fid_res = fopen(filename_res,'at');
            %Ver si hacer así, o dentro de matlab2gid_res poner condicional para que imprima la variables
            %adicionales que se quiere que se imprima. La primera permitiría sacar la impresión de la
            %deformación fluctuante, macro aplicada, etc., de matlab2gid_res, haciéndolo más simple. El
            %segundo tendría la ventaja de mantener el código más unificado.
            c_NormalesMicro = hvar_newMacro(iPGMacro,iElemSetMacro).c_NormalesMicro;
            for iSet = 1:e_VG.nSet
               %nomGaussSet = ['GP_Set_',num2str(iSet)]; 
               nomGaussUnicoSet = ['GP_Unico_Set_',num2str(iSet)];
               nElemSet = e_DatSet(iSet).nElem;
               %nPG = e_DatSet(iSet).e_DatElem.npg;
               m_NumElem = e_DatSet(iSet).m_NumElem;
               m_NormalesMicro = c_NormalesMicro{iSet,1};
               m_lMicro = c_NormalesMicro{iSet,2};
               %Se guarda a partir que se activa el modelo cohesivo multiescala, ya que previo la matriz
               %m_NormalesMicro está vacía y por lo tanto se escapa la impresión.
               if ~isempty(m_NormalesMicro)
                  %m_NormalesMicro = reshape([m_NormalesMicro(1,1,:,:);m_NormalesMicro(2,2,:,:)],...
                  %   ndime*nPG,nElemSet);
                  m_NormalesMicro = reshape([m_NormalesMicro(1,1,:);m_NormalesMicro(2,2,:)],...
                    ndime,nElemSet);
                  fprintf(fid_res,['Result "Multiescala//Normales Micro seleccionadas" ',...
                     '"Load Analysis" %d Vector OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
                  if ndime==2
                     %En este caso el Gid le pone un nombre por defecto a la componente Z y le asigna un valor 0.
                     s_NomComp = 'ComponentNames "nX" "nY"\n';
                  elseif ndime==3
                     s_NomComp = 'ComponentNames "nX" "nY" "nZ"\n';
                  end
                  fprintf(fid_res,s_NomComp);
                  fprintf(fid_res,'Values\n');
                  %format = ['%d',repmat([repmat(' %.15g',1,ndime),'\n'],1,nPG)];
                  format = ['%d',repmat(' %.15g',1,ndime),'\n'];
                  fprintf(fid_res,format,[m_NumElem;m_NormalesMicro]);
                  fprintf(fid_res,'End Values\n');
               end
               if ~isempty(m_lMicro)
                  %Se acepta la convención de que si es infinito no se ingresó ningún valor válido para
                  %longitud micro.
                  indlMicroVal = ~isinf(m_lMicro);
                  if any(indlMicroVal)
                     fprintf(fid_res,['Result "Multiescala//Longitudes Micro seleccionadas" ',...
                        '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
                     fprintf(fid_res,'Values\n');
                     format = '%d %.15g\n';
                     fprintf(fid_res,format,[m_NumElem(indlMicroVal);m_lMicro(indlMicroVal)]);
                     fprintf(fid_res,'End Values\n');
                  end
               end
               %
               fprintf(fid_res,['Result "Multiescala//Proyección de deformación micro fluctuante" ',...
                  '"Load Analysis" %d Scalar OnGaussPoints "',nomGaussUnicoSet,'"\n'],iStep);
               fprintf(fid_res,'Values\n');
               format = '%d %.15g\n';
               m_ProyIncrEps = c_ProyIncrEps{iSet};
               fprintf(fid_res,format,[m_NumElem;m_ProyIncrEps]);
               fprintf(fid_res,'End Values\n');
            end
            %
            fclose(fid_res);
         end
      end
      
      %Impresión gráficas
      f_SaveDatGraf(u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,e_DatSet,m_SetElem,[],[],e_VG) %JLM
      
   end
         
end