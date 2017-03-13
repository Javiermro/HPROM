function [] = f_OperPosConv_MEClasBif(m_F,m_CT,e_DatSet,e_VG)

nElem = e_DatSet.nElem;
nPG = e_DatSet.e_DatElem.npg;
ntens = e_VG.ntens;
%ndime = e_VG.ndime;
%nomArchCompl = e_VG.fileCompleto;
%
%sihvarpg = e_DatMatSet.sihvarpg;
%siavarpg = e_DatMatSet.siavarpg;
esImplex = e_DatSet.e_DatMat.esImplex;
%m_ElemPGImpr = e_DatMatSet.m_ElemPGImpr;
iStep = e_VG.istep;
m_NumElem = e_DatSet.m_NumElem;
%Variables micro
%e_VGMicro       = e_DatMatSet.e_VG;
%e_DatSetMicro   = e_DatMatSet.e_DatSet;
%m_NumElem       = e_DatSet.m_NumElem;

%Nombre de archivo que se guarda los datos de bifurcación.
nomArchBif = [e_VG.fileCompleto,'.elemBif'];
%Se abre el archivo.
fId = fopen(nomArchBif,'at');

%Se elije un punto de gauss para hacer la bifurcación.
iPG = 1;

%parfor iElem = 1:nElem
for iElem = 1:nElem
   
   %Se toma los gradientes de deformación del elemento
   m_FElem = reshape(m_F(:,iElem),ntens,nPG);   
   
   %condBif =  m_VarAuxElem(p_condBif,iElem);
   condBif = 0;

   e_VG.iElemSet = iElem;
   e_VG.iElemNum = m_NumElem(iElem);

   % Operaciones que se realiza hasta que se detecta la bifurcación.
   %if ~condBif

      % evaluar condicion de bifurcacion en el PG macro    
      if esImplex
         m_CTBif = m_CT(:,:,iPG*2,iElem);
      else
         m_CTBif = m_CT(:,:,iPG,iElem);
      end
      if e_VG.struhyp==20
         [condBif,m_AngBif] = f_CondBifctNoSim(m_CTBif,e_VG);
      else
         [condBif,m_AngBif] = f_CondBifct(m_CTBif,e_VG);
      end
      
      %El vector que se está calculando es el tensor normal en la cofiguración de referencia (N), por lo que
      %se calcula el en la configuración de referencia (n).
      m_NBif = [cos(m_AngBif);sin(m_AngBif)];
      if e_VG.struhyp==20
         %% Grandes deformaciones
         m_F2D = [m_FElem(1,iPG),m_FElem(4,iPG);m_FElem(5,iPG),m_FElem(2,iPG)];
         m_InvF2D = m_F2D\eye(2,2);
         %J = Fzz*det(m_F2D);   
         %m_nBif = zeros(size(m_NBif));
         %m_nBif(1,:) = m_InvF2D(1,1)*m_NBif(1,:)+m_InvF2D(1,2)*m_NBif(2,:);
         %m_nBif(2,:) = m_InvF2D(2,1)*m_NBif(1,:)+m_InvF2D(2,2)*m_NBif(2,:);
         %Al vector hay que normalizarlo, pero para el cálculo del ángulo no es necesario, ya que se usa
         %tangente.
         m_nBif = m_InvF2D'*m_NBif;
         m_AngBifDef = atan(m_nBif(2,:)./m_nBif(1,:));
      else
         %% Pequeñas deformaciones (deformación y tensión plana)
         m_nBif = m_NBif;
         m_AngBifDef = m_AngBif;
      end

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

      %m_VarAuxElem(p_condBif,iElem) = condBif;
      %Se almacena los ángulos de bifurcación.
      %m_VarAuxElem(pAngBif,iElem) = m_AngBif;

   %end
   %
   %Para realizar los cálculos inmediatamente después que se detecta bifurcación.
   if condBif      

      %Se ejecuta entre la bifurcación y la activación de la SD, y mientras no se active la SD.

%          %Como criterio se utiliza los datos del análisis de bifurcación previo a detectar bifurcación
%          %(detQ<0), estos han mostrado ser más cercanos a una normal de bifurcación más coherente con el
%          %problema.
%          %En los problemas monoescala no se vio lo mismo (VERIFICAR!)
%          if ~condBifPrev
%             m_VarAuxElem(pAngBif,iElem) = m_AngBifPrev;
%             %Esta línea siguiente no es necesario ya que en la función Normal_vector_selection se selecciona
%             %una nueva normal con los dos ángulos posibles.
%             %m_VarAuxElem(p_n_tens,iElem) = n_tensPrev(:);
%          end
%          %
%          m_AngBif = m_VarAuxElem(pAngBif,iElem)';
%          m_NormBif = [cos(m_AngBif);sin(m_AngBif)];


         %if ~condBifPrev
            %Impresión en pantalla
            fprintf('** Se detectó bifurcación en el elemento %d.\n',m_NumElem(iElem))
            %Impresión en archivo
            %Texto de los ángulos, dependiendo si hay uno o dos.
            if length(m_AngBif)>1
               s_AngBif = sprintf('Ang(N): %f y %f. Ang(n): %f y %f.',m_AngBif,m_AngBifDef);
            else
               s_AngBif = sprintf('Ang(N): %f. Ang(n): %f.',m_AngBif,m_AngBifDef);
            end
            fprintf(fId,['Fin de Paso=%d (t=%f), Elemento=%d: Se detecta bifurcación. ',s_AngBif,'\n'],...
               iStep,e_VG.Dtime*iStep,m_NumElem(iElem));
         %end

   end
%      
end
%Se cierra archivo de impresión de datos de bifurcación.
fclose(fId);

end

