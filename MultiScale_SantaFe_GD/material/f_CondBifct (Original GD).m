function [condBif,m_thetaBif] = f_CondBifct(m_ct,varargin)

   if nargin==2
      %f_CondBifct(m_ct,e_VG)
      impr = 0;
      e_VG = varargin{1};
   elseif nargin==3
      %f_CondBifct(m_ct,nomArch,e_VG)      
      impr = 1;
      nomArch = varargin{1};
      e_VG = varargin{2};      
   else
      error('Análisis de bifurcación: Número de argumentos de llamada no definido.')
   end
   
   delTheta = 0.0025*pi;
   %delTheta = 0.00001*pi;
   facLimBusq = 0.10;

   theta = (-pi/2):delTheta:(pi/2-delTheta/2);
   nTheta = size(theta,2);
   detQ = zeros(1,nTheta);

   %Versión semi-vectorizada.
   m_nVoigt = f_Vec2Voigt2DVect([cos(theta)',sin(theta)'],e_VG);
   for ind = 1:nTheta
      %Tensor de localizacion
      detQ(ind) = det(m_nVoigt(:,:,ind)'*m_ct*m_nVoigt(:,:,ind));
   end

   if impr
      %nomArch = [nomArch,'.analisisBif'];
      fId = fopen(nomArch,'at');
      fprintf(fId,'#Paso: %d\n',e_VG.istep);
      fprintf(fId,'#Ángulo Det(Q)\n');
      fprintf(fId,'%f %f\n',[theta;detQ]);
      fprintf(fId,'\n\n');
      fclose(fId);
   end

   [mindetQ1,m_IndBif1] = min(detQ);
   %Se utiliza un detQ auxiliar para determinar los próximos mínimos, donde los mínimos encontrados
   %se los reeamplaza con infinito, y así ignorarlos. Y el detQ original se utiliza para calcular
   %las pendientes (si no se hace esto, al reemplazar con inf, se puede detectar mínimos en forma
   %errónea).
   detQAux = detQ;
   detQAux(m_IndBif1) = inf;
   %Se hace la búsqueda del segundo mínimo. Itera hasta encontrarlo (podría llegar a encontrar un mínimo
   %local, pero se descarta con facLimBusq). Notar que ignora los máximos locales por el control de las
   %derivadas.
   %Se limita la búsqueda del segundo límite hasta quede en un porcentaje del codominio total de detQ por
   %debajo del mínimo absoluto.
   sdoMin = 0;
   iiTer = 1;
   %Para asegurar que entre al loop se impone un valor inicial que asegure que la resta sea menor al 
   %codBusq2Min
   mindetQ2 = -inf;
   codBusq2Min = facLimBusq*(max(detQ)-mindetQ1);
   while (mindetQ2-mindetQ1)<codBusq2Min&&iiTer<nTheta
   %for i = 1:nTheta-1
      [mindetQ2,m_IndBif2] = min(detQAux);
      %Para considerar el anterior y el siguiente, del primero y el último de la matriz, respectivamente,
      %se ingresa condiciones adicionales.
      if (m_IndBif2==1&&detQ(1)-detQ(nTheta)<0&&detQ(2)-detQ(1)>0)||...
            (m_IndBif2==nTheta&&detQ(nTheta)-detQ(nTheta-1)<0&&detQ(1)-detQ(nTheta)>0)||...
            (m_IndBif2>1&&m_IndBif2<nTheta&&detQ(m_IndBif2)-detQ(m_IndBif2-1)<0&&...
            detQ(m_IndBif2+1)-detQ(m_IndBif2)>0) 
         sdoMin = 1;
         break
      end
      detQAux(m_IndBif2) = inf;
      iiTer = iiTer+1;
   end
   if ~sdoMin
      m_IndBif2 = m_IndBif1;
      mindetQ2 = mindetQ1;
   end
   m_thetaBif = theta([m_IndBif1,m_IndBif2]);
   m_mindetQ = [mindetQ1,mindetQ2];
   %
   if mindetQ1<0
      condBif = 1;
      %Imprensión en pantalla
      %fprintf('Análisis de bifurcación:\n')
      %fprintf('theta=%f Det(Q)=%f\n',[m_thetaBif;m_mindetQ])
      %Impresión en archivo de paso que bifurca
      if impr
         %nomArch = [e_VG.fileCompleto,'.analisisBif'];
         fId = fopen(nomArch,'at');
         %Se usa # porque es el comentario del GNUPlot.
         fprintf(fId,'#################################################\n');
         fprintf(fId,'#Bifurcación, paso: %d (t=%f)\n',e_VG.istep,e_VG.Dtime*e_VG.istep);
         fprintf(fId,'#theta=%f Det(Q)=%f\n',[m_thetaBif;m_mindetQ]);
         fprintf(fId,'#################################################\n');
         fclose(fId);
      end
   else
      condBif = 0;
%       m_thetaBif = zeros(1,2);
   end

end