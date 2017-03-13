% *********************************************
% *  POST-PROCESO BASICO: CURVAS DE RESPUESTA *
% * *******************************************

%Sobreescribir gráficas existentes
esSobreEs = 0;
%
struhyp = e_VG.struhyp;
dofrs = e_VG.dofrs;
%A parte de Read_data, pareciera que es el único lugar donde se usa DOFRS.
DOFRS = e_VG.DOFRS;
fileCompleto = e_VG.fileCompleto;
%Dtime = e_VG.Dtime; %Está comentado el código donde se usa Dtime
max_time = e_VG.max_time;

% CURVA DE RESPUESTA CARGA-DESPLAZAMIENTO
textoEjes = {'X','Y','Z'};
for iFig = 1:size(dofrs,1)
   uPlot = u_dofrs(:,iFig)*DOFRS(iFig,4);
   fPlot = f_dofrs(:,iFig)*DOFRS(iFig,5);
   nodo = DOFRS(iFig,1);
   %% Se guarda en archivo cur los resultados
   fid = fopen([fileCompleto,'.cur',num2str(iFig,'%03d')],'wt');
   fprintf(fid,'%-.8f %12.8f\n',[uPlot,fPlot]');
   fclose(fid);
   %% Se plotea los resultados
   %Sobreescribe las gráficas, usando el handle igual al número de figuras (puede sobreescribir
   %sobre figuras que no tienen nada que ver con el programa)
   if esSobreEs==1
      figure(iFig);
   else
      figure;
   end
   hold on
   plot(uPlot,fPlot,'-ok','markersize',2);
   grid on
   title('Grafica Fuerza vs Desplazamiento');
   xtit = ['Desplazamiento ',textoEjes{DOFRS(iFig,2)},', Nodo Nº ',num2str(nodo)];
   ytit = ['Fuerza ',textoEjes{DOFRS(iFig,3)},', Nodo Nº ',num2str(nodo)];
   xlabel(xtit)
   ylabel(ytit)
   hold off
end

% CURVA DE RESPUESTA CARGA-TIEMPO
% for i = 1:length(dofrs)
%     figure
%     nodo = DOFRS(i,1);
%     if (DOFRS(i,2) == 1)
%         dof_direction = 'X';
%     else
%         dof_direction = 'Y';
%     end
%     plot([0:Dtime:max_time],f_dofrs(:,i),'-ok','markersize',2);
%     grid on;
%     title('Grafica Fuerza vs Tiempo');
%     xtit = strcat('Tiempo [seg]');
%     ytit = strcat('Fuerza ',dof_direction,', Nodo Nº ',num2str(nodo));
%     xlabel(xtit);
%     ylabel(ytit);
% end

% DEFORMADA PARA EL CASO DE ELEMENTOS DE BARRA
if struhyp==5
   plot_reticulado(conec,xx,u,hvar_new,1);
end