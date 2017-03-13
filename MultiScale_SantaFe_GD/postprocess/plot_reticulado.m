function plot_reticulado (elem,node,u,hvar_new,scale)

figure
nelem = size(elem,1);
inn = 1:size(node,1);
hold off;
for i=1:nelem
   nodo1=elem(i,1);
   nodo2=elem(i,2);
   gdl_global = [(nodo1-1)*2+1, (nodo1-1)*2 + 2, (nodo2-1)*2 + 1, (nodo2-1)*2 + 2]; 
   ue_global = u(gdl_global);
   plot([node(nodo1,1),node(nodo2,1)],[node(nodo1,2),node(nodo2,2)],'-ob','markerfacecolor',[0.75 0.75 0.75]);
   hold on;
   despx = scale*[ue_global(1) ue_global(3)];
   despy = scale*[ue_global(2) ue_global(4)];
   plot([node(nodo1,1),node(nodo2,1)]+despx,[node(nodo1,2),node(nodo2,2)]+despy,'-or','markerfacecolor',[0.75 0.75 0.75]);
end
axis equal;
axis off;
hold on
for i = 1:length(inn)
   text(node(i,1),node(i,2),num2str(inn(i)),'verticalAlignment','bottom');
end

barra_plast = find(hvar_new(2:2:length(hvar_new)) ~= 0);
for i = 1:length(barra_plast)
   nbarra = barra_plast(i);
   ni = elem(nbarra,1);
   nf = elem(nbarra,2);
   hold on
   plot([node(ni,1)+u(2*ni-1) node(nf,1)+u(2*nf-1)],[node(ni,2)+u(2*ni) node(nf,2)+u(2*nf)],...
       '-or','linewidth',3,'markerfacecolor',[0.75 0.75 0.75]);
end
title('Configuracion deformada final');