function plot_contour (xx,conec,in,eltype,tline)

inn = zeros(max(in),1);
inn(in(:)) = 1:length(in);

[nElem,npe] = size(conec);

if (eltype == 5);
   npe = npe-1;
end

X = zeros(npe+1,nElem);
Y = zeros(npe+1,nElem);
x = zeros(size(xx,1),1);
y = zeros(size(xx,1),1);

for k = 1:npe,
   X(k,:) = xx(inn(conec(:,k)),1)';
   Y(k,:) = xx(inn(conec(:,k)),2)';
end
X(npe+1,:) = X(1,:);
Y(npe+1,:) = Y(1,:);

% Grafica elementos
plot(X,Y,tline);

% Imprime número de nodos:
nnodstr = num2str(in);
for i = 1:size(in,1);
   x(i) = xx(inn(i),1);
   y(i) = xx(inn(i),2);
end
hold on
plot(x,y,'ko','linewidth',1,'markerfacecolor','w','markersize',4);
text(x,y,nnodstr,'fontsize',7,'verticalalignment','bottom','horizontalalignment','left','color','k');

% Configuración de la figura:
set(gcf,'color',[0.75 0.75 0.75]);
set(gca,'color',[0.75 0.75 0.75],'fontsize',8);
axis equal ; hold off ; zoom on ; axis off;
title('Malla de elementos finitos');
