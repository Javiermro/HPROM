function plot_mapa (xx,conec,inn,eltype,var,color_range,font_size)

in = (1:size(find(inn ~= 0),1))';

[nElem,npe] = size(conec);

if (eltype == 5);
   npe = npe-1;
end

X = zeros(npe,nElem);
Y = zeros(npe,nElem);
Z = zeros(npe,nElem);
x = zeros(size(xx,1),1);
y = zeros(size(xx,1),1);

for k = 1:npe,
   X(k,:) = xx(inn(conec(:,k)),1)';
   Y(k,:) = xx(inn(conec(:,k)),2)';
   Z(k,:) = var(inn(conec(:,k)))';
end

% Grafica elementos
patch(X,Y,Z);
%colormap(gray);
%color = colormap;
%color = ones(size(color,1),3)-color;colormap(color);


% Configuracion
if (nargin >= 6) 
   caxis([0 0.004]);
end
colorbar;
if (nargin >= 7) 
   a = get(gcf,'children');
   set(a(2),'fontsize',14)
end
axis equal;
hold off;
axis off;
grid on;
