function [vol] = volume_bbar_hexahedral(coord_n,e_VG)

%global npg wg xg
npg = e_VG.npg;
wg = e_VG.wg;
xg = e_VG.xg;

vol = 0;

for i = 1:npg;
   for j = 1:npg;
      for k = 1:npg;
 
          E = xg(i);
          n = xg(j);
          z = xg(k);

          % Derivadas de Ni respecto "E" y "n":
          % ииииииииииииииииииииииииииииииииии
          dN_E = [ -1/8*(-1+n)*(-1+z),  1/8*(-1+n)*(-1+z),  -1/8*(1+n)*(-1+z),   1/8*(1+n)*(-1+z),...
                    1/8*(-1+n)*(1+z),  -1/8*(-1+n)*(1+z),    1/8*(1+n)*(1+z),   -1/8*(1+n)*(1+z)];
 
          dN_n = [ -1/8*(-1+E)*(-1+z),   1/8*(1+E)*(-1+z),  -1/8*(1+E)*(-1+z),  1/8*(-1+E)*(-1+z),...
                    1/8*(-1+E)*(1+z),   -1/8*(1+E)*(1+z),    1/8*(1+E)*(1+z),  -1/8*(-1+E)*(1+z)];
 
          dN_z = [ -1/8*(-1+E)*(-1+n),   1/8*(1+E)*(-1+n),   -1/8*(1+E)*(1+n),   1/8*(-1+E)*(1+n),...
                    1/8*(-1+E)*(-1+n),  -1/8*(1+E)*(-1+n),    1/8*(1+E)*(1+n),  -1/8*(-1+E)*(1+n)];
 
          dN_Enz = [dN_E ; dN_n ; dN_z];

          % Matriz Jacobiana:
          % ииииииииииииииии
          J11 = coord_n(1,:)*dN_E';
          J12 = coord_n(2,:)*dN_E';
          J13 = coord_n(3,:)*dN_E';
          J21 = coord_n(1,:)*dN_n';
          J22 = coord_n(2,:)*dN_n';
          J23 = coord_n(3,:)*dN_n';
          J31 = coord_n(1,:)*dN_z';
          J32 = coord_n(2,:)*dN_z';
          J33 = coord_n(3,:)*dN_z';

          J    = [J11 J12 J13; J21 J22 J23; J31 J32 J33];
          detJ = det(J);
          
          vol = vol + wg(i)*wg(j)*wg(k)*detJ;
          
      end
   end
end