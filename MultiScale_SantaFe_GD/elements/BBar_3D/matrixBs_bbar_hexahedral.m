function [Bs,detJ] = matrixBs_bbar_hexahedral(coord_n,xg,e_VG)

%global ndn dofpe ntens
ndn = e_VG.ndn;
dofpe = e_VG.dofpe;
ntens = e_VG.ntens;

E = xg(1);
n = xg(2);
z = xg(3);

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

% Matriz B standard:
% ииииииииииииииииии
Bs = zeros(ntens,dofpe);
dN_xyz = J\dN_Enz;
Bs(1,1:ndn:dofpe) = dN_xyz(1,:);
Bs(2,2:ndn:dofpe) = dN_xyz(2,:);
Bs(3,3:ndn:dofpe) = dN_xyz(3,:);
Bs(4,1:ndn:dofpe) = dN_xyz(2,:);
Bs(4,2:ndn:dofpe) = dN_xyz(1,:);
Bs(5,1:ndn:dofpe) = dN_xyz(3,:);
Bs(5,3:ndn:dofpe) = dN_xyz(1,:);
Bs(6,2:ndn:dofpe) = dN_xyz(3,:);
Bs(6,3:ndn:dofpe) = dN_xyz(2,:);