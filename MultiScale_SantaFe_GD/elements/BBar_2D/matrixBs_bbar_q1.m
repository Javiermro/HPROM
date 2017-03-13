function [Bs,detJ] = matrixBs_bbar_q1(coord_n,xg,e_DatElemSet,e_VG)

ndn = e_VG.ndn;
ntens = e_VG.ntens;
dofpe = e_DatElemSet.dofpe;

E = xg(1);
n = xg(2);

% Derivadas de Ni respecto "E" y "n":
% ииииииииииииииииииииииииииииииииии
dN_E = [(1/4*n-1/4)  (-1/4*n+1/4)  (1/4*n+1/4)  (-1/4*n-1/4)];
dN_n = [(1/4*E-1/4)  (-1/4*E-1/4)  (1/4*E+1/4)  (-1/4*E+1/4)];
dN_En = [dN_E ; dN_n];

% Matriz Jacobiana:
% ииииииииииииииии
J11 = coord_n(1,:)*dN_E';
J12 = coord_n(2,:)*dN_E';
J21 = coord_n(1,:)*dN_n';
J22 = coord_n(2,:)*dN_n';
J = [J11 J12 ; J21 J22];
detJ = det(J);

% Matriz B:
% ииииииии
Bs = zeros(ntens,dofpe);
dN_xy = J\dN_En;
Bs(1,1:ndn:dofpe) = dN_xy(1,:);
Bs(2,2:ndn:dofpe) = dN_xy(2,:);
Bs(4,1:ndn:dofpe) = dN_xy(2,:);
Bs(4,2:ndn:dofpe) = dN_xy(1,:);

