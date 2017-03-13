function [le_elem,fii] = le_quad_q1_epd(coord_n,new_coord_n,bif_angle_rad)

%*********************************************************************************
%*  function [B,detJ] = bmatrix_CL4N_EP (coord_n,xg,npe)                         *
%*  Obtiene la matriz B del cuadrแngulo lineal de 4 nodos para EPD ๓ EPT         *
%*  evaluada en xg.                                                              *
%*                                                                               *
%*  ARGUMENTOS DE ENTRADA:                                                       *
%*  จจจจจจจจจจจจจจจจจจจจจ                                                        *
%*  coord_n : Coordenadas de los nodos que componen el elemento.                 *
%*            coord_n = [x1  x2  x3 x4                                           *
%*                       y1  y2  y3 y4]                                          *
%*  numeraci๓n de nodos: 4-------------3                                         *
%*                       |             |                                         *
%*                       |             |                                         *
%*                       |             |                                         *
%*                       1-------------2                                         *
%*  xg      : Vector de posici๓n del punto de Gauss en estudio.                  *
%*            xg = [E n]                                                         *
%*  npe     : N๚mero de nodos por elemento.                                      *
%*                                                                               *
%*  LISTA DE VARIABLES:                                                          *
%*  จจจจจจจจจจจจจจจจจจ                                                           *
%*  E       : coordenada natural                                                 *
%*  n       : coordenada natural                                                 *
%*  dN_E    : Derivadas de las funciones de forma respecto a "E"                 *
%*            dN_E = [dN1_E  dN2_E  dN3_E dN4_E]                                 *
%*  dN_n    : Derivadas de las funciones de forma respecto a "n"                 *
%*            dN_n = [dN1_n  dN2_n  dN3_n dN4_n]                                 *
%*  dN_En   : Derivadas de las funciones de forma respecto a "E" y a "n"         *
%*  J       : matriz jacobiana de la transformaci๓n                              *
%*            J = [dx_E dy_E ; dx_n dy_n]                                        *
%*  dN_xy   : Derivadas de las funciones de forma respecto a "x" y a "y"         *
%*                                                                               *
%*  ARGUMENTOS DE SALIDA:                                                        *
%*  จจจจจจจจจจจจจจจจจจจจ                                                         *
%*  B       : Matriz B para elemento cuadrangular de 4 nodos en estado plano.    *
%*  detJ    : determinante de la matriz jacobiana.                               *
%*********************************************************************************

%global ndn dofpe ntens
%C ndn = e_VG.ndn;
%C dofpe = e_VG.dofpe;
%C ntens = e_VG.ntens;

%C E = xg(1);
%C n = xg(2);
E = 0;
n = 0;
% le_elem = 0;
fii = zeros(size(coord_n,2),1);

% Derivadas de Ni respecto "E" y "n":
% จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
dN_E = [(1/4*n-1/4)  (-1/4*n+1/4)  (1/4*n+1/4)  (-1/4*n-1/4)];
dN_n = [(1/4*E-1/4)  (-1/4*E-1/4)  (1/4*E+1/4)  (-1/4*E+1/4)];
dN_En = [dN_E ; dN_n];

% Matriz Jacobiana:
% จจจจจจจจจจจจจจจจ
J11 = coord_n(1,:)*dN_E';
J12 = coord_n(2,:)*dN_E';
J21 = coord_n(1,:)*dN_n';
J22 = coord_n(2,:)*dN_n';
J = [J11 J12 ; J21 J22];
%C detJ = det(J);

for inode = 1:size(new_coord_n,2)
    if (new_coord_n(1,inode)>0.0)
        fii(inode) = 1.0;
    end
end

% Matriz dN_xy:
% จจจจจจจจ
dN_xy = J\dN_En;
dN_xy(1,:) = dN_xy(1,:)*cos(bif_angle_rad);
dN_xy(2,:) = dN_xy(2,:)*sin(bif_angle_rad);

le_elem = 1/(sum(dN_xy)*fii);
