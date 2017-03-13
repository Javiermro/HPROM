function [B,detJ] = matrixB_barra2D (coord_n,xg,e_VG)

%******************************************************************************************
%*  EVALUACION DE LA MATRIZ B PARA UNA BARRA 2D                                           *
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

%global ntens dofpe ndn
ndn = e_VG.ndn;
dofpe = e_VG.dofpe;
ntens = e_VG.ntens;

% Longitud del elemento
% иииииииииииииииииииии
dx = coord_n(1,2) - coord_n(1,1);
dy = coord_n(2,2) - coord_n(2,1);
L  = sqrt(dx^2+dy^2);

% Matriz B
% ииииииии
B    = zeros(ntens,dofpe);
dN_x = [-1/L 1/L];
B(ntens,1:ndn:dofpe) = dN_x;

% Determinante del Jacobiano de la transformacion
% иииииииииииииииииииииииииииииииииииииииииииииии
detJ = L;