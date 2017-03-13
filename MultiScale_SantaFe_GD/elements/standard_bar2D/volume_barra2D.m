function [vol] = volume_barra_2D(coord_n,Area)

%******************************************************************************************
%*  EVALUACION DEL VOLUMEN ARA UNA BARRA 2D                                               *
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

% Longitud del elemento
% иииииииииииииииииииии
dx  = coord_n(1,2) - coord_n(1,1);
dy  = coord_n(2,2) - coord_n(2,1);
L   = sqrt(dx^2+dy^2);
vol = L*Area;