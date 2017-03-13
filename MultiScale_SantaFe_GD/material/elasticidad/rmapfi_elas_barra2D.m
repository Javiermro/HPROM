function [sigma_new,hvar_new,aux_var] = rmapfi_elas_barra2D (eps_new,hvar_old,Eprop,ce)

%******************************************************************************************
%*  RETTURN-MAPPING (TENSOR DE TERNSIONES) PARA BARRAS EN 2D                              *
%*  MODELO ELASTICO                                                                       *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

% calculo de la tension
% иииииииииииииииииииии
sigma_new = ce*eps_new;

% Variables historicas
% ииииииииииииииииииии
hvar_new  = hvar_old;

% Variables auxiliares
% ииииииииииииииииииии
aux_var = 0;