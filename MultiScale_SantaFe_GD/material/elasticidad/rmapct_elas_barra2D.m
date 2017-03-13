function [ct] = rmapct_elas_barra2D (eps_new,hvar_old,aux_var,Eprop,ce)

%******************************************************************************************
%*  RETTURN-MAPPING MATRIZ TANGENTE PARA BARRAS EN 2D                                     *
%*  MODELO ELASTICO                                                                       *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

% Calculo de la matriz tangente
% иииииииииииииииииииииииииииии
ct = ce;