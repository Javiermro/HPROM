function [ct,sigma_n1,sigma_impl_n1,hvar_n1] = ...    %,aux_var
    rmap_damage_elastic (eps_n1,hvar_n,e_DatMatSet)

%*********************************************************************************
%*                     CURSO DE MODELOS CONSTITUTIVOS - UPC -                    *
%*           Algoritmo de integración para el modelo de daño isótropo            *
%*                                                                               *
%* Llamada :: [sigma_n1,hvar_n1,aux_var] = rmap_dano (eps_n1,hvar_n,Eprop,ce)    *
%*                                                                               *
%* ARGUMENTOS ENTRADA: eps_n1(4)   deformación total , paso n+1                  *
%*                                 vector R4    (exx eyy exy ezz)                *
%*                     hvar_n(5)   variables internas , paso n                   *
%*                                 hvar_n(1:4) deformacion total paso n          *
%*                                 hvar_n(5) = r  ; hvar_n(6)=q                  *
%*                     Eprop(4)    vector de propiedades de material             *
%*                                 Eprop(1)=  E     : modulo de Young            *
%*                                 Eprop(2)=  nu    : modulo de Poisson          *
%*                                 Eprop(3)=  H     : modulo de Softening/hard.  *
%*                                 Eprop(4)=sigma_u : tensión última             *
%*                     ce(4,4)     tensor constitutivo elástico                  *
%*                                                                               *
%* ARGUMENTOS SALIDA:  sigma_n1(4) tensión de cauchy , paso n+1                  *
%*                     hvar_n(5)   variables internas , paso n+1                 *
%*                     aux_var(3)  variables auxiliares para calculo tensor tang.*
%*********************************************************************************

ce       = e_DatMatSet.ce;
%r0       = e_DatMatSet.r_0;
esImplex = e_DatMatSet.esImplex;

% Recuperacion de variables internas del paso anterior
% ****************************************************
r_n  = hvar_n(5,1);
q_n  = hvar_n(6,1);

% Inicializacion de las variables
% *******************************
%Estas variables se asumen inicializadas con r0 (ver función f_eVarEstInic).
% if(r_n<=0.d0)
%     r_n=r0;
%     q_n=r0;
% end

% ESTADO ELASTICO (O DE DESCARGA ELASTICA)
% ****************************************
%tension_efectiva=ce*eps_n1;
ct = q_n/r_n*ce;
sigma_n1 = ct*eps_n1;
sigma_impl_n1 = sigma_n1;

% Actualización de variables historicas
% *************************************
hvar_n1  = hvar_n ;
hvar_n1(3,1)    = 0  ; % r_n1-r_n ;  

if esImplex
   ct = struct('Implex',ct,'Impli',ct);
end

end