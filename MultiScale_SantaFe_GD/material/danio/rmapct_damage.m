function [ct] = rmapct_damage(eps_n1,aux_var,ce)

%*********************************************************************************
%*  cálculo del tensor constitutivo tangente ct                                  *
%*  en el modelo de daño simétrico                                               *
%*********************************************************************************

%global FODPT SONT

% Variables
% ¨¨¨¨¨¨¨¨¨
fload= aux_var(1) ;
factor1 = aux_var(2);
factor2 = aux_var(3);

% Módulo elastoplástico consistente
% ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
sigma_efectiva= ce*eps_n1;
ct = factor1*ce;
if fload > 0;
  ct=ct - factor2*sigma_efectiva* sigma_efectiva' ;    
end
return
