function [beta_impli, mod_beta,factor_i,der_chi] = newton_cohesive_ds(qbifu,...
                                    mod_beta,gfval,qfi_0,tracc, beta_impli)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ***  Newton  method for solving the cohesive force in the model DS  ***
% input variables:
%
%
% output vbariables:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
iterMax = 50;
facRel = 1.e-10;
% Iteración cero
iter = 0;
factor_i   = qbifu/mod_beta*exp(-qbifu*mod_beta/gfval);
der_factor = (-qbifu^2/gfval* exp(-qbifu*mod_beta/gfval));
der_chi    = (der_factor  - factor_i ) / mod_beta^2;
qfi        = factor_i * eye(2,2)+qfi_0;
res_b      = qfi*beta_impli -tracc;
mod_res=norm(res_b);
normt      = norm(tracc);
%
%disp('--------------------------------------')
%fprintf('SDA - Exponencial: Iter %d |res| = %e\n',iter,mod_res)
%
while(mod_res>facRel*normt&&iter<iterMax)
   iter = iter+1;
   k          = qfi+der_chi*(beta_impli* beta_impli');
   delta_b    = -k\res_b;
   beta_impli = beta_impli+delta_b;
   mod_beta   = norm(beta_impli);
   expone     = exp(-qbifu*mod_beta/gfval);
   factor_i   = (qbifu/mod_beta)*expone; 
   der_factor = (-qbifu^2/gfval)*expone; 
   der_chi    = (der_factor  - factor_i  ) / mod_beta^2;
   qfi        = factor_i * eye(2,2)+qfi_0;
   res_b      = qfi*beta_impli -tracc; 
   mod_res    = norm(res_b);
   %
   %fprintf('SDA - Exponencial: Iter %d |res| = %e\n',iter,mod_res)
end

if(mod_res>facRel*normt)   
   %fprintf('Modelo de fuerzas centradas: Newton sobre salto: Módulo constitutivo no converge\n');
   warning('Modelo de fuerzas centradas: Newton sobre salto: Módulo constitutivo no converge') %#ok<WNTAG>
end
%
%disp('--------------------------------------')