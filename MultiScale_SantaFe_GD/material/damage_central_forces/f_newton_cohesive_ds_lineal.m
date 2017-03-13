function [beta_impli,mod_beta,factor_i,der_chi] = f_newton_cohesive_ds_lineal(qbifu,...
                                    mod_beta,gfval,qfi_0,tracc,beta_impli)  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ***  Newton  method for solving the cohesive force in the model DS  ***
% input variables:
%
%
% output vbariables:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
iterMax = 100;
facRel = 1.e-10;
% Iteración cero
iter=0;
factor_i = (1/mod_beta-qbifu/2/gfval)*qbifu;
der_chi = -qbifu/mod_beta^3;
qfi = factor_i*eye(2,2)+qfi_0;
res_b = qfi*beta_impli-tracc;
mod_res = norm(res_b);
normt = norm(tracc);
%
%disp('--------------------------------------')
%fprintf('SDA - Lineal: Iter %d |res| = %e\n',iter,mod_res)
%
while(mod_res>facRel*normt&&iter<iterMax)
   iter = iter+1;
   k = qfi+der_chi*(beta_impli*beta_impli');
   delta_b = -k\res_b;
   beta_impli = beta_impli+delta_b;
   mod_beta = norm(beta_impli);
   factor_i = (1/mod_beta-qbifu/2/gfval)*qbifu; 
   der_chi = -qbifu/mod_beta^3;
   qfi = factor_i*eye(2,2)+qfi_0;
   res_b = qfi*beta_impli-tracc; 
   mod_res = norm(res_b);
   %
   %fprintf('SDA - Lineal: Iter %d |res| = %e\n',iter,mod_res)
end
%if iter>1
%   fprintf('Hizo más de una iteración (%d).\n',iter)
%   pause(1)
%end
if(mod_res>facRel*normt)
   %fprintf('Modelo de fuerzas centradas: Newton sobre salto: Módulo constitutivo no converge\n');
   warning('Modelo de fuerzas centradas: Newton sobre salto: Módulo constitutivo no converge') %#ok<WNTAG>
end
%
%disp('--------------------------------------')