function [sigma_new,hvar_new,aux_var] = rmapfi_plasJ2_barra2D (eps_new,hvar_old,Eprop,ce,e_VG)

%******************************************************************************************
%*  RETTURN-MAPPING FUERZA INTERNA PARA BARRAS EN 2D                                      *
%*  MODELO DE PLASTICIDAD J2 CON ENDURECIMIENTO ISOTROPO LINEAL                           *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

%global ntens sihvarpg
ntens = e_VG.ntens;
sihvarpg = e_VG.sihvarpg;

% Variables
% จจจจจจจจจ
eps_old_plast  = hvar_old(1:ntens,1);
alpha_old      = hvar_old(ntens+1,1);
E              = Eprop(4);
sigmay         = Eprop(6);
tita           = Eprop(8);
hb             = Eprop(9);
kinfb          = Eprop(10);
kcerob         = Eprop(11);
delta          = Eprop(12);

% Cแlculo del estado trial
% จจจจจจจจจจจจจจจจจจจจจจจจ
eps_trial_plast = eps_old_plast;
alpha_trial     = alpha_old;
sigma_trial     = ce*(eps_new - eps_old_plast);

% Evaluaci๓n del criterio de fluencia
% จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
k         = hb;
f_trial   = abs(sigma_trial) - (sigmay + k*alpha_trial);

if f_trial <= 0
   
   % Paso elแstico
   % จจจจจจจจจจจจจ
   sigma_new      = sigma_trial;
   eps_new_plast  = eps_trial_plast;
   alpha_new      = alpha_trial;
   
else
   
   % Paso plแstico
   % จจจจจจจจจจจจจ
   dgamma        = f_trial/(E+k);
   sigma_new     = (1-dgamma*E/abs(sigma_trial))*sigma_trial;
   eps_new_plast = eps_old_plast + dgamma*sign(sigma_trial);
   alpha_new     = alpha_old + dgamma;
   
end

% Variables historicas
% จจจจจจจจจจจจจจจจจจจจ
hvar_new            = zeros(sihvarpg,1);
hvar_new(1:ntens,1) = eps_new_plast;
hvar_new(ntens+1,1) = alpha_new;

% Variables auxiliares
% จจจจจจจจจจจจจจจจจจจจ
aux_var = 0;