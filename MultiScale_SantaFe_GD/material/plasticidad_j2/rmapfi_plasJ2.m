function [sigma_new,hvar_new,aux_var] = rmapfi_plasJ2 (eps_new,hvar_old,Eprop,ce,e_VG)

%******************************************************************************************
%*  RETTURN-MAPPING PARA COMPUTO DE FUERZA INTERNA: PLANE STRAIN - 3D                     *
%*  MODELO DE PLASTICIDAD J2 CON ENDURECIMIENTO ISOTROPO                                  *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

%global sihvarpg FODPT SONT SOIT ntens
sihvarpg = e_VG.sihvarpg;
FODPT = e_VG.FODPT;
SOIT = e_VG.SOIT;
SONT = e_VG.SONT;
ntens = e_VG.ntens;
sihvarpg = e_VG.sihvarpg;

% Variables
% จจจจจจจจจ
eps_old_plast = hvar_old(1:ntens,1);
alpha_old     = hvar_old(ntens+1,1);
E             = Eprop(4);
poiss         = Eprop(5);
sigmay        = Eprop(6);
tita          = Eprop(8);
hb            = Eprop(9);
kinfb         = Eprop(10);
kcerob        = Eprop(11);
delta         = Eprop(12);

% Constantes elasticas
% จจจจจจจจจจจจจจจจจจจจ
mu     = E/(2*(1+poiss));
lambda = (poiss*E)/((1+poiss)*(1-2*poiss)); 
cappa  = E/(3*(1-2*poiss));

% Cแlculo del estado trial
% จจจจจจจจจจจจจจจจจจจจจจจจ
sigma_trial = ce*(eps_new - eps_old_plast);
alpha_trial = alpha_old;
s_trial     = FODPT*sigma_trial;

% Evaluaci๓n del criterio de fluencia
% จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
norm_s_trial = sqrt(s_trial.'*SONT*s_trial);
k_old        = sigmay + (tita*hb*alpha_old) + (kinfb-kcerob)*(1-exp(-delta*alpha_old));
f_trial      = norm_s_trial - sqrt(2/3)*k_old;

if f_trial <= 0;
   
   % Paso elแstico
   % จจจจจจจจจจจจจ
 
   % Actualiza sigma_new, eps_new_plast, alpha_new
   sigma_new     = sigma_trial;
   eps_new_plast = eps_old_plast;
   alpha_new     = alpha_trial;
   
else
   
   % Paso plแstico
   % จจจจจจจจจจจจจ
  
   % Cแlculo de N_new
   N_new = s_trial/norm_s_trial;
   
   % Cแlculo de dgamma, h_old, h_new, hp_new (solo para hardening-softening isotropo lineal)
   dgamma = f_trial/(1+(hb/3/mu))/2/mu;
   h_old  = 0; 
   h_new  = 0; 
   hp_new = 0;
   
   % Actualizar sigma_new, eps_new_plast, alpha_new
   sigma_new     = cappa*(SOIT.'*eps_new)*SOIT + s_trial - 2*mu*dgamma*N_new;
   eps_new_plast = eps_old_plast + dgamma*SONT*N_new;
   alpha_new     = alpha_old + sqrt(2/3)*dgamma;
   
end

% Variables historicas
% จจจจจจจจจจจจจจจจจจจจ
hvar_new            = zeros(sihvarpg,1);
hvar_new(1:ntens,1) = eps_new_plast;
hvar_new(ntens+1,1) = alpha_new;

% Variables auxiliares
% จจจจจจจจจจจจจจจจจจจจ
aux_var = 0;