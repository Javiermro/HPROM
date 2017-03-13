function [ct,sigma_new,eps_new,hvar_new,aux_var] = rmap_plasJ2bl (eps_new,hvar_old,Eprop,ce,e_VG)

%******************************************************************************************
%*  RETTURN-MAPPING PARA COMPUTO DE FUERZA INTERNA Y MATRIZ TANGENTE: PLANE STRAIN - 3D   *
%*  MODELO DE PLASTICIDAD J2 CON ENDURECIMIENTO/ABLANDAMIENTO ISOTROPO                    *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

% Variables globales
FODPT           = e_VG.FODPT;
SOIT            = e_VG.SOIT;
SSOIT           = e_VG.SSOIT;
SONT            = e_VG.SONT;
ntens           = e_VG.ntens;
Omega_micro     = e_VG.Omega_micro;
Omega_micro_loc = e_VG.Omega_micro_loc;
factor          = 1e-6;

% Recupera propiedades del material
E             = Eprop(4);
poiss         = Eprop(5);
sigmay        = Eprop(6);
tita          = Eprop(8);
hb            = Eprop(9);
alpha_lim     = Eprop(10); % kinfb
hbb           = Eprop(11); % kcerob
chi_factor    = Eprop(12); % delta
mu            = E/(2*(1+poiss));
cappa         = E/(3*(1-2*poiss));
chi           = 0.0;

% Recupera variables internas
eps_old_plast = hvar_old(1:ntens);
alpha_old     = hvar_old(ntens+1);

% Cálculo del estado trial
sigma_trial = ce*(eps_new - eps_old_plast);
alpha_trial = alpha_old;
s_trial     = FODPT*sigma_trial;

% Evaluación del criterio de fluencia
norm_s_trial = sqrt(s_trial.'*SONT*s_trial);
if (alpha_old > alpha_lim)
   chi = chi_factor;
end
k_old = sigmay + (tita*hb*alpha_old) + (chi*hbb*(alpha_old-alpha_lim));

if (k_old < factor*sigmay)
    sigma_new     = cappa*(SOIT.'*eps_new)*SOIT + factor*s_trial;
    eps_new_plast = eps_old_plast;
    alpha_new     = alpha_old;
    ct            = (cappa*SSOIT) + factor*FODPT*ce;
    fload_new     = 1;
    hvar_new      = [eps_new_plast ; alpha_new ; fload_new ; sigma_new];
    aux_var       = fload_new;
    return;
end

f_trial = norm_s_trial - sqrt(2/3)*k_old;

condition = 0;
%if ((e_VG.iElem == 17 || e_VG.iElem == 18 || e_VG.iElem == 19 || e_VG.iElem == 20) && ...
%        (e_VG.istep_bif ~= 0 && e_VG.istep >= e_VG.istep_bif))
%if ((e_VG.iElem == 5) && (e_VG.istep_bif ~= 0 && e_VG.istep >= e_VG.istep_bif))
%  condition = 1;
%end

if (f_trial <= 0 | condition == 1)
   
   % Paso elástico
 
   % Actualiza sigma_new, eps_new_plast, alpha_new
   sigma_new     = sigma_trial;
   eps_new_plast = eps_old_plast;
   alpha_new     = alpha_trial;
   
   % Computar tensor tangente
   ct = ce;
   
   % Indice de carga
   fload_new = 0;
   
else
   
   % Paso plástico
  
   % Cálculo de N_new
   N_new = s_trial/norm_s_trial;
   
   % Cálculo de dgamma, h_old, h_new, hp_new (solo para hardening-softening isotropo lineal)
   dgamma = f_trial/(1+(tita*hb+chi*hbb)/3/mu)/(2*mu);
   
   % Actualizar sigma_new, eps_new_plast, alpha_new
   sigma_new     = cappa*(SOIT.'*eps_new)*SOIT + s_trial - 2*mu*dgamma*N_new;
   eps_new_plast = eps_old_plast + dgamma*SONT*N_new;
   alpha_new     = alpha_old + sqrt(2/3)*dgamma;   
   
   % Actualizar derivada del modulo de hardening-softening
   kp_new = (tita*hb)+(chi*hbb);
   hp_new = 0;

   % Computar tensor tangente
   ct = ctpg_plasJ2(dgamma,norm_s_trial,kp_new,hp_new,N_new,mu,cappa,e_VG);
   
   % Indice de carga
   fload_new = 1;
   
end

% Variables historicas
hvar_new = [eps_new_plast ; alpha_new ; fload_new ; sigma_new];

% Variables auxiliares
aux_var = fload_new;