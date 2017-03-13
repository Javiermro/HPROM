function [ct,sigma_new,eps_new,hvar_new,aux_var] = rmap_plasJ2(eps_new,hvar_old,e_DatMatSet,e_VG)

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
factor          = 1e-6;

% Recupera propiedades del material
% E             = Eprop(4);
% poiss         = Eprop(5);
% sigmay        = Eprop(6);
% tita          = Eprop(8);
% hb            = Eprop(9);
% kinfb         = Eprop(10);
% kcerob        = Eprop(11);
% delta         = Eprop(12);
E = e_DatMatSet.young;
poiss = e_DatMatSet.poiss;
sigmay = e_DatMatSet.ftult;
tita = e_DatMatSet.tit;
hb = e_DatMatSet.hba;
kinfb = e_DatMatSet.kin;
kcerob = e_DatMatSet.kce;
delta = e_DatMatSet.del;
ce = e_DatMatSet.ce;

mu            = E/(2*(1+poiss));
cappa         = E/(3*(1-2*poiss));

% Recupera variables internas
eps_old_plast = hvar_old(1:ntens);
alpha_old     = hvar_old(ntens+1);

% Cálculo del estado trial
sigma_trial = ce*(eps_new - eps_old_plast);
alpha_trial = alpha_old;
s_trial     = FODPT*sigma_trial;

% Evaluación del criterio de fluencia
norm_s_trial = sqrt(s_trial.'*SONT*s_trial);
k_old        = sigmay + (tita*hb*alpha_old) + (kinfb-kcerob)*(1-exp(-delta*alpha_old));

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

if (f_trial <= 0 || e_VG.elast>0)
   
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
   if ((kinfb~=kcerob) && (delta~=0))
      % Softening exponencial
      [dgamma] = get_dgamma (alpha_old,norm_s_trial,e_DatMatSet,e_VG);
   else
      % Softening lineal
      dgamma = f_trial/(1+(hb/3/mu))/2/mu; 
   end
   
   % Actualizar sigma_new, eps_new_plast, alpha_new
   sigma_new     = cappa*(SOIT.'*eps_new)*SOIT + s_trial - 2*mu*dgamma*N_new;
   eps_new_plast = eps_old_plast + dgamma*SONT*N_new;
   alpha_new     = alpha_old + sqrt(2/3)*dgamma;   
   
   % Actualizar derivada del modulo de hardening-softening
   kp_new = (tita*hb)+(kinfb-kcerob)*(delta*exp(-delta*alpha_new));
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