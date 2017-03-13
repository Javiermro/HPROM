function ct = rmapct_plasJ2 (eps_new,hvar_old,aux_var,Eprop,ce,e_VG)

%******************************************************************************************
%*  RETTURN-MAPPING PARA COMPUTO DEL TENSOR TANGENTE: PLANE STRAIN - 3D                   *
%*  MODELO DE PLASTICIDAD J2 CON ENDURECIMIENTO ISOTROPO                                  *                 
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

%global FODPT SONT ntens
FODPT = e_VG.FODPT;
SONT = e_VG.SONT;
ntens = e_VG.ntens;

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

% M๓dulo elastoplแstico consistente
% จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
if f_trial <= 0;
   
   ct = ce;
   
else
   
   % Cแlculo de N_new
   N_new = s_trial/norm_s_trial;
   
   % Cแlculo de dgamma, h_old, h_new, hp_new (solo para hardening-softening isotropo lineal)
   dgamma = f_trial/(1+(hb/3/mu))/2/mu;
   h_old  = 0; 
   h_new  = 0; 
   hp_new = 0;
   
   % Actualizar alpha_new
   alpha_new = alpha_old + sqrt(2/3)*dgamma;

   % Actualizar derivada del modulo de hardening-softening
   kp_new = (tita*hb)+(kinfb-kcerob)*(delta*exp(-delta*alpha_new));

   % Computar tensor tangente
   ct = ctpg_plasJ2(dgamma,norm_s_trial,kp_new,hp_new,N_new,mu,cappa,e_VG);
   
end