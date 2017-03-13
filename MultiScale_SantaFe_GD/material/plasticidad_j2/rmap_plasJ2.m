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
esImplex = e_DatMatSet.esImplex;

mu            = E/(2*(1+poiss));
cappa         = E/(3*(1-2*poiss));

% Recupera variables internas
eps_old_plast = hvar_old(1:ntens);
alpha_old     = hvar_old(ntens+1);

% C�lculo del estado trial
sigma_trial = ce*(eps_new - eps_old_plast);
alpha_trial = alpha_old;
s_trial     = FODPT*sigma_trial;

% Evaluaci�n del criterio de fluencia
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
   
   % Paso elastico
 
   % Actualiza sigma_new, eps_new_plast, alpha_new
   sigma_new     = sigma_trial;
   eps_new_plast = eps_old_plast;
   alpha_new     = alpha_trial;
   
   % Computar tensor tangente
   ct = ce;
   
   % Indice de carga
   fload_new = 0;
   
else
   
   % Paso plastico
  
   % Calculo de N_new
   N_new = s_trial/norm_s_trial;
   
   % C�lculo de dgamma, h_old, h_new, hp_new (solo para hardening-softening isotropo lineal)
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

% for computing the internal energy
%k_new = sigmay + (tita*hb*alpha_new) + (kinfb-kcerob)*(1-exp(-delta*alpha_new));
%k_new = sigmay*alpha_new+0.5*(tita*hb)*alpha_new^2+(kinfb-kcerob)*alpha_new+(1/delta)*(exp(-delta*alpha_new)-1);
%k_new = sigmay*alpha_new+0.5*(tita*hb)*alpha_new^2+(kinfb-kcerob)*(alpha_new+(1/delta)*(exp(-delta*alpha_new)-1));
%k_new = 0.5*(tita*hb)*alpha_new^2;
%k_new = sigmay*alpha_new + 0.5*(tita*hb)*alpha_new^2;

k_new = -(tita*hb*alpha_new) - (kinfb-kcerob)*(1-exp(-delta*alpha_new));

%
s_new     = FODPT*sigma_new;
norm_s_new = sqrt(3/2)*sqrt(s_new.'*SONT*s_new);

%Computo variables explicitas (Esquema IMPLEX)
if esImplex
   % bandwidth = ksd; %if bandwidth==0; bandwidth=1; end
    Dtime    = e_VG.Dtime;
    if ~e_VG.elast
    %if (f_trial > 0 && ~e_VG.elast)
        
        %dgamma_old    = hvar_old(12);
        dgamma_old  = hvar_old(13);  if dgamma_old==0; dgamma=0; end
        %bandwidth_pre = hvar_old(6); if bandwidth_pre==0; bandwidth_pre=bandwidth; end
        Dtime_n        = hvar_old(14); if Dtime_n==0; Dtime_n=Dtime; end
        
        %dgamma_Implex=(bandwidth_pre/bandwidth)*(Dtime/Dtime_n)*delta_r_old;
        dgamma_Implex = (Dtime/Dtime_n)*dgamma_old;
        
        % Actualizar derivada del modulo de hardening-softening
        %kp_old = (tita*hb)+(kinfb-kcerob)*(delta*exp(-delta*alpha_old));        
        %q_tilde = k_old + kp_old*dgamma_Implex;
               
        % Actualizar sigma_new, eps_new_plast, alpha_new
        %sigma_newImplex     = cappa*(SOIT.'*eps_new)*SOIT + s_trial - 2*mu*dgamma_Implex*N_new;
        
        %rad_Mises = sqrt(2/3)*(q_tilde);
        rad_Mises = sqrt(2/3)*(k_old);        
        fact = (rad_Mises/(rad_Mises+2*mu*dgamma_Implex));

        % Computar tensiones IMPLEX
        sigma_newImplex = cappa*(SOIT.'*eps_new)*SOIT + fact*s_trial;
        
        % Computar tensor tangente IMPLEX
        ctImplex = ctpg_plasJ2_implex(mu,cappa,fact,e_VG);
        
        % Se reemplaza ya que no se utiliza las sigma impl�cita por las implex en el argumento de salida.
        sigma_new = sigma_newImplex;
        
        %Se almacena los resultados.
        ct = struct('Implex',ctImplex,'Impli',ct);
    else
        dgamma=0.0;
        ct = struct('Implex',ct,'Impli',ct);
    end
   
   % Variables historicas
   hvar_new = [eps_new_plast ; alpha_new ; fload_new ; sigma_new; norm_s_new; k_new; dgamma; Dtime];
   
else
   % Variables historicas
   hvar_new = [eps_new_plast ; alpha_new ; fload_new ; sigma_new; norm_s_new; k_new];   
end

% Variables historicas
%hvar_new = [eps_new_plast ; alpha_new ; fload_new ; sigma_new; norm_s_new];

% Variables auxiliares
aux_var = fload_new;
