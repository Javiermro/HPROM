function [sigma_new,hvar_new,aux_var] = rmapfi_danio_plasdp (eps_new,hvar_old,Eprop,ce,e_VG)

%******************************************************************************************
%*  RETTURN-MAPPING PARA COMPUTO DE FUERZA INTERNA: PLANE STRAIN - 3D                     *
%*  MODELO DE:                                                                            *
%*  - DAัO ISOTROPO CON ABLANDAMIENTO PARA TENSION MEDIA POSITIVA                         *
%*  - PLASTICIDAD PERFECTA DE DUCKER PRAGER PARA TENSION MEDIA NEGATIVA                   *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

%global sihvarpg FODPT SONT SOIT ntens istep
sihvarpg = e_VG.sihvarpg;
FODPT = e_VG.FODPT;
SOIT = e_VG.SOIT;
SONT = e_VG.SONT;
ntens = e_VG.ntens;
sihvarpg = e_VG.sihvarpg;
istep = e_VG.istep;

% Propiedades del material
% จจจจจจจจจจจจจจจจจจจจจจจจ
E      = Eprop(4);
poiss  = Eprop(5);
sigmay = Eprop(6);
tita   = Eprop(8);     % Factor de sobre-resistencia a compresion
hb     = Eprop(9);     % Modulo ablandamiento para da๑o
kinfb  = Eprop(10);
kcerob = Eprop(11);
delta  = Eprop(12);

% Variables internas
% จจจจจจจจจจจจจจจจจจ
if (istep==1)
   hvar_old(ntens+1,1) = sigmay/sqrt(E);
   hvar_old(ntens+2,1) = sigmay/sqrt(E);
end
eps_old_plast = hvar_old(1:ntens,1);
q_old         = hvar_old(ntens+1,1);
r_old         = hvar_old(ntens+2,1);
alpha_old     = hvar_old(ntens+3,1);
alpha_max     = hvar_old(ntens+4,1);

% Constantes elasticas
% จจจจจจจจจจจจจจจจจจจจ
mu        = E/(2*(1+poiss));
lambda    = (poiss*E)/((1+poiss)*(1-2*poiss)); 
cappa     = E/(3*(1-2*poiss));
chi_max   = 3*(sqrt(2/3)-sqrt(2*mu/E)/tita);
chi_0     = 3*(sqrt(2/3)-sqrt(2*mu/E)/tita*3);
%alpha_max = 4*(5/3)*tita*sigmay/E;
if (alpha_max == 0)
    theta = 0;
else
    theta     = (chi_max-chi_0)/(alpha_max)^2;
end
if (alpha_old < alpha_max)
    chi       = chi_0 + 2*theta*alpha_max*alpha_old - theta*alpha_old^2;
else
    chi = chi_max - chi_max*(alpha_old-alpha_max)/(100*alpha_max-alpha_max);
end

if (alpha_max == 0)
    chi       = chi_0;
end

if (chi <= 0.0)
   chi = 0.0;
end

% Cแlculo del estado trial y tensi๓n efectiva
% จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
eps_elas      = (eps_new - eps_old_plast);
sigma_eff     = ce * eps_elas;
sigma_trial   = (q_old/r_old) * sigma_eff;
sigma_m_trial = (SOIT.' * sigma_trial)/3;

% Variable auxiliar
% จจจจจจจจจจจจจจจจจ
fload = 0;

if (sigma_m_trial >= 0)
   
   % Modelo de da๑o
   
   tau_epsilon = sqrt(eps_elas.' * ce * eps_elas);

   if (tau_epsilon > r_old)
      r_new = tau_epsilon;
      q_new = q_old + hb * (r_new - r_old);
      fload = 1;
   else
      r_new = r_old;
      q_new = q_old;
   end
   
   eps_new_plast = eps_old_plast;
   sigma_new     = (q_new/r_new) * sigma_eff;
   
else

    % Modelo de plasticidad de Drucker-Prager

    s_trial      = FODPT*sigma_trial;
    norm_s_trial = sqrt(s_trial.'*SONT*s_trial);
    f_trial      = chi * sigma_m_trial + norm_s_trial - sqrt(2*mu) * q_old;
    
    if (f_trial > 0)
        sigma_m_new          = sigma_m_trial;
        norm_s_new           = sqrt(2*mu) * q_old - chi * sigma_m_new;
        s_new                = norm_s_new * (s_trial/norm_s_trial);
        sigma_new            = s_new + sigma_m_new * SOIT;
        delta_eps_plast      = (s_trial - s_new)/(2*mu);
        norm_delta_eps_plast = sqrt(delta_eps_plast.'*SONT*delta_eps_plast);
        eps_new_plast        = eps_old_plast + delta_eps_plast;
        q_new                = q_old;
        r_new                = r_old;
        alpha_new            = alpha_old + norm_delta_eps_plast;
        if (alpha_max == 0)
            alpha_max = (5/3)*norm_s_trial/2/mu;
%           alpha_max = (5/3)*tita*sigmay/E + 0.0*abs(sigma_m_trial)/cappa;
        end
        fload                = -1;
    else
        sigma_new       = sigma_trial;
        eps_new_plast   = eps_old_plast;
        q_new           = q_old;
        r_new           = r_old;
        alpha_new       = alpha_old;
    end

end

% Actualizaci๓n variables internas
% จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
hvar_new            = zeros(sihvarpg,1);
hvar_new(1:ntens,1) = eps_new_plast;
hvar_new(ntens+1,1) = q_new;
hvar_new(ntens+2,1) = r_new;
hvar_new(ntens+3,1) = alpha_new;
hvar_new(ntens+4,1) = alpha_max;


% Variables auxiliares
% จจจจจจจจจจจจจจจจจจจจ
aux_var(1,1) = fload;
aux_var(2,1) = q_new;
aux_var(3,1) = r_new;
aux_var(4,1) = alpha_new;


