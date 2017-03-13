function [sigma_new,hvar_new,aux_var] = rmapfi_danio_plasdp(eps_new,hvar_old,Eprop,ce,e_VG)

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
   hvar_old(ntens+3,1) = sigmay/sqrt(E);
   hvar_old(ntens+2,1) = sigmay/sqrt(E);
end
eps_old_plast = hvar_old(1:ntens,1);
alpha_old     = hvar_old(ntens+1,1);
r_old         = hvar_old(ntens+2,1);
q_old         = hvar_old(ntens+3,1);
alpha_0       = hvar_old(ntens+4,1);

% Constantes elasticas
% จจจจจจจจจจจจจจจจจจจจ
mu        = E/(2*(1+poiss));
lambda    = (poiss*E)/((1+poiss)*(1-2*poiss)); 
cappa     = E/(3*(1-2*poiss));


kapa  = 2*mu*sigmay^2/E;
chi_max = (2*tita-6*mu/E/tita)*sigmay;
chi_0   = chi_max/3;

% este nivel de hardening es para un :: chi_max = chi_0 
% (fc corresponde a chi_max)
coe1 = 100/3;  
coe2 = 1000/3;

 if (alpha_old <= 0)
     a= 0;
     b= 0;
     c= chi_0;
 elseif (alpha_old <= coe1*alpha_0);
     a= -2*chi_0 / (coe1*alpha_0)^2; 
     b=  4*chi_0 / (coe1*alpha_0)  ;
     c=  chi_0;
 elseif (alpha_old < coe2*alpha_0)
     a= -3*chi_0 /(coe2-coe1)^2 /alpha_0^2;
     b=  6*chi_0*coe1/(coe2-coe1)^2 /alpha_0;
     c=  3*chi_0*coe2*(coe2-2*coe1) / (coe2-coe1)^2;
 else
     a= 0;
     b= 0;
     c= 0;
 end
 chi_n = a* alpha_old ^2 + b*alpha_old + c ;    

 
 
if (chi_n <= 0)
    chi_n = 0;
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
rho   = 0;
fact1 = 0;
fact2 = 0;
fact3 = 0;
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
   alpha_new     = alpha_old;
   
else

    % Modelo de plasticidad de Drucker-Prager
    s_trial      = FODPT*sigma_trial;
    norm_s_trial = sqrt(s_trial.'*SONT*s_trial);
    f_trial = norm_s_trial^2  + chi_n * sigma_m_trial - kapa ;
    
    if (f_trial > 0)
        
        xi1= (2*mu*q_old/r_old)^2+a*sigma_m_trial;
        xi2= -4*mu*q_old/r_old*norm_s_trial+ (2*a*alpha_old+b)*sigma_m_trial ;
        xi3= f_trial;
        discri= xi2^2-4*xi1*xi3;
        if discri >= 0
            gamma= (-xi2- sqrt(discri))/2/xi1;
            gamma1= (-xi2 + sqrt(discri))/2/xi1;
        else
              cccccc=22222
        end
        
        rho                  = (1- q_old/r_old*2*mu*gamma/norm_s_trial);
        s_new                =  rho * s_trial;
        sigma_new            = s_new + sigma_m_trial * SOIT;
        delta_eps_plast      = gamma *s_trial/norm_s_trial;  %acordarse de corregir el termino cortante
        norm_delta_eps_plast = sqrt(delta_eps_plast.'*SONT*delta_eps_plast);
        eps_new_plast        = eps_old_plast + delta_eps_plast;
        q_new                = q_old;
        r_new                = r_old;
        alpha_new            = alpha_old + norm_delta_eps_plast;
        if (alpha_0 == 0)
%            alpha_0 = sqrt(eps_elas.'*SONT*eps_elas); %ojo
            alpha_0 = (sigmay/E)/ (tita/3) *((3*abs(sigma_m_trial)/(sigmay*(tita/3))))^1.5 ; %ojo
%            alpha_0 = sigmay/E ; %ojo
        end
        fload  = -1;
        den    =  2*gamma*xi1+xi2;
        coe0   = (q_old/r_old);
        fact1  = (1- 2*mu*gamma*coe0/norm_s_trial);
        fact2  = 2*mu*coe0*( gamma/norm_s_trial + 2*(norm_s_trial-2*mu*coe0*gamma)/den ) ;
        fact3  = (a*gamma^2+ gamma*(2*a*alpha_old+b)+ chi_n) /den ;
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
hvar_new(ntens+1,1) = alpha_new;
hvar_new(ntens+2,1) = r_new;
hvar_new(ntens+3,1) = q_new;
hvar_new(ntens+4,1) = alpha_0;


% Variables auxiliares
% จจจจจจจจจจจจจจจจจจจจ
aux_var(1,1) = fload;
aux_var(2,1) = q_new;
aux_var(3,1) = r_new;
aux_var(4,1) = fact1;
aux_var(5,1) = fact2;
aux_var(6,1) = fact3;


