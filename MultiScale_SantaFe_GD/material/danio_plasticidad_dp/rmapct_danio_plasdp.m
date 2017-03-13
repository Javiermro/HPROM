function ct = rmapct_danio_plasdp(eps_new,hvar_old,aux_var,Eprop,ce,e_VG)

%******************************************************************************************
%*  RETTURN-MAPPING PARA COMPUTO DEL TENSOR TANGENTE: PLANE STRAIN - 3D                   *
%*  MODELO DE:                                                                            *
%*  - DAัO ISOTROPO CON ABLANDAMIENTO PARA TENSION MEDIA POSITIVA                         *
%*  - PLASTICIDAD PERFECTA DE DUCKER PRAGER PARA TENSION MEDIA NEGATIVA                   *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

%global FODPT FOAT2 SOIT SONT SSOIT ntens istep
FODPT = e_VG.FODPT;
FOAT2 = e_VG.FOAT2;
SOIT = e_VG.SOIT;
SONT = e_VG.SONT;
SSOIT = e_VG.SSOIT;
ntens = e_VG.ntens;
istep = e_VG.istep;

% Propiedades del material
% จจจจจจจจจจจจจจจจจจจจจจจจ
E      = Eprop(4);
poiss  = Eprop(5);
sigmay = Eprop(6);
tita   = Eprop(8);
hb     = Eprop(9);
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

% Variables auxiliares
% จจจจจจจจจจจจจจจจจจจจ
fload     = aux_var(1,1);
q_new     = aux_var(2,1);
r_new     = aux_var(3,1);
fact1     = aux_var(4,1);
fact2     = aux_var(5,1);
fact3     = aux_var(6,1);

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

if (fload == 0)

   % Elasticidad
   ct = (q_new/r_new) * ce;

elseif (fload == 1)

   % Modelo de da๑o
   coe1 = (q_new/r_new);
   coe2 = (hb*r_old-q_old)/(r_new^3);
   coe3 = 0;
   ct   = coe1*ce + (coe2+coe3)*(sigma_eff*sigma_eff.');
   
elseif (fload == -1)

   % Modelo de plasticidad de Drucker-Prager
   sigma_m_trial = (SOIT.' * sigma_trial)/3;
   s_trial       = FODPT*sigma_trial;
   norm_s_trial  = sqrt(s_trial.'*SONT*s_trial);
   ntrial        = s_trial/norm_s_trial;

   coe0 = (q_old/r_old);
   coe1 =  coe0*cappa;
   coe2 =  2*mu*(q_old/r_old) ;
      
   ct   = coe1*SSOIT + coe2*fact1*FOAT2 + coe2*fact2*(ntrial*ntrial.')+coe2*fact3*coe1*(ntrial*SOIT.');
   
end