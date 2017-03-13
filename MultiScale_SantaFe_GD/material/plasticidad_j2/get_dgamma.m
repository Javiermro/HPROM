function [dgamma] = get_dgamma (alpha_old,norm_s_trial,e_DatMatSet,e_VG)

% Variables globales
toldg = e_VG.toldg;

% Recupera propiedades del material
%E             = Eprop(4);
%poiss         = Eprop(5);
%sigmay        = Eprop(6);
%tita          = Eprop(8);
%hb            = Eprop(9);
%kinfb         = Eprop(10);
%kcerob        = Eprop(11);
%delta         = Eprop(12);
E = e_DatMatSet.young;
poiss = e_DatMatSet.poiss;
sigmay = e_DatMatSet.ftult;
tita = e_DatMatSet.tit;
hb = e_DatMatSet.hba;
kinfb = e_DatMatSet.kin;
kcerob = e_DatMatSet.kce;
delta = e_DatMatSet.del;

mu            = E/(2*(1+poiss));

% Inicializaciones
dgamma         = 0.0;
alpha_new      = alpha_old;
norm_gfunction = 1.0;

while (norm_gfunction > toldg)
    
    k_new          = sigmay + (tita*hb*alpha_new) + (kinfb-kcerob)*(1-exp(-delta*alpha_new));
    kp_new         = (tita*hb)+(kinfb-kcerob)*(delta*exp(-delta*alpha_new));
    gfunction      = -sqrt(2/3)*k_new + norm_s_trial - 2*mu*dgamma;
    dgfunction     = -2*mu*(1+kp_new/3/mu);
    dgamma         = dgamma - gfunction/dgfunction;
    alpha_new      = alpha_old + sqrt(2/3)*dgamma;
    norm_gfunction = abs(gfunction);

end

