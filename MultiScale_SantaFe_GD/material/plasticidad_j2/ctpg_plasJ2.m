function ct = ctpg_plasJ2(dgamma,norm_s_trial,kp_new,hp_new,N_new,mu,capa,e_VG)

%******************************************************************************************
%*  RETTURN-MAPPING PARA COMPUTO DEL TENSOR TANGENTE: PLANE STRAIN - 3D                   *
%*  MODELO DE PLASTICIDAD J2 CON ENDURECIMIENTO ISOTROPO                                  *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

%global SSOIT FOAT2
SSOIT = e_VG.SSOIT;
FOAT2 = e_VG.FOAT2;

tita_new   = 1 - (2*mu*dgamma)/(norm_s_trial);
tita_new_b = 1/(1+(kp_new+hp_new)/(3*mu)) - (1-tita_new);
ct         = (capa*SSOIT) + (2*mu*tita_new*FOAT2) - (2*mu*tita_new_b*(N_new*N_new.'));