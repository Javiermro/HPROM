function [ct,sigma_n1,sigma_n1_impl,hvar_n1,aux_var] = ...
   rmap_damage_IMPLEX (eps_n1,hvar_n,e_DatMatSet,Eprop_qSD,e_VG)

%*********************************************************************************
%*                     CURSO DE MODELOS CONSTITUTIVOS - UPC -                    *
%*           Algoritmo de integraciï¿½n para el modelo de daï¿½o isï¿½tropo            *
%*                                                                               *
%* Llamada :: [sigma_n1,hvar_n1,aux_var] = rmap_dano (eps_n1,hvar_n,Eprop,ce)    *
%*                                                                               *
%* ARGUMENTOS ENTRADA: eps_n1(4)   deformaciï¿½n total , paso n+1                  *
%*                                 vector R4    (exx eyy exy ezz)                *
%*                     hvar_n(5)   variables internas , paso n                   *
%*                                 hvar_n(1:4) deformacion total paso n          *
%*                                 hvar_n(5) = r  ; hvar_n(6)=q                  *
%*                     Eprop(4)    vector de propiedades de material             *
%*                                 Eprop(1)=  E     : modulo de Young            *
%*                                 Eprop(2)=  nu    : modulo de Poisson          *
%*                                 Eprop(3)=  H     : modulo de Softening/hard.  *
%*                                 Eprop(4)=sigma_u : tensiï¿½n ï¿½ltima             *
%*                     ce(4,4)     tensor constitutivo elï¿½stico                  *
%*                                                                               *
%* ARGUMENTOS SALIDA:  sigma_n1(4) tensiï¿½n de cauchy , paso n+1                  *
%*                     hvar_n(5)   variables internas , paso n+1                 *
%*                     aux_var(3)  variables auxiliares para calculo tensor tang.*
%*********************************************************************************
%global sihvarpg

% Variables - Propiedades mecï¿½nicas
% *********************************
E        = e_DatMatSet.young;
nu       = e_DatMatSet.poiss;
%sigma_u  = e_DatMatSet.ftult;
tita     = e_DatMatSet.tit;
Gf       = e_DatMatSet.gfv;
ce       = e_DatMatSet.ce;
r0       = e_DatMatSet.r_0;
ksd      = Eprop_qSD.ksd;
esImplex = e_DatMatSet.esImplex;

% Parametros de localizacion - energia de fractura
% ************************************************
%ago   bandwidth = e_VG.bandwidth;
% Dtime     = e_VG.Dtime;

% Recuperacion de variables internas del paso anterior
% ****************************************************
%Estas variables se asumen inicializadas con r0 (ver función f_eVarEstInic).
r_n = hvar_n(5);
q_n = hvar_n(6);
% NOTAR: ASUMIMOS QUE LA VARIABLE DTIME ES CONSTANTE !!!!!
%%%%  Dtime_n       = hvar_n(12); if Dtime_n==0; Dtime_n=Dtime; end
%alf Delta_r_n     = hvar_n(13);
%alf bandwidth_pre = hvar_n(14); if bandwidth_pre==0; bandwidth_pre=1; end
bandwidth_pre = hvar_n(2);
Delta_r_n     = hvar_n(3);   %  se cambio por hvar_n(13)
%if bandwidth_pre==0
%   bandwidth_pre=1;
%end  %  se cambio por hvar_n(14)

% Inicializacion de las variables
% *******************************
hvar_n1 = zeros(length(hvar_n),1);
%q0 = r0;
zero_q = 1.d-5*r0;
%Se inicializa hasta que comienza a dañar.
% if r_n<=0
%    r_n = r0;
%    q_n = r0;
% end

% Determinacion del parametro de ablandamiento H regularizado
% ***********************************************************
bandwidth = ksd; %if bandwidth==0; bandwidth=1; end
%Sebastian: La regularización depende del tipo de ablandamiento.
if tita==0 %Lineal
   %Hbar = -(sigma_u^2)/(2*E*Gf);
   Hbar = -r0^2/Gf/2;
else  %Exponencial
   %Hbar = -(sigma_u^2)/(E*Gf);
   %En el caso de ablandamiento exponencial y donde se produzcan cambios del espesor de regularización 
   %en el tiempo, es necesario cambiar el modo que se regulariza para poder capturar correctamente la energía
   %fractura.
   %Se asume que la variable hvar(7) viene inicializada con q0 (ver función f_eVarEstInic).
   %Solo se cambia la qInic cuando se produce un cambio del espesor de regularización.
   if abs(bandwidth-bandwidth_pre)>0
%       fprintf('Cambió la longitud de regularización en el elemento %d de %f a %f.\n',e_VG.iElemNum,...
%          bandwidth_pre,bandwidth)
      rInic = r_n;
      qInic = q_n;
   else
      rInic = hvar_n(end-1);
      qInic = hvar_n(end);
   end
   hvar_n1(end-1) = rInic;
   hvar_n1(end) = qInic;
   %rInic = r0;
   %qInic = r0;
   Hbar = -qInic*r0/Gf;
end
H = ksd*Hbar;

%Tension elastica
tension_efectiva = ce*eps_n1;

%----------------   IMPLEX STRESS   --------------------
% ******************************************************
if esImplex  % IMPLICIT stresses -------------
   
   Delta_r_tilde = (bandwidth_pre/bandwidth)*Delta_r_n;
   r_tilde = r_n+Delta_r_tilde;
   if (tita == 0) % Comportamiento lineal
      q_tilde = q_n+H*Delta_r_tilde;
   elseif (tita == 1) % Comportamiento exponencial
      %q_tilde = q0 * exp(H*(r_tilde-r0)/q0);
      %q_tilde = q_n+H*exp(H*(r_tilde-r0)/q0)*Delta_r_tilde;
      %El H cambia en el tiempo, por lo que la evolución del daño debería ser incremental. En el caso
      %explícito se ve que oscila menos si se hace en forma incremental, pero el implícito es mejor si se hace
      %en totales.
      q_tilde = q_n+H*exp(H*(r_tilde-rInic)/qInic)*Delta_r_tilde;
      %q_tilde = qInic*exp(H*(r_tilde-rInic)/qInic);
   end
   % Limite inferior para ablandamiento
   if q_tilde<zero_q
%       fprintf('Se alcanzó el límite inferior de daño en el elemento %d para el IMPLEX (q = %g).\n',...
%          e_VG.iElemNum,q_tilde)
      q_tilde = zero_q;
   end
   % Tension tilde (aproximada)
   beta2tilde = q_tilde/r_tilde;
   sigma_n1 = beta2tilde*tension_efectiva;
   c_tilde  = beta2tilde*ce ;
   
end
% ******************************************************
%--------------- END IMPLEX STRESS ---------------------
% tension efectiva
MDtype = 2; % MDtype =  1 : SYMMETRIC  % MDtype =  2 : ONLY TENSION  % MDtype =  3 : NON-SYMMETRIC
n = 0;
[r_trial,tension_efectiva2,theta] = rtrial_damage(MDtype,ce,eps_n1,E,nu,n);

% Return mapping del modelo de danio
if(r_trial > r_n)
   %Estado de carga
   fLoad = 1;
   % ############## ESTADO INELASTICO (O DE CARGA INELASTICA)
   r_n1 = r_trial;
   
   if (tita == 0) % Comportamiento lineal
      cal_H = H;
      delta_r = r_n1-r_n;
      q_n1 = q_n+H*delta_r;
   elseif (tita == 1) % Comportamiento exponencial
      %cal_H = H*exp(H*(r_n1-r0)/q0);
      %q_n1 = q0*exp(H*(r_n1-r0)/q0);
      %q_n1 = q_n+cal_H*delta_r;
      expFq = exp(H*(r_n1-rInic)/qInic);
      cal_H = H*expFq;
      %Se vio que en totales con el implícito se obtiene una respuesta que oscila menos y la curva se acerca
      %más a cero (degrada completamente).
      q_n1 = qInic*expFq;
      %delta_r = r_n1-r_n;
      %q_n1 = q_n+cal_H*delta_r;
   end
  
   % Limite inferior para ablandamiento
   if(q_n1<zero_q)
      fLoad = 2;
      %fprintf('Se alcanzó el límite inferior de daño en el elemento %d para el IMPLÍCITO (q = %g).\n',...
      %   e_VG.iElemNum,q_n1)
      q_n1 = zero_q;
      %H = 0;
      cal_H = 0;
   end
   % Parametros para el calculo del tensor constitutivo tangente
   %beta_1 = 1.0;
   beta_2 = q_n1/r_n1;
   dano_n1 = 1-beta_2;
   
   %    if tita==0 % Comportamiento lineal
   %       beta_3 = (H*r_n1-q_n1)/r_n1^3;
   %    elseif tita==1 % Comportamiento exponencial
   %       beta_3 = -(q_n1/r_n1^3)*(-cal_H*beta_2+1);
   %    end
   beta_3 = (cal_H*r_n1-q_n1)/r_n1^3;
   
   sigma_n1_impl = beta_2*tension_efectiva;
   %
   if MDtype==1||MDtype==2
      ct_implicit = beta_2*ce+beta_3*(tension_efectiva*tension_efectiva2');
   elseif MDtype==3
      ct_implicit = beta_2*ce+beta_3*((r_trial^2/theta)*(tension_efectiva*tension_efectiva2')+...
         theta^2*(tension_efectiva*tension_efectiva'));
   end
   
else
   % ############# ESTADO ELASTICO (O DE DESCARGA ELASTICA)
   r_n1     = r_n;
   q_n1     = q_n;
   %beta_1   = 0.0;
   beta_2   = q_n1/r_n1;
   %beta_3   = 0.0;
   dano_n1 = 1-beta_2;
   %Estados de carga
   if dano_n1==0
      %Estado de carga o descarga elástica sin haber nunca dañado
      fLoad = 0;
   else
      %Estado de carga o descarga elástica después de haberse producirdo daño.
      fLoad = -1;
   end
   %
   sigma_n1_impl = beta_2*tension_efectiva;
   ct_implicit = beta_2*ce;
end

% Matriz tangente de danio consistente
% if MDtype==1||MDtype==2
%    ct_implicit = beta_2*ce+beta_1*beta_3*(tension_efectiva*tension_efectiva2');
% elseif (MDtype == 3)
%    ct_implicit = beta_2*ce+beta_1*beta_3*((r_trial^2/theta)*(tension_efectiva*tension_efectiva2') + theta^2*(tension_efectiva*tension_efectiva'));
% end

% Actualizaciï¿½n de variables historicas
% *************************************
hvar_n1(1)    = dano_n1;
hvar_n1(2)    = bandwidth;     % bandwidth de la banda de localizacion macro  %  se cambio por hvar_n(14)
hvar_n1(3)    = r_n1-r_n;    % Delta_r_n  %  se cambio por hvar_n(13)
%Variable de suavizado del crack path field.
hvar_n1(4)    = hvar_n(4)+ksd*(r_n1-r_n) ;              % ro en punto de gauss
hvar_n1(5)    = r_n1;
hvar_n1(6)    = q_n1;
hvar_n1(7)    = r_n1-r0;     % Delta_r_n  %  se cambio por hvar_n(13)
hvar_n1(9) = r_trial;

if ~esImplex  % IMPLICIT stresses -------------
   %Criterio de carga y descarga.
   hvar_n1(8) = r_n1-r_n;
   % Tension tilde (aproximada)
   sigma_n1 = sigma_n1_impl;
   ct = ct_implicit;
else
   %Criterio de carga y descarga.
   hvar_n1(8) = Delta_r_tilde;
   %hvar_n1(8) = r_trial-r_n;
   ct = struct('Implex',c_tilde,'Impli',ct_implicit);
end

% Variables auxiliares
aux_var = fLoad;


end