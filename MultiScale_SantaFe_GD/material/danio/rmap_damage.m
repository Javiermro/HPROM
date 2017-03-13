function [ct,sigma_new,eps_new,hvar_new,aux_var] = rmap_damage(eps_new,hvar_old,e_DatMatSet,...
   e_VG)

%******************************************************************************************
%*  RETTURN-MAPPING PARA COMPUTO DE FUERZA INTERNA Y MATRIZ TANGENTE: PLANE STRAIN - 3D   *
%*  MODELO DE DAÑO ESCALAR CON ABLANDAMIENTO                                              *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

% Variables globales
factor          = 1e-6;

% Recupera propiedades del material
E = e_DatMatSet.young;
%sigma_u = e_DatMatSet.ftult;
H = e_DatMatSet.hba;
tita = e_DatMatSet.tit;
ce = e_DatMatSet.ce;
r_0 = e_DatMatSet.r_0;

% Recupera variables internas
r_old     = hvar_old(1);
q_old     = hvar_old(2);

% Inicializacion de variables
%r_0 = sigma_u/sqrt(E);
q_0 = r_0;
zero_q = factor*r_0;
if(r_old <= 0.0)
   %Se inicializa los valores iniciales de las variables de estado.
   r_old = r_0;
   q_old = r_0;
end

% Tension efectiva
effective_stress = ce*eps_new;

% Estado Trial
r_trial = sqrt(eps_new'*effective_stress);

if (r_trial > r_old && ~e_VG.elast)          
  
   % Daño
   fload_new = 1.0; 
   delta_r   = r_trial - r_old; 
   r_new     = r_trial;
   
   if (tita == 0)
      q_new     = q_old + H*delta_r;
   elseif (tita == 1)
      q_new     = q_0 * exp(H*(r_new-r_0)/q_0);
   end
   
   if(q_new < zero_q) 
      q_new = zero_q;
      H     = 0;
   end
   
   beta_1     = 1.0;
   beta_2     = q_new/r_new;
   if (tita == 0)
      beta_3  = (H*r_new-q_new)/r_new^3;
   elseif (tita == 1)
      beta_3  = -q_new/r_new^3 * (-H*r_new/q_0 + 1);
   end
   
   damage_new = 1.0 - (q_new/r_new);
   sigma_new  = beta_2*effective_stress;
    
else
    
   % Descarga elastica
   if (r_old > r_0)
      % Esto es para descarga del beta
      fload_new  = -1.0; 
   else
      fload_new  = 0.0; 
   end
   r_new      = r_old;
   q_new      = q_old;
   if E~=0
      damage_new = 1.0 - (q_new/r_new);
   else
      damage_new = 0;
   end
   sigma_new  =(1.d0 - damage_new) * effective_stress;
   beta_1     = 0.0;
   if E~=0
      beta_2     = q_new/r_new;
   else
      beta_2 = 1.0;
   end
   beta_3     = 0.0;

end

% Modulo elastoplastico consistente
ct = beta_2*ce + beta_1*beta_3*(effective_stress*effective_stress');

% Variables historicas
%hvar_new = [r_new ; q_new ; damage_new ; fload_new];
hvar_new = [r_new;q_new];

% Variables auxiliares
%aux_var = fload_new;
aux_var = [damage_new;fload_new];