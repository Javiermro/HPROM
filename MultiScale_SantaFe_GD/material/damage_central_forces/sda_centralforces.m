function [ct,sigma_new,hvar_new,aux_var] = sda_centralforces(eps_new,hvar_old,e_DatMatSet,B,e_VG)

%******************************************************************************************
%*  RETTURN-MAPPING PARA COMPUTO DE FUERZA INTERNA                                        *
%*  MODELO DE DAÑO ESCALAR DE FZAS CENTRALES (SANCHO-PLANAS)                              *
%*                                                                                        *
%*  ESQUEMA DE INTEGRACION:                                                               *
%*  IMPLI = 1 --> IMPLICITO                                                               *
%*  IMPLI = 0 --> IMPL-EX                                                                 *
%*                                                                                        *
%*  TIPO DE ABANDAMIENTO:                                                                 *
%*  MODFIS = 1 --> LINEAL                                                                 *
%*  MODFIS = 1 --> EXPONENCIAL                                                            *
%******************************************************************************************

% Esquema de integración
% impli = 1 Implicito
% impli = 0 Impl-Ex
impli = 0;

% Modelo de ablandamiento
% modFis = 0 ablandamiento lineal
% modFis = 1 ablandamiento exponencial
modFis = 1;

% Recupera propiedades del material
% young  = Eprop(4);
% %nu     = Eprop(5);
% ftult  = Eprop(6);
% gfval  = Eprop(7);
young = e_DatMatSet.young;
ftult = e_DatMatSet.ftult;
gfval = e_DatMatSet.gfv;
ce = e_DatMatSet.ce;


% Variables historicas del paso previo
%El vector salto del paso previo (no se utiliza, se guarda para imprensión)
%beta_pre = hvar_old(1:2);
%En xi se almacena el máximo valor histórico del módulo del salto.
xi_pre   = hvar_old(3);
qbifu    = hvar_old(4);
angle    = hvar_old(5);
first_t  = hvar_old(6);
n        = hvar_old(7:8);
%m_SigmaCrit = hvar_old(10:13);
%nsoli = hvar_old(10);

%Se fuerza que la normal sea horizontal
%n = [1;0];

% Tensor de elasticidad
dmatx = ce;

% Tension de referencia (tensión elástica trial, con el que se obtiene el vector tracc = n·E:Eps^a)
sigma_ref = dmatx*eps_new;

if (qbifu == 0)
   % La condicion de bifurcacion "NO" se ha verificado en pasos previos   

   %Tensión elástica en este paso
   sigma_new = sigma_ref;
   %Se obtiene el ángulo del plano de la tensión principal y las tensiones principales
   [angle,sigmacri] = cribis_sp(sigma_new);
   % Vector unitario normal a la fisura (actualzado a este paso, pero que no se utiliza en el mismo)
   n = [cos(angle);sin(angle)];
   %Para que grafique la normal en la zona elástica correctamente, se corrije para que apunte al 
   %nodo solitario (aunque para el cálculo no es necesario, solo cuando bifurca).
   %Si es por velocidad se puede comentar y descomentar la que está dentro del condicional siguiente
   n = get_solitary_node(n,B,e_VG);
   
   %Criterio de bifurcación, tensión principal máxima(Sigma 1) mayor a la tensión última
   if sigmacri>ftult
      % Busqueda del nodo solitario (activarlo en el caso no se haga en la parte elástica).
      %[n,n_m] = get_solitary_node(n,B,e_VG);

      % Vector traccion
      %qbifu toma el valor del módulo de vector tracción que se produce en el paso que se supera con
      %la tensión principal 1 a ftult (este valor puede ser un poco mayor a ftult). 
      qbifu = sigmacri;
      %first_t no se utiliza.
      %first_t = 1;
      
      %Abertura de fisura inicial
      %El valor xi inicial tiene que ser lo más chico posible pero que no haga que la matriz
      %tangente K del newton del modelo cohesivo no esté mal condicionada.
      %Con las siguientes operaciones se intenta buscar el valor más chico teniendo en cuenta que:
      %- Que en K hay que calcular 1/xi^3 o 1/xi^2, y por lo tanto que xi no sea tan chico que se 
      %supere el valor del doble máximo.
      %- Que el número de condición de K de no sea mayor a 1/eps(1) (notar que a menores valores de
      %xi la matriz K está peor condicionada). Por ello mediante algunas simplificaciones se
      %obtiene una aproximación a este mínimo xi.
      %- Que el valor de xi sea tal que la menor variación del vector tracción (eps(tu)) se detecte
      %como una mínima variación del xi (eps(xi)).
      xi = 10*eps(1)*qbifu/(cos((1-1e-3)*pi/2)*young);
      xi = max([xi;10*eps(xi)*qbifu/eps(qbifu);(qbifu/(realmax*0.1))^(1/3)]);
      %El vector beta no es variable histórica, y no se lo utiliza como tal en el cálculo (solo
      %se lo utiliza para ploteo de los resultados), por lo tanto se lo sigue manteniendo como
      %un vector nulo, aunque por la aproximación anterior ya no lo es más.
      %Ver si no poner beta = xi*n;
      beta = [0;0];
      
      %Se almacena el tensor crítico, donde bifurcó (previo a eso es cero, luego se mantiene
      %constante)
      %m_SigmaCrit = sigma_new;
   else
      beta = [0;0];
      xi = 0;
   end
   fload = 0;
   %Módulo tangente
   ct = dmatx;

else
   % La condicion de bifurcacion "SI" se ha verificado en pasos previos   

   % Busqueda del nodo solitario y cálculo de grad_phi
   %(también se modifica el sentido de la normal n, para que apunte al nodo solitario).
   [n,grad_phi] = get_solitary_node(n,B,e_VG);
   n_m = f_Vec2Voigt2D(n,e_VG);
   grad_phi = f_Vec2Voigt2D(grad_phi,e_VG);
   
   % Vector traccion
   tracc = n_m'*sigma_ref;
   
   Eb = dmatx*grad_phi;
   nEb = n_m'*Eb;
   
   % Integracion explicita a partir de la solucion implícita previa.
   % Tambien puede utilizarse como estado trial para la integracion implícita.
   %(sirve como condición inicial del newton de la solución implícita y como valor en caso de 
   %descarga al avanzar con xi_pre fijo).
   switch modFis
      case 0
         factor_expli = (1/xi_pre-qbifu/2/gfval)*qbifu;
      case 1
         factor_expli = qbifu/xi_pre*exp(-qbifu*xi_pre/gfval);
      otherwise
         error('Modelo de fuerzas centradas: Modelo de fisura no definido.')
   end

   qfi_expli = factor_expli*eye(2,2)+nEb;
   %Se resuelve el sistema de equilibrio local en forma explícita (las variables explicitas
   %son las históricas, es decir xi_pre y la normal n, las otras provienen directamente de la 
   %deformación del paso actual a resolver). Y tracc sale del producto n·E:Eps, con la normal
   %n del paso previo.
   beta_expli = qfi_expli\tracc;
   mod_beta_expli = norm(beta_expli);
   
   % Integración implícita   
   if mod_beta_expli>xi_pre   %beta_expli'*n>=0 
   % Carga en régimen de discontinuidad fuerte
      if impli==1         
         switch modFis
            case 0
               [beta_impli,mod_beta_impli,factor_impli,der_chi_impli] = f_newton_cohesive_ds_lineal(...
                  qbifu,mod_beta_expli,gfval,nEb,tracc,beta_expli);
            case 1
               [beta_impli,mod_beta_impli,factor_impli,der_chi_impli] = newton_cohesive_ds(...
                  qbifu,mod_beta_expli,gfval,nEb,tracc,beta_expli);
            otherwise
               error('Modelo de fuerzas centradas: Modelo de fisura no definido.')
         end
            
         % Se verifica nuevamente la condición de carga
         if mod_beta_impli<=xi_pre  %beta_impli'*n<0
         % IMPLÍCITO: La condición de carga "NO" se verifica
            fload = -1;
            % Aca creo que hay que controlar que beta no cambie de sentido
            xi = xi_pre;
            beta = beta_expli;            
            %factor = factor_expli;
            qfi = qfi_expli;
            der_chi = 0;
         else
         % IMPLÍCITO: La condicion de carga "SÍ" se verifica
            fload = 1;
            %Se verifica que mod_beta_impli>xi_pre
            xi = mod_beta_impli;           
            beta = beta_impli;            
            factor = factor_impli;            
            qfi = factor*eye(2,2)+nEb;
            der_chi = der_chi_impli;
         end         
      else 
      % IMPLEX: La condicion de carga "SÍ" se verifica.
         fload = 1;
         %En este caso se verifica que mod_beta_expli>xi_pre
         xi = mod_beta_expli;
         beta = beta_expli;
         %factor = factor_expli;
         qfi = qfi_expli;
         der_chi = 0;
      end
   else
   % IMPLÍCITO Y IMPLEX: Descarga en régimen de discontinuidad fuerte
      fload = -1;
      % Aca creo que hay que controlar que beta no cambie de sentido
      xi = xi_pre;
      %Aún en descarga el vector salto varía (no es una variable histórica), lo que no varía es el
      %xi (módulo de salto histórico), es decir el vector puede girar dentro de un círculo de radio 
      %xi.
      beta = beta_expli;      
      %factor = factor_expli;
      qfi = qfi_expli;
      der_chi = 0;
   end
   
   %Se modifica la deformación considerando el valor de beta
   eps_new = eps_new-grad_phi*beta;
   % Evaluacion del tensor de tensiones
   sigma_new = dmatx*eps_new;
   %der_chi es nula en el caso de ser explicito o estar en descarga.
   qfi = qfi+der_chi*(beta*beta');
   nEbq = Eb/qfi*n_m';
   % Tensor tangente
   ct = dmatx-nEbq*dmatx;
   
   %Se indica que ya no es no es el primer paso después de la bifurcación
   %first_t = 0;   
end

% Limite para el valor de traccion
if (modFis == 0)
   % si el modelo es con ablandamiento lineal hay que restringir el vector
   % de tracciones
end

% Actualizacion de variables internas
hvar_new = zeros(e_DatMatSet.sihvarpg,1);
hvar_new(1:2) = beta;
hvar_new(3) = xi;  
hvar_new(4) = qbifu;
hvar_new(5) = angle;
hvar_new(6) = first_t;
hvar_new(7:8) = n;
hvar_new(9) = fload;
%hvar_new(10:13) = m_SigmaCrit;
%hvar_new(10) = nsoli;

% Almancenamiento de variables auxiliares
aux_var(1) = 0;

end