function [i_bifurc,bif_angle,leq_element_MACRO,n_bifurc,n_tensor] = ...
   bifurcation_condition_MACRO(ct,...
   leq_elem_old,i_bifurc_old,e_VG,xx_MACRO,...
   conec_MACRO,u_new,e_DatElemSet,dN_xy,e_DatSet_iSet)

%m_NumElem = e_DatSet_iSet.m_NumElem;

%[~,DD]=eig(ct);
%VAL=min(diag(DD));
% if VAL<0.0
% fprintf('ELEMENT: %i - EIGENVALUE: %f \n',m_NumElem(iElem_MACRO),VAL)
% end

% ************************************
% RAYLEIGH QUOTIENT METHOD ***********a
% ************************************

% Inicializaciones de las variables del procedimiento
% ***************************************************
sym = true;
%p=[ 0; 1];
%Se elige un ángulo arbitrario como valor inicial.
ang = 0.2345432;
p = [cos(ang);sin(ang)];
n_i1 = [p(1),0;0,p(2);0,0;p(2),p(1)];
Qpp = n_i1'*ct*n_i1;
if sym
   Qpp = (Qpp+Qpp')/2;
end
[V,lambda] = eig(Qpp);
%a = diag(lambda);
a = lambda([1;4]);
lambda_zero = max(a);
[lambda_i1,x1] = min(a);

%La función eig devuelve los autovectores normalizados, así que no es necesario hacerlo explícitamente.
N_i1 = V(:,x1);
lambda_i = 100*lambda_i1;

iter = 0;
tol = 1.e-8;

% Algoritmo Iterativo
% *******************
%Al usar iter==0 ó lambda_i = 100*lambda_i1; (esta última es dependiente de la tolerancia y el factor que
%multiplica, pero es muy difícil que no lo sea) se asegura que al menos entre una vez.
while abs((lambda_i1-lambda_i)/lambda_zero)>tol&&iter<30||iter==0
   iter=iter+1;
   lambda_i=lambda_i1;
   %  detq_i=detq_i1;
   %  N_i2=N_i;
   N_i = N_i1;
   n_i1 = [N_i1(1),0;0,N_i1(2);0,0;N_i1(2),N_i1(1)];
   Qpp = n_i1'*ct*n_i1;
   
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% simetrizada  %%%%%%%%%%%%%%%%
   %Por hipótesis del esquema de solución numérico se simetriza el tensor Q.
   if sym
      Qpp = (Qpp+Qpp')/2;
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   [V,lambda] = eig(Qpp);
   [lambda_i1,x1] = min(lambda([1;4]));

   %La función eig devuelve los autovectores normalizados, así que no es necesario hacerlo explícitamente.
   N_i1 = V(:,x1);
end

% if iter>28
%    warning('off','backtrace')
%    warning('Análisis de bifurcación: Elemento %d: bifurcation analysis does not converge.',e_VG.iElemNum) %#ok<WNTAG>
%    warning('on','backtrace')
% end

if norm(N_i1-N_i)<1e-5
   fprintf('Análisis de bifurcación: Elemento %d: Se obtuvo el mismo ángulo de bifurcación.\n',e_VG.iElemNum)
   ct %#ok<NOPRT>
end

eltype_MACRO = e_DatElemSet.eltype;
%delta_u       = u_new(dof)-u_old(dof);
coord_n_MACRO = f_CoordElem(xx_MACRO,conec_MACRO);

% Calculo vector gradiente ||u|| normalizado
switch eltype_MACRO
   case {4,10,21,22,23} % Elementos cuadrilaterales - indistinto de si es enriquecido o no.
      %grad_norm_u_vector = matrix_grad_u_quad_q1_epd(coord_n_MACRO,delta_u);
      grad_norm_u_vector = matrix_grad_u_quad_q1_epd(u_new,dN_xy);
   otherwise
      error('This case is not implemented yet!');
end

% DISPLACEMENT NORMAL OF THE ELEMENT
n_tensor = zeros(e_VG.ntens,e_VG.ndime);

if lambda_i1<0.0  % CASO BIFURCACION ##############
   
   i_bifurc = 1; % FLAG DE BIFURCACION
   
   n_bifurc(1,1) =  N_i1(1);
   n_bifurc(2,1) =  N_i1(2);
   n_bifurc(1,2) =  N_i(1);
   n_bifurc(2,2) =  N_i(2);
   
   if abs(grad_norm_u_vector'*N_i1)>abs(grad_norm_u_vector'*N_i)
      bif_angle = atan(N_i1(2)/N_i1(1))*180/pi;
   else
      bif_angle = atan(N_i(2)/N_i(1))*180/pi;
   end
   
   % Calculo longitud de Oliver (Oliver X./89)
   npe = e_DatElemSet.npe ;
   %  npe = e_VG.npe;
   new_coord_n_MACRO = zeros(size(coord_n_MACRO,1),size(coord_n_MACRO,2));
   xcg = sum(coord_n_MACRO(1,:))/npe;
   ycg = sum(coord_n_MACRO(2,:))/npe;
   bif_angle_rad=bif_angle*pi/180;
   
   for inode = 1:npe
      new_coord_n_MACRO(1,inode) = (coord_n_MACRO(1,inode)-xcg)*cos(bif_angle_rad)+(coord_n_MACRO(2,inode)-ycg)*sin(bif_angle_rad);
      new_coord_n_MACRO(2,inode) = -(coord_n_MACRO(1,inode)-xcg)*sin(bif_angle_rad)+(coord_n_MACRO(2,inode)-ycg)*cos(bif_angle_rad);
   end
   
   % Calculo longitud equivalente para los elementos macro (Oliver/89)
   switch eltype_MACRO
      case {4,10,21,22,23} % Elementos cuadrilaterales - indistinto de si es enriquecido o no.
         leq_element_MACRO = le_quad_q1_epd(coord_n_MACRO,new_coord_n_MACRO,bif_angle_rad);
         %[leq_element_MACRO,fii] = le_quad_q1_epd(coord_n_MACRO,new_coord_n_MACRO,bif_angle_rad);
      otherwise
         error('This case is not implemented yet!');
   end
   
else % CASO DE NO BIFURCACION ##############
   
   i_bifurc = i_bifurc_old;
   leq_element_MACRO = leq_elem_old;
   bif_angle = 0;
   
   n_bifurc(1,1) =  N_i1(1);
   n_bifurc(2,1) =  N_i1(2);
   n_bifurc(1,2) =  N_i(1);
   n_bifurc(2,2) =  N_i(2);
   
end

% por si salta de no bifurcar a haber inyectado... caso que no deberia pasar
% NORMAL TO THE LOCALIZATION BAND
%grad_norm_u_vector = [-0.866025403784439 0.5]';
if(abs(grad_norm_u_vector'*N_i1)> abs(grad_norm_u_vector'*N_i))
   n_tensor(1,1) = N_i1(1);
   n_tensor(2,2) = N_i1(2);
   n_tensor(4,1) = N_i1(2);
   n_tensor(4,2) = N_i1(1);
else
   n_tensor(1,1) = N_i(1);
   n_tensor(2,2) = N_i(2);
   n_tensor(4,1) = N_i(2);
   n_tensor(4,2) = N_i(1);
end

end