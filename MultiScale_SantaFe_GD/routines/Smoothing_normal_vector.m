function n_smoothing = Smoothing_normal_vector(xx,conec,u,dN_xy,eltype,e_VG)

%coord_n_MACRO = f_CoordElem(xx,conec);
% Calculo vector gradiente ||u|| normalizado
switch eltype
    case {4,10,11,21,22,23} % Elementos cuadrilaterales - indistinto de si es enriquecido o no.
        %grad_norm_u_vector = matrix_grad_u_quad_q1_epd(coord_n_MACRO,delta_u);
        grad_norm_u_vector = matrix_grad_u_quad_q1_epd(u,dN_xy);
    otherwise
        error('This case is not implemented yet!');
end

% DISPLACEMENT GRADIENT NORM

%Normal de referencia
%Sebastian: es la que se ingresa en el archivo de datos en N_REF_SMOOTHING, y pareciera que solo sirve para
%darle sentido a n_smoothing.
nRefSmoothing = e_VG.nRefSmoothing;

prod = nRefSmoothing'*grad_norm_u_vector;

%if (grad_norm_u_vector(1) < 0) || ((grad_norm_u_vector(1) == 0) && (grad_norm_u_vector(2) < 0))
if prod < 0 
    %n_smoothing = -n_tensor;
    n_smoothing(1,1) = -grad_norm_u_vector(1);
    n_smoothing(2,2) = -grad_norm_u_vector(2);
    n_smoothing(4,1) = -grad_norm_u_vector(2);
    n_smoothing(4,2) = -grad_norm_u_vector(1);
else
    %n_smoothing = n_tensor;
    n_smoothing(1,1) = grad_norm_u_vector(1);
    n_smoothing(2,2) = grad_norm_u_vector(2);
    n_smoothing(4,1) = grad_norm_u_vector(2);
    n_smoothing(4,2) = grad_norm_u_vector(1);
end

%     n_smoothing(1,1) = 0;
%     n_smoothing(2,2) = -1;
%     n_smoothing(4,1) = -1;
%     n_smoothing(4,2) = 0;
