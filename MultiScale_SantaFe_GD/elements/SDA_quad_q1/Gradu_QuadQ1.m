function [n_vector,Grad_u] = Gradu_QuadQ1(coord_n,u)

%*********************************************************************************
%*  function [B,detJ] = bmatrix_CL4N_EP (coord_n,xg,npe)                         *
%*  Obtiene la matriz B del cuadr�ngulo lineal de 4 nodos para EPD � EPT         *
%*  evaluada en xg.                                                              *
%*                                                                               *
%*  ARGUMENTOS DE ENTRADA:                                                       *
%*  ���������������������                                                        *
%*  coord_n : Coordenadas de los nodos que componen el elemento.                 *
%*            coord_n = [x1  x2  x3 x4                                           *
%*                       y1  y2  y3 y4]                                          *
%*  numeraci�n de nodos: 4-------------3                                         *
%*                       |             |                                         *
%*                       |             |                                         *
%*                       |             |                                         *
%*                       1-------------2                                         *
%*  xg      : Vector de posici�n del punto de Gauss en estudio.                  *
%*            xg = [E n]                                                         *
%*  npe     : N�mero de nodos por elemento.                                      *
%*                                                                               *
%*  LISTA DE VARIABLES:                                                          *
%*  ������������������                                                           *
%*  E       : coordenada natural                                                 *
%*  n       : coordenada natural                                                 *
%*  dN_E    : Derivadas de las funciones de forma respecto a "E"                 *
%*            dN_E = [dN1_E  dN2_E  dN3_E dN4_E]                                 *
%*  dN_n    : Derivadas de las funciones de forma respecto a "n"                 *
%*            dN_n = [dN1_n  dN2_n  dN3_n dN4_n]                                 *
%*  dN_En   : Derivadas de las funciones de forma respecto a "E" y a "n"         *
%*  J       : matriz jacobiana de la transformaci�n                              *
%*            J = [dx_E dy_E ; dx_n dy_n]                                        *
%*  dN_xy   : Derivadas de las funciones de forma respecto a "x" y a "y"         *
%*                                                                               *
%*  ARGUMENTOS DE SALIDA:                                                        *
%*  ��������������������                                                         *
%*  B       : Matriz B para elemento cuadrangular de 4 nodos en estado plano.    *
%*  detJ    : determinante de la matriz jacobiana.                               *
%*********************************************************************************

E = 0;
n = 0;
% Derivadas de Ni respecto "E" y "n":
dN_E = [(1/4*n-1/4)  (-1/4*n+1/4)  (1/4*n+1/4)  (-1/4*n-1/4)];
dN_n = [(1/4*E-1/4)  (-1/4*E-1/4)  (1/4*E+1/4)  (-1/4*E+1/4)];
dN_En = [dN_E ; dN_n];
% Matriz Jacobiana:
J11 = coord_n(1,:)*dN_E';
J12 = coord_n(2,:)*dN_E';
J21 = coord_n(1,:)*dN_n';
J22 = coord_n(2,:)*dN_n';
J = [J11 J12 ; J21 J22];
% Matriz dN_xy:
% ��������
dN_xy = J\dN_En;

Grad_u = dN_xy*[u(1:2:end) u(2:2:end)];

flag =2;

switch flag
    
    case 1
        % Funciones de forma respecto "E" y "n":
        N_E = [(1/4)*(1-E)*(1-n) (1/4)*(1+E)*(1-n) (1/4)*(1+E)*(1+n) (1/4)*(1-E)*(1+n)];
        N_n = N_E;
        
        u_mid(1) = N_E*u(1:2:end);
        u_mid(2) = N_n*u(2:2:end);
        
        if (norm(u_mid)==0)
            norm_despl=1;
        else
            norm_despl=norm(u_mid);
        end
        
        % estimacion de Grad(||u||)
        Grad_norm_u = (1/norm_despl)*(Grad_u*u_mid');
        
        if (norm(Grad_norm_u)==0)
            norm_Grad=1;
        else
            norm_Grad=norm(Grad_norm_u);
        end
        
        % estimacion de Grad(||u||)/norm(Grad(||u||))
        n_vector = Grad_norm_u/norm_Grad;
        
    case 2
        
        norm_grad_x = norm(Grad_u(:,1));
        norm_grad_y = norm(Grad_u(:,2));
        
        if norm_grad_x > norm_grad_y
            n_vector = Grad_u(:,1)/norm_grad_x;
        else
            n_vector = Grad_u(:,2)/norm_grad_y;
        end
        
end