function [kt,m_TensorTang] = elem_stiff_barra2D(coord_n,u,hvar_old,aux_var,Eprop,ce,...
   m_DefInicPaso,e_VG)

%******************************************************************************************
%*  EVALUACION DE LA MATRIZ TANGENTE PARA UN ELEMENTO DE BARRA 2D                         *
%*  SEGUN EL TIPO DE MODELO CONSTITUTIVO UTILIZADO                                        *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

%global_var
conshyp = e_VG.conshyp;
dofpe = e_VG.dofpe;
npg = e_VG.npg;
xg = e_VG.xg;
wg = e_VG.wg;
ntens = e_VG.ntens;
sihvarpg = e_VG.sihvarpg;
%siavarpg = e_VG.siavarpg;
%sitvare = e_VG.sitvare;

% Inicializaciones
% ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
kt = zeros(dofpe,dofpe);
eps_new = zeros(ntens,npg);
%Para la homogenización del tensor tangente
m_TensorTang = zeros(ntens,ntens,npg);

% Area elemento
% ¨¨¨¨¨¨¨¨¨¨¨¨¨
Area = Eprop(2);

hvar_old = reshape(hvar_old,sihvarpg,[]);

for iPG = 1:npg
   
   % Matriz de rotacion
   % ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
   dx = coord_n(1,2)-coord_n(1,1);
   dy = coord_n(2,2)-coord_n(2,1);
   tita = atan2(dy,dx);
   R = [cos(tita),-sin(tita);sin(tita),cos(tita)];
   Z = zeros(2,2);
   Q = [R,Z;Z,R];
   
   % Deformación
   % ¨¨¨¨¨¨¨¨¨¨¨
   %No se utiliza el punto de Gauss dentro de la función matrixB_barra2D
   [B,detJ] = matrixB_barra2D(coord_n,xg(iPG,:),e_VG);
   %m_DefInicPaso viene como tensor de 4 componentes (xx, yy, zz y xy) en notación Voigt.
   eps_new(:,iPG) = B*(Q'*u)+[cos(tita)^2,sin(tita)^2,0,sin(tita)*cos(tita)]*m_DefInicPaso;
   
   % Retturn-mapping segun el modelo constitutivo
   % ¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨¨
   switch conshyp
      case 1
         ct = rmapct_elas_barra2D(eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),Eprop,ce);
      case 2
         ct = rmapct_plasJ2_barra2D(eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),Eprop,ce,e_VG);
      case 3
         ct = rmapct_viscoplasJ2_barra2D(eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),Eprop,ce);
      case 4
         ct = rmapct_danio_barra2D(eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),Eprop,ce);
      case 5
         ct = rmapct_danio_stracc_barra2D(eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),Eprop,ce);
      case 6
         ct = rmapct_viscodanio_barra2D(eps_new(:,iPG),hvar_old(:,iPG),aux_var(:,iPG),Eprop,ce);
   end
   
   %Se almacena para la homogenización del tensor tangente
   m_TensorTang(:,:,iPG) = ct;
   
   % Cálculo de kt
   % ¨¨¨¨¨¨¨¨¨¨¨¨¨¨
   kt = kt+Q*(Area*B.'*ct*B*wg(iPG)*detJ)*Q';
   
end