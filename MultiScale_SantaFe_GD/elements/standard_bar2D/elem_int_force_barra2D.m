function [fint,sigma_new,eps_new,hvar_new,aux_var,sigma_nodo] = elem_int_force_barra2D(coord_n,...
    u,hvar_old,Eprop,ce,m_DefInicPaso,e_VG)

%******************************************************************************************
%*  EVALUACION DE LA FUERZA INTERNA PARA UN ELEMENTO DE BARRA 2D                          *
%*  SEGUN EL TIPO DE MODELO CONSTITUTIVO UTILIZADO                                        *                  
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

%global_var
conshyp = e_VG.conshyp;
npe = e_VG.npe;
dofpe = e_VG.dofpe;
npg = e_VG.npg;
xg = e_VG.xg;
wg = e_VG.wg;
ntens = e_VG.ntens;
%sihvare = e_VG.sihvare;
sihvarpg = e_VG.sihvarpg;
%siavare = e_VG.siavare;
siavarpg = e_VG.siavarpg;
%sitvare = e_VG.sitvare;

% Inicializaciones
% จจจจจจจจจจจจจจจจ
fint = zeros(dofpe,1);
sigma_new = zeros(ntens,npg);
eps_new = zeros(ntens,npg);
hvar_new = zeros(sihvarpg,npg);
aux_var = zeros(siavarpg,npg);
sigma_nodo = zeros(ntens,npe);

% Area elemento
% จจจจจจจจจจจจจ
Area = Eprop(2);

hvar_old = reshape(hvar_old,sihvarpg,[]);

for iPG = 1:npg;
   
   % Matriz de rotacion
   % จจจจจจจจจจจจจจจจจจ
   dx = coord_n(1,2)-coord_n(1,1);
   dy = coord_n(2,2)-coord_n(2,1);
   tita = atan2(dy,dx);
   R = [cos(tita),-sin(tita);sin(tita),cos(tita)];
   Z = zeros(2,2);
   Q = [R,Z;Z,R];

   % Deformaci๓n
   % จจจจจจจจจจจ
   %No se utiliza el punto de Gauss dentro de la funci๓n matrixB_barra2D
   [B,detJ] = matrixB_barra2D(coord_n,xg(iPG,:),e_VG);
   %m_DefInicPaso viene como tensor de 4 componentes (xx, yy, zz y xy) en notaci๓n Voigt.
   eps_new(:,iPG) = B*(Q'*u)+[cos(tita)^2,sin(tita)^2,0,sin(tita)*cos(tita)]*m_DefInicPaso;

   % Retturn-mapping segun el modelo constitutivo
   % จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
   switch conshyp
      case 1
         [sigma_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmapfi_elas_barra2D(...
            eps_new(:,iPG),hvar_old(:,iPG),Eprop,ce);
      case 2
         [sigma_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmapfi_plasJ2_barra2D(...
            eps_new(:,iPG),hvar_old(:,iPG),Eprop,ce,e_VG);
      case 3
         %Esta funci๓n no estแ.
         [sigma_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmapfi_viscoplasJ2_barra2D(...
            eps_new(:,iPG),hvar_old(:,iPG),Eprop,ce);
      case 4
         [sigma_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmapfi_danio_barra2D(...
            eps_new(:,iPG),hvar_old(:,iPG),Eprop,ce);
      case 5
         [sigma_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmapfi_danio_stracc_barra2D(...
            eps_new(:,iPG),hvar_old(:,iPG),Eprop,ce);
      case 6
         [sigma_new(:,iPG),hvar_new(:,iPG),aux_var(:,iPG)] = rmapfi_viscodanio_barra2D(...
            eps_new(:,iPG),hvar_old(:,iPG),Eprop,ce);
   end

   % Cแlculo de fint
   % จจจจจจจจจจจจจจจ
   fint = fint+Q*(Area*B.'*sigma_new(:,iPG)*wg(iPG)*detJ);
   
end

%Se ordena las matrices como vectores columnas
sigma_new = sigma_new(:);
eps_new = eps_new(:);
hvar_new = hvar_new(:);
aux_var = aux_var(:);
sigma_nodo = sigma_nodo(:);