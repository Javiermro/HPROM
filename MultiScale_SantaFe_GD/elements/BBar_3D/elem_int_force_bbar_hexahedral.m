function [fint,sigma_new,eps_new,hvar_new,aux_var,sigma_nodo] = elem_int_force_bbar_hexahedral(...
    coord_n,u,hvar_old,Eprop,ce,e_VG)
      
%global_var
struhyp = e_VG.struhyp;
conshyp = e_VG.conshyp;
npe = e_VG.npe;
dofpe = e_VG.dofpe;
npg = e_VG.npg;
xg = e_VG.xg;
wg = e_VG.wg;
ntens = e_VG.ntens;
sihvare = e_VG.sihvare;
sihvarpg = e_VG.sihvarpg;
siavare = e_VG.siavare;
siavarpg = e_VG.siavarpg;
sitvare = e_VG.sitvare;

% Inicializaciones
% จจจจจจจจจจจจจจจจ
fint       = zeros(dofpe,1);
sigma_new  = zeros(sitvare,1);
eps_new    = zeros(sitvare,1);
hvar_new   = zeros(sihvare,1);
aux_var    = zeros(siavare,1);
sigma_nodo = zeros(ntens,npe);

% Cแlculo de la matriz Bh
% จจจจจจจจจจจจจจจจจจจจจจจ

BH = matrixBh_bbar_hexahedral(coord_n,e_VG);

for i = 1:npg;
    for j = 1:npg;
        for k = 1:npg;
            
            % Cแlculo de ํndices
            % จจจจจจจจจจจจจจจจจจ
            ini =    ntens*(npg*i-npg+j*(j-1)+i*(i-1)+k-1)+1 ; inf =    ntens*(npg*i-npg+j*(j-1)+i*(i-1)+k) ; itv = ini:inf;
            ini = sihvarpg*(npg*i-npg+j*(j-1)+i*(i-1)+k-1)+1 ; inf = sihvarpg*(npg*i-npg+j*(j-1)+i*(i-1)+k) ; ihv = ini:inf;
            ini = siavarpg*(npg*i-npg+j*(j-1)+i*(i-1)+k-1)+1 ; inf = siavarpg*(npg*i-npg+j*(j-1)+i*(i-1)+k) ; iav = ini:inf;
            
            % Deformaci๓n
            % จจจจจจจจจจจ
            if struhyp == 3;
                [B,detJ] = matrixB_bbar_hexahedral(coord_n,[xg(i),xg(j),xg(k)],BH,e_VG);
            end
            eps_new(itv) = B*u;
            
            % Retono a la superficie de fluencia
            % จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
            if (conshyp == 2);
                [sigma_new(itv),hvar_new(ihv),aux_var(iav)] = rmapfi_plasJ2(eps_new(itv),...
                    hvar_old(ihv),Eprop,ce,e_VG);
            elseif (conshyp == 7)
                [sigma_new(itv),hvar_new(ihv),aux_var(iav)] = rmapfi_danio_plasdp_new(...
                    eps_new(itv),hvar_old(ihv),Eprop,ce,e_VG);
            end
            
            % Cแlculo de fint
            % จจจจจจจจจจจจจจจ
            fint = fint + B.'*sigma_new(itv)*wg(i)*wg(j)*wg(k)*detJ;
            
            % Cแlculo de la tension en los nodos
            % จจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจจ
            % E  = xg(i);
            % n  = xg(j);
            % z  = xg(k);
            % N1 = 1/8*(1-E)*(1-n)*(1-z);
            % N2 = 1/8*(1+E)*(1-n)*(1-z);
            % N3 = 1/8*(1+E)*(1+n)*(1-z);
            % N4 = 1/8*(1-E)*(1+n)*(1-z);
            % N5 = 1/8*(1-E)*(1-n)*(1+z);
            % N6 = 1/8*(1+E)*(1-n)*(1+z);
            % N7 = 1/8*(1+E)*(1+n)*(1+z);
            % N8 = 1/8*(1-E)*(1+n)*(1+z);
            % N  = [N1 N2 N3 N4 N5 N6 N7 N8];
            
            % for m = 1:length(N)
                % sigma_nodo(:,m) = sigma_nodo(:,m) + N(m)*sigma_new(itv)*wg(i)*wg(j)*wg(k)*detJ;
            % end
            
        end
    end
end