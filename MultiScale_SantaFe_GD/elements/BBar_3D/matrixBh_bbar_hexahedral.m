function [BH] = matrixBh_bbar_hexahedral(coord_n,e_VG)

%global ndn dofpe npg wg
ndn = e_VG.ndn;
dofpe = e_VG.dofpe;
npg = e_VG.npg;
wg = e_VG.wg;
xg = e_VG.xg;

Bh  = zeros(1,dofpe);
BH  = zeros(3,dofpe);
vol = 0;

for i = 1:npg
   for j = 1:npg
      for k = 1:npg
         [Bs,detJ] = matrixBs_bbar_hexahedral(coord_n,[xg(i),xg(j),xg(k)],e_VG);
         Bh  = Bh + sum(Bs(1:3,:))*wg(i)*wg(j)*wg(k)*detJ;
         vol = vol + wg(i)*wg(j)*wg(k)*detJ;
      end
   end
end

Bh = Bh/vol;

for i = 1:3
    BH(i,:) = Bh;
end