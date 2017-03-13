function BH = matrixBh_bbar_q1(coord_n,e_DatElemSet,e_VG)

dofpe = e_DatElemSet.dofpe;
npg = e_DatElemSet.npg;
wg = e_DatElemSet.wg;
xg = e_DatElemSet.xg;

Bh = zeros(1,dofpe);
BH = zeros(3,dofpe);
vol = 0;

for iPG = 1:npg
   [Bs,detJ] = matrixBs_bbar_q1(coord_n,xg(iPG,:),e_DatElemSet,e_VG);
   Bh = Bh+sum(Bs(1:3,:))*wg(iPG)*detJ;
   vol = vol+wg(iPG)*detJ;
end

Bh = Bh/vol;

for i = 1:3
    BH(i,:) = Bh;
end