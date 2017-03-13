function [PGLOC] = get_fload_plasJ2bl (hvar_old,e_VG)

sihvarpg = e_VG.sihvarpg;
sihvare  = e_VG.sihvare;
ntens    = e_VG.ntens;
npg      = e_VG.npg;
nElem    = e_VG.nElem;

PGLOC              = zeros(npg,nElem);
fload              = hvar_old(ntens+2:sihvarpg:sihvare,:);
efload             = sum(fload);
PGLOC(:,efload~=0) = 1;