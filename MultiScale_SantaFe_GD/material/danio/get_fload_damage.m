function [PGLOC] = get_fload_damage(hvar_old,e_VG)

sihvarpg = e_VG.sihvarpg;
sihvare  = e_VG.sihvare;
npg      = e_VG.npg;
nElem    = e_VG.nElem;

%Se devuelve los PG que localizaron y siguen cargando.
PGLOC = hvar_old(sihvarpg:sihvarpg:sihvare,:)>0;
