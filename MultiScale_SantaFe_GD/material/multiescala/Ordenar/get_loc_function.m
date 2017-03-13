function [loc_function] = get_loc_function (fload,Omega_micro,Omega_micro_loc,e_VG)

if (e_VG.new_formulation == 1 && e_VG.istep_bif ~= 0 && e_VG.istep > e_VG.istep_bif+1 && fload == 1)
   loc_function = Omega_micro/Omega_micro_loc;
else
   loc_function = 0;
end