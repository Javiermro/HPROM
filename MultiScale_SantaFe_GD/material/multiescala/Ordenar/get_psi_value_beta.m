function [psi_value_beta,psi_value_beta_old,psi_value] = get_psi_value_beta...
            (fun_beta,time,e_VG,psi_value,psi_value_old,psi_value_beta_old,NMF,paso0)

if (NMF~=0)
   time_correct = time - paso0*e_VG.Dtime;
   psi_value_beta = get_psi_value(fun_beta,time_correct,e_VG);
   psi_value = psi_value_old-0.0001/1000;
   if (psi_value < 0)
       psi_value = 0;
   end
else
   psi_value_beta = 0;
end