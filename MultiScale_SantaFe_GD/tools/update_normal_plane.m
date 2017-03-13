function [du_iter,dlambda,iconv] = update_normal_plane(res,KT,Fext,dofl,Du_step,Du_prev,e_VG)

%global istep iter ds
ds = e_VG.ds;
%dofl = e_VG.dofl;
istep = e_VG.istep;
iter = e_VG.iter;
 
iconv = 1;

u   = KT(dofl,dofl)\[-(res) Fext(dofl)];
du1 = u(:,1);
du2 = u(:,2);

if ((istep*iter) == 1);
   dlambda = 0.0;
else
   if (iter == 1); 
      dlambda = ds^2/dot(du2,du2);
      if (dot(du2,Du_prev(dofl)) < 0)
         dlambda = -dlambda;
      end
   else
      dlambda =-(DOT(Du_step(dofl),du1)/DOT(Du_step(dofl),du2));
   end
end

du_iter = du1 + dlambda*du2;

return