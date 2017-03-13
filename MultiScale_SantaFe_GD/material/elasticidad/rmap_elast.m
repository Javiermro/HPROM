function [ct,sigma_new] = rmap_elast(eps_new,ce)

   sigma_new = ce*eps_new;
   ct = ce;
   
end