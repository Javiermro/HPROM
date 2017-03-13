function [ct,sigma_new,eps_new,hvar_new,aux_var] = rmap_damage_reg(...
   eps_new,hvar_old,e_DatMatSet,e_VG)

   %Se regulariza el modelo constitutivo según el tamaño de elemento.
   tita = e_DatMatSet.tit;
   %Ver si no pasar la matriz hReg de otra forma, ya que depende explícitamente que se pase en el
   %e_VG el número de elemento que le corresponde en el set, algo se había tratado de evitar a este
   %nivel (solo se usaba para debug o imprensión).
   %También se podría interpretar que es una propiedad que depende de la posición donde se evalúa, y
   %donde se asume que la posición se indica como el número de elemento finito.
   hElemReg = e_DatMatSet.m_hReg(e_VG.iElemSet);
   switch tita
      case 0
         % Lineal            
         if e_DatMatSet.young~=0
            %H = -Sigmau^2/2/E/Gf*h;
            e_DatMatSet.hba = -e_DatMatSet.ftult^2/2/e_DatMatSet.young/e_DatMatSet.gfv*hElemReg;
         else
            e_DatMatSet.hba = 0;
         end
         [ct,sigma_new,eps_new,hvar_new,aux_var] = rmap_damage(eps_new,hvar_old,e_DatMatSet,...
            e_VG);
      case 1
         % Exponencial         
         %Para considerar que en el caso se utilice E=0, y el modelo constitutivo funcione
         %(igualmente hay que resolver que no queden nodos sin rigidez)
         if e_DatMatSet.young~=0
            %H = -Sigmau^2/E/Gf*h;
            e_DatMatSet.hba = -e_DatMatSet.ftult^2/e_DatMatSet.young/e_DatMatSet.gfv*hElemReg;
         else
            e_DatMatSet.hba = 0;
         end
         [ct,sigma_new,eps_new,hvar_new,aux_var] = rmap_damage(eps_new,hvar_old,e_DatMatSet,...
            e_VG);
      otherwise
         error('Modelo de daño regularizado: Modelo de ablandamiento no definido.')
   end
   
end