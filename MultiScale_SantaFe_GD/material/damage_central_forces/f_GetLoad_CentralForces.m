function PGLOC = f_GetLoad_CentralForces(hvar_old,e_VG)

   % Determinación de los elementos localizados
   %Se asume que si los fload son distintos de cero están localizados, aunque esté en descarga.
   
   %Se recupera variables
   sihvarpg = e_VG.sihvarpg;
   sihvare  = e_VG.sihvare;
   %npg      = e_VG.npg;
   %nElem    = e_VG.nElem;
   
   %En este caso de CentralForces generalmente no se va se usar más de un punto de gauss, ya que
   %está hecho para triángulos lineales.
   %%Si se considera que el elemento se considera localizado si uno de los puntos de gauss lo está.
   %Luego a cada punto de gauss lo pone como localizado en la matriz PGLOC, repetiendo lo que le
   %corresponde en el elemento.
   %PGLOC = false(npg,nElem);
   %PGLOC(:,any(hvar_old(sihvarpg:sihvarpg:sihvare,:))) = true;
   %%Para obtener los PG localizados, individualmente
   PGLOC = hvar_old(sihvarpg:sihvarpg:sihvare,:)>0;

end