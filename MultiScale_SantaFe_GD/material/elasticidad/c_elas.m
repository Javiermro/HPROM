function Ce = c_elas(E,poiss,e_VG)

   %global_var;
   struhyp = e_VG.struhyp;
   SSOIT = e_VG.SSOIT;
   FOAT1 = e_VG.FOAT1;

   mu    = E/(2*(1+poiss));
   lamda = (poiss*E)/((1+poiss)*(1-2*poiss)); 
   capa  = E/(3*(1-2*poiss));

   switch struhyp
      case 1                 % Estado plano de deformacion   
         Ce = lamda*SSOIT+2*mu*FOAT1;
      case 2                 % Estado plano de tension
         %Buscar como es expresado con mu y lambda, y si es posible con los tensores.
         Ce = E/(1-poiss^2)*[1,poiss,0,0;poiss,1,0,0;0,0,0,0;0,0,0,(1-poiss)/2];
      case 3                 % Estado tridimensional
         Ce = lamda*SSOIT+2*mu*FOAT1;
      case 4                 % Estado axisimetrico (no se si esta bien)
         Ce = lamda*SSOIT+2*mu*FOAT1;
      case 5                 % Barras
         Ce = E;
      otherwise
         error('Tensor tangente elástico: Hipótesis estructural no implementada');
   end
   
end