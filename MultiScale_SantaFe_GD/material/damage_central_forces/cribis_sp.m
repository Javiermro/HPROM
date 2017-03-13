function [angle,sigmacri] = cribis_sp(strsg)
%***********************************************************************
%
%           THIS ROUTINE DETERMINES THE CRITICAL CONDITIONS 
%
%    Input variables:
%       strsg                 Total strains at Gauss point NGAUL
%
%    Output variables: 
%       angle                 Localization angle
%       sigmacri              Critical hardening modulus
%***********************************************************************
%La verificación de que si es tensor de tensiones es nulo y imponer que angle=0,sigmacri=0 no es
%necesario ya que atan2 si es 0/0 (todas las tensiones nulas) devuelve angle=0.
%Tampoco es necesario considerar el caso en que las dos tensiones principales son iguales
%(TAU1==TAU2), ya que atan2 con 0/0 devuelve angle 0.
%Es necesario usar atan2 en lugar de atan porque debe ser unívoca la solución entre los ángulos 
%[-pi/2,pi/2] para abarcar todas las direcciones posibles de las normales (considerando que la 
%normal en una dirección abarca ambos sentidos), ya que atan2 devuelve entre [-pi,pi], y al 
%dividir por 2, se obtiene el rango [-pi/2,pi/2] buscado. 
   
   %[stpos,TAU1P,TAU2P]=posi2d(strsg);
   TAUXX = strsg(1);
   TAUYY = strsg(2);
   TAUXY = strsg(4);
   %TAUZZ = strsg(3);

   ADDTAU = 0.5*(TAUXX+TAUYY);
   REDTAU = TAUXX-TAUYY;
   %Es más estable numéricamente usar hypot, pero no sé cuánto más lento es.
   %raiz = sqrt(0.25*REDTAU*REDTAU+TAUXY*TAUXY);
   raiz = hypot(REDTAU/2,TAUXY);
   %Tensiones principales Sigma 1 y Sigma 2 en el plano xy (en el caso que de tensión plana, SigmaZ
   %es nula por definición, así que no influye).
   TAU1 = ADDTAU+raiz;
   %TAU2 = ADDTAU-raiz;
   %
   angle = 0.5*atan2(2*TAUXY,REDTAU);
   sigmacri = TAU1;
    
end
         
