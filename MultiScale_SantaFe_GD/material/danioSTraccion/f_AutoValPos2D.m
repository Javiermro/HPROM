function m_SigmaPos = f_AutoValPos2D(m_Sigma,e_VG)

   %La verificación de que si es tensor de tensiones es nulo y imponer que angle=0,sigmacri=0 no es
   %necesario ya que atan2 si es 0/0 (todas las tensiones nulas) devuelve angle=0.
   %Tampoco es necesario considerar el caso en que las dos tensiones principales son iguales
   %(TAU1==TAU2), ya que atan2 con 0/0 devuelve angle 0.
   %Es necesario usar atan2 en lugar de atan porque debe ser unívoca la solución entre los ángulos 
   %[-pi/2,pi/2] para abarcar todas las direcciones posibles de las normales (considerando que la 
   %normal en una dirección abarca ambos sentidos), ya que atan2 devuelve entre [-pi,pi], y al 
   %dividir por 2, se obtiene el rango [-pi/2,pi/2] buscado. 

   TAUXX = m_Sigma(1);
   TAUYY = m_Sigma(2);
   TAUZZ = m_Sigma(3);
   TAUXY = m_Sigma(4);

   ADDTAU = 0.5*(TAUXX+TAUYY);
   REDTAU = TAUXX-TAUYY;
   %Es más estable numéricamente usar hypot, pero no sé cuánto más lento es.
   %raiz = sqrt(0.25*REDTAU*REDTAU+TAUXY*TAUXY);
   raiz = hypot(REDTAU/2,TAUXY);
   %Tensiones principales Sigma 1 y Sigma 2 en el plano xy.
   %(en el caso que de tensión plana, SigmaZ es nula por definición, así que no influye, pero en deformación
   %plana sigmaZ es una tensión principal ya que los esfuerzos de corte sigmaXZ y sigmaYZ son nulos).
   %Lo que no se sabe correctamente cuál es la relación TAU1>TAU2>TAU3 para todos los modelos constitutivos.
   TAU1 = ADDTAU+raiz;
   TAU2 = ADDTAU-raiz;
   TAU3 = TAUZZ;
   %
   %Se utiliza para pasar de la base de los autovectores a la base global del problema.
   angle = 0.5*atan2(2*TAUXY,REDTAU);
   cosAng = cos(angle);
   senAng = sin(angle);
   m_SigmaPos = zeros(e_VG.ntens,1);
   if TAU1>0
      m_SigmaPos(1) = TAU1*cosAng^2;
      m_SigmaPos(2) = TAU1*senAng^2;
      m_SigmaPos(4) = TAU1*(senAng*cosAng);
   end
   if TAU2>0
      m_SigmaPos(1) = m_SigmaPos(1)+TAU2*senAng^2;
      m_SigmaPos(2) = m_SigmaPos(2)+TAU2*cosAng^2;
      m_SigmaPos(4) = m_SigmaPos(4)-TAU2*(senAng*cosAng);
   end
   if TAU3>0
      m_SigmaPos(3) = TAU3;
   end
    
end
         
