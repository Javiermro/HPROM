function hvar_new = f_ReCalcNormal(sigma_new,hvar_new,e_DatMatSet,e_VG)

   %% REVISAR ESTA FUNCIÓN POR LOS CAMBIOS EN EL MANEJO DE LOS SET.

   %Luego que bifurca y de la convergencia del Newton se recalcula la normal de la fisura (es decir
   %su dirección) para que considerar una cierta adaptación de la misma hasta cierto criterio (en
   %este caso abertura de la fisura). Recordar que a la normal igualmente hay que calcularla dentro
   %del modelo constitutivo, ya que angle viene como variable histórico. Esto se debe realizar hasta
   %que bifurca, luego se mantiene inamovible, excepto para los casos considerados en esta función.
   
   sihvarpg = e_VG.sihvarpg;
   sihvare = e_VG.sihvare;
   npg = e_VG.npg;
   ntens = e_VG.ntens;
   
   ftult = e_DatMatSet.ftult;
   gfval = e_DatMatSet.gfv;

   %Posición en las matrices históricas
   posFload = 9;
   posXi = 3;
   posNormX = 7;
   posNormY = 8;
   posAng = 5;
   %posQBifu = 4;
   %posSigmCrit = 10:13;
   
   ratio = 0.0038;
   %Límite de la variación de la normal, según el módulo del salto histórico.
   m_LimXi = ratio*(gfval./ftult);
   %Puntos de gauss y elementos donde hay que recalcular la normal.
   m_PGNorVar = hvar_new(posFload:sihvarpg:sihvare,:)==1&...
      bsxfun(@le,hvar_new(posXi:sihvarpg:sihvare,:),m_LimXi');
      
   %En el caso de no tener ninguna bifurcación que esté en la etapa de normal variable, no se
   %realiza ninguno de los cálculos siguientes.
   if any(m_PGNorVar)
      %Determinación de las direcciones principales
      m_IndTens = 1:ntens:npg*ntens;
      %atan2 devuelve 0, si es atan2(0,0), y pi/2 si es atan2(·/0) con · distinto de cero
      m_Angulo = atan2(2*sigma_new(m_IndTens+(ntens-1),m_PGNorVar),...
         (sigma_new(m_IndTens,m_PGNorVar)-sigma_new(m_IndTens+1,m_PGNorVar)))/2;
      %Lo adelante está pensado para un solo punto de gauss y ntens igual a 4.
      %Actualización de las variables históricas
      %Normal
      hvar_new(posNormX,m_PGNorVar) = cos(m_Angulo);
      hvar_new(posNormY,m_PGNorVar) = sin(m_Angulo);
      %Ángulo
      hvar_new(posAng,m_PGNorVar) = m_Angulo;
      %Determinación de un nuevo límite del modelo constitutivo, al variar la normal de la fisura.
      %Se utilizar el tensor en el momento de la bifurcación, con la nueva normal.      
%       m_n1 = hvar_new(posNormX,m_PGNorVar);
%       m_n2 = hvar_new(posNormY,m_PGNorVar);
%       m_SigmaCrit = hvar_new(posSigmCrit,m_PGNorVar);      
%       hvar_new(posQBifu,m_PGNorVar) = hypot(m_n1.*m_SigmaCrit(1,:)+m_n2.*m_SigmaCrit(4,:),...
%             m_n2.*m_SigmaCrit(1,:)+m_n1.*m_SigmaCrit(4,:));      
   end

end

