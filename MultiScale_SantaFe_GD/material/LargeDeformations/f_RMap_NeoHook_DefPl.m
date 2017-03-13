function [m_DPF,m_P] = f_RMap_NeoHook_DefPl(m_F,e_DatMatSet,e_VG)

   %Se está considerando el caso especial de deformación plana, donde Fzz es 1, pero para considerar el
   %ruido que agregan los elementos BBar, se considera Fzz no es necesariamente 1.
   %Los términos de corte en z (Fxz,Fzx,Fyz,Fzy) se siguen considerando nulo (por la hipótesis de deformación
   %plana Fzz=1,Fzx=0,Fyz=0,Fzy=0).
   
   %
   %nTens = e_VG.ntens;
   
   %Constantes
   %Módulo de compresibilidad logarítmico.
   K = e_DatMatSet.modCompr;
   %Módulo de corte.
   G = e_DatMatSet.modCorte;
   
   %Se transforma la matriz de gradiente (F) en notación de Voigt a matriz normal, pero solo los términos x-y.
   m_F2D = [m_F(1),m_F(4);m_F(5),m_F(2)];
   Fzz = m_F(3);
   %m_F3D = [m_F2D,zeros(2,1);zeros(1,2),Fzz]
   
   %Inversa de la matriz (F)
   %m_InvF2D = inv(m_F2D);
   m_InvF2D = m_F2D\eye(2,2);
   %Si no se produce deformación en Fzz viene la identidad.
   InvF2Dzz = 1/Fzz;
   
   %Determinante del gradiente de deformación.
   %J = m_F(3)*(m_F(1)*m_F(2)-m_F(4)*m_F(5));
   J = Fzz*det(m_F2D);
   
   %Tensor izquierdo de deformación de Cauchy-Green (la componente zz es Fzz^2).
   m_B2D = m_F2D*m_F2D';
   B2Dzz = Fzz^2;
   
   %Tensor isocórico del tensor de deformación de Cauchy-Green (la componente zz es J^(-2/3)*Fzz^2)
   m_B2DIso = J^(-2/3)*m_B2D;
   B2DIsozz = J^(-2/3)*B2Dzz;   
   
   %Tensor desviador del tensor isocórico del tensor de deformación de Cauchy-Green.
   m_DesvB2DIso = [(2*m_B2DIso(1,1)-m_B2DIso(2,2)-B2DIsozz)/3,m_B2DIso(1,2);...
      m_B2DIso(2,1),(2*m_B2DIso(2,2)-m_B2DIso(1,1)-B2DIsozz)/3];
   DesvB2DIsozz = (2*B2DIsozz-m_B2DIso(1,1)-m_B2DIso(2,2))/3;
   
   %Tensor de tensiones de Kirchhoff.
   m_Tau = G*m_DesvB2DIso+K*log(J)*eye(2,2);
   Tauzz = G*DesvB2DIsozz+K*log(J);
   
   %Tensor de tensiones Primero de Piola-Kirchhoff.
   m_P = m_Tau*m_InvF2D';
   Pzz = Tauzz*InvF2Dzz;
   m_P = [m_P(1,1);m_P(2,2);Pzz;m_P(1,2);m_P(2,1)];
   
   %Tensor de tensiones Segundo de Piola-Kirchoff
   %m_S2D = m_InvF2D*m_P;
   %Szz = InvF2Dzz*Pzz;
   
   %Tensor tangente constitutivo
   m_DPF = [];
   
end