function  m_VecVoigt = f_Vec2Voigt2DVect(m_Vector,e_VG)

   % Vector en formato matricial, asumiendo notación de Voigt con el que el tensor Sigma y
   % Deformaciones teniendo la siguiente forma [Sxx;Syy;Szz;Sxy] (caso 2D, de tensión y deformación
   % plana).
   %Esta función está vectorizada, donde en dim=2 de m_Vector se indica las componentes x e y 
   %del vector, y en dim=1, se indica todos los vectores que se quiere obtener expresados en
   %notación de Voigt.
   %Devuelve una matriz, donde cada (:,:,i) es el Vector Voigt del vector m_Vector(i,:).   

   m_VecVoigt = zeros(e_VG.ntens,e_VG.ndime,size(m_Vector,1));
   m_VecVoigt(1,1,:) = m_Vector(:,1);  
   m_VecVoigt(4,1,:) = m_Vector(:,2);
   m_VecVoigt(2,2,:) = m_Vector(:,2); 
   m_VecVoigt(4,2,:) = m_Vector(:,1);  

end