function m_VarAuxElem = f_OperPosConv_MixStrInj_quad_q1(m_VarAuxElem,m_VarAuxGP,e_DatMatSet,m_CT,e_VG)

   %parfor iElem = 1:nElem
%    for iElem = 1:nElem
%       
%       condBif = m_VarAuxElem(1);
%       if ~condBif
%          %Para el análisis de bifurcación se toma el tensor de punto de gauss 5 para el implícito o 10 para el
%          %caso del implex (se asume que este elemento finito siempre tiene 5 puntos de gauss).
%          if esImplex
%             %Tensores tangentes constitutivos para bifurcación (se toma el tensor implícito).
%             m_CTBif = m_CT(:,:,10,iElem);
%          else
%             m_CTBif = m_CT(:,:,5,iElem);
%          end
%          condBif = f_CondBifct(m_CTBif,e_VG);
%       end
%       m_VarAuxElem(1) = condBif;
%       
%    end

   %Por ahora se realiza este análisis usando el fload.
   conshyp = e_DatMatSet.conshyp;
   switch conshyp
      case {4,5,10,12}
         siavarpg = e_DatMatSet.siavarpg;
         %Se considera que este elemento finito siempre tiene 5 puntos de gauss, y que este se utiliza para
         %determinar si el elemento completo bifurcó. El fload se considera que está en la posición 2 de estos
         %modelos constitutivos.
         %Se convierte automáticamente a double los valores lógicos m_VarAuxGP(4*siavarpg+2,:)>0 al ser
         %m_VarAuxElem una matriz de doubles.
         %Notar que se considera que pertenece al dominio mixto con cualquier valor positivo, y cuando hay
         %carga elástica (0) y descarga elástica (-1) queda excluido. Esta determinación se debe realizar en
         %todos los pasos, aún después de haber detectado la bifurcación ya que determinar la descarga.
         m_VarAuxElem(1,:) = m_VarAuxGP(4*siavarpg+2,:)>0;
      otherwise
         error(['Operaciones de fin de paso: Elemento finito quadrilátero Q1 mixto con inyección ',...
            'de deformación: Modelo constitutivo no definido'])
   end

end