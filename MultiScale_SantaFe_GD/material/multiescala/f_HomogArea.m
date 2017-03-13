function m_valHomog = f_HomogArea(c_valElem,nVar,volTotal,c_DetJT,e_DatSET,e_VG)

   %El espesor o el área viene con los pesos de gauss, y por lo tanto no es necesario considerar
   %en la integración (ya viene como peso "volumétrico").

   %%
   %Función que devuelve los valores homogenizados a partir de la integración en toda la superficie
   %de la estructrura.
   %En la celda c_valElem se ingresa las matrices separadas de cada set. Esta matriz debe ser de dimensiones
   %(nVar,npg,nElem) o (nVar*npg,nElem), donde nVar la cantidad de variables o la dimensión de la variable que
   %se quiere homogenizar. Esta nVar debe ser la misma para todos los set.
%    nElem = e_VG.nElem;
   nSet = e_VG.nSet;

   m_valHomog = zeros(nVar,1);
   for iSet = 1:nSet
      %nElem = e_DatSET(iSet).nElem;
      wg = e_DatSET(iSet).e_DatElem.wg;
      npg = e_DatSET(iSet).e_DatElem.npg;
      m_valElemSet = reshape(c_valElem{iSet},nVar,npg,[]);
      %
      m_valHomog = m_valHomog+sum(sum(bsxfun(@times,m_valElemSet,...
         reshape(bsxfun(@times,c_DetJT{iSet},wg),1,npg,[])),2),3);
   end
   m_valHomog = m_valHomog/volTotal;
      
end
     