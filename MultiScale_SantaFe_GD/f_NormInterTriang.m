function m_NormElem = f_NormInterTriang(s_archDat,m_setInterf)
   %Función para el cálculo exacto de las normales de las interfaces formadas por elementos triangulares.
   %Tiene que agregado todos los path del programa.
   %La lógica que se usa es que el lado que tiene como elemento vecino un material distinto al de interfaz o
   %no tiene elementos (en la matriz de elemento vecinos debido al hueco que existe en las interfaces, es
   %decir aparece NaN).

   % Lectura de los datos de la malla
   [pathArchDat,nomArchDat,ext] = fileparts(s_archDat);
   nomArchDat = [nomArchDat,ext];
   [in,m_Coord,m_SetElem,~,~,e_DatSet,e_VG] = read_data(nomArchDat,pathArchDat);
   nSet = e_VG.nSet;
   m_Coord = m_Coord(:,1:2);
   % Determinación de los elementos vecinos
   m_ElemVec = f_ElemVecino(e_DatSet,e_VG);
   
   %Determinación de las normales por set
   c_NormElemSet = cell(nSet,1);
   %
   close all
   hold on
   for iSet = 1:nSet      
      if any(iSet==m_setInterf)
         
         m_ElemVecSet = m_ElemVec(e_DatSet(iSet).m_IndElemSet,:);
         m_SetVecSet = nan(size(m_ElemVecSet));
         m_SetVecSet(~isnan(m_ElemVecSet)) = m_SetElem(m_ElemVecSet(~isnan(m_ElemVecSet)));
         %
         %Se traspone para indexar correctamente con la matriz lógica.
         m_ConecSet = e_DatSet(iSet).conec';
         %
         %Matriz de índices del nodo inicial del lado que está en contacto con un set distinto al de interface
         %(y no es NaN). Se trasponde para indexar correctamente con la matriz lógica (es decir que el 
         %resultado salga ordenado según la numeración de elementos del set).
         m_IndNodInic = (m_SetVecSet~=iSet&~isnan(m_SetVecSet))';
         %Verificación:
         %- Debería haber un solo lado donde ocurra eso, sino hay error en la malla ya que no
         %coincide con la hipótesis de este programa.
         m_Verif = sum(m_IndNodInic,1);
         if any(m_Verif~=1)
            error(['Normales de interfaces triangulares: Hay un elemento de interfaz que tiene más ',...
               'de un lado en contacto con materiales con sets distintos o no tiene ninguno.'])
         end
         %
         m_VecLad = m_Coord(m_ConecSet(circshift(m_IndNodInic,[1,0])),:)-m_Coord(m_ConecSet(m_IndNodInic),:);
         %Gráfica de vectores de los lados.
         quiver(m_Coord(m_ConecSet(m_IndNodInic),1),m_Coord(m_ConecSet(m_IndNodInic),2),...
            m_VecLad(:,1),m_VecLad(:,2),0)
         axis equal
         %Longitud del lado
         m_LongLado = hypot(m_VecLad(:,1),m_VecLad(:,2));
         %Se normaliza los vectores.
         m_VecLad = bsxfun(@rdivide,m_VecLad,m_LongLado);
         %Se determina los vectores normales.
         %Se ordena la matriz de salida de la siguiente forma: NroElem CompXNorm CompYNorm LongElem
         m_NormElemSet = [e_DatSet(iSet).m_NumElem',m_VecLad(:,2),-m_VecLad(:,1),m_LongLado];
         %Gráfica de los vectores normales
         %Se adopta la dirección de vectores salientes al elemento.
         m_PosMediaLado = (m_Coord(m_ConecSet(circshift(m_IndNodInic,[1,0])),:)+...
            m_Coord(m_ConecSet(m_IndNodInic),:))/2;
         quiver(m_PosMediaLado(:,1),m_PosMediaLado(:,2),m_NormElemSet(:,2),m_NormElemSet(:,3),0)
         
      else
            
         m_NormElemSet = zeros(0,4);
         
      end
      %
      c_NormElemSet{iSet} = m_NormElemSet;
   end
   hold off
   %
   %
   m_NormElem = cat(1,c_NormElemSet{:});

end