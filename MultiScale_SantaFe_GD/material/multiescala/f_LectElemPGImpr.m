function m_ElemPGImpr = f_LectElemPGImpr(imprResCU,dirDat,m_NumElemSet)

   if ~isempty(imprResCU)&&~any(isnan(imprResCU))
      %Si se ingrensa el nombre del archivo a leer (tiene que tener extensión imprResCU).
      fiD = fopen(fullfile(dirDat,[imprResCU,'.imprResCU']));
      %En el archivo se debe indicar el número de elemento (el número de elemento debe corresponder
      %al utilizado en el archivo de datos, sin importar su posición en la lista de conectividades)
      %y el número de PG (esta numeración según el orden de las coordenadas indicadas en la matriz
      %xg) del problema macro. Una fila por PG.
      %No es necesario crear un archivo por cada set si se tiene varios set multiescala en la misma
      %estructura, ya los elementos se descartan si no pertenece al set.
      %Cuidado que no se descarta en el caso que se indique un punto de gauss que no está definido en el
      %elemento.
      m_ElemPGImpr = fscanf(fiD,'%f %f\n',[2,inf]);
      fclose(fiD);
      %Se obtiene la posición en la lista de elementos (denominación) indicados para
      %imprimir. Es transformar como se denomina los elementos en el archivo de datos (una lista no
      %necesariamente ordenada y completa, pero única) en su posición en la lista (matriz) de 
      %elementos del set.
      m_IndElem = bsxfun(@eq,m_NumElemSet',m_ElemPGImpr(1,:));
      [m_ElemImpr,~] = find(m_IndElem);
      %Se produce un error cuando se tiene un solo elemento en el set y se imprime ambos puntos de gauss, ya
      %que find en m_ElemImpr devuelve una matriz fila, distintas a más elementos donde devuelve una matriz
      %columna. Por ello se usa m_ElemImpr(:)' en lugar de solo m_ElemImpr'.
      m_ElemPGImpr = [m_ElemImpr(:)';m_ElemPGImpr(2,any(m_IndElem,1))];
   else
      %En el caso no estar definida la variable (NaN) o se indique con valor vacío (""), no se imprime ningún
      %elemento y punto de gauss.
      m_ElemPGImpr = [];
   end

end