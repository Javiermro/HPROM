function m_hReg = f_LongReg(hReg,nomArchHReg,dirDat,m_NumElemSet,nElemSet)

   %Se utiliza NaN para verificar que todos los elementos tiene indicados como proceder para el cálculo de la
   %distancia de regularización.
   %automáticamente, según algún criterio posterior y impuesto en el archivo de datos.
   m_hReg = nan(1,nElemSet);
   if ~isnan(hReg)
      %En hReg se puede indicar la distancia de regularización, ó en valor negativo el método a usar para el
      %cálculo automático. Este valor se aplica a todos los elementos.
      m_hReg(:) = hReg;
   end
   %En la propiedad ARCHHREG se puede indicar "", quedando la variable como un string vacío.
   if ~isempty(nomArchHReg)&&~any(isnan(nomArchHReg))
      %Se lee un archivo que se indica las distancias de regularización. Esta sobreescriben las indicadas por
      %el hReg.
      %Si se ingresa el nombre del archivo a leer (tiene que tener extensión hReg).
      fiD = fopen(fullfile(dirDat,[nomArchHReg,'.hReg']));
      %En el archivo se debe indicar el número de elemento (el número de elemento debe corresponder
      %al utilizado en el archivo de datos, sin importar su posición en la lista de conectividades)
      %y el distancia de regularización (una fila por cada número de elemento finito).
      %No es necesario crear un archivo por cada set si se tiene varios set de daño en la misma
      %estructura, ya los elementos se descartan si no pertenece al set.
      m_ElemLeido = fscanf(fiD,'%f %f\n',[2,inf]);
      fclose(fiD);
      %Se obtiene la posición en la lista de elementos del set de los elementos (denominación) 
      %indicados.
      %Es transformar como se denomina los elementos en el archivo de datos (una lista no
      %necesariamente ordenada y completa, pero única) en su índice posición dentro de la lista de 
      %elementos del set (que es lo que se utiliza dentro del programa para indentificar al
      %elemento).
      m_IndElem = bsxfun(@eq,m_NumElemSet',m_ElemLeido(1,:));
      [m_NroIndhReg,~] = find(m_IndElem);
      %Se reemplaza únicamente los elementos que se indicó en el archivo ARCHHREG, en el resto se impone el
 
     %valor de HREG o por defecto. Cuidado no tira ninguna notificación si en el archivo ARCHHREG se olvida
      %de ingresar algún elemento.
      m_hReg(m_NroIndhReg) = m_ElemLeido(2,any(m_IndElem,1));
   end
   %
   if any(isnan(m_hReg))
      error(['Lectura de datos: Propiedades Materiales: Modelo constitutivo ',...
         'regularizado: Longitud de regularización: Hay elementos finitos que no tiene definido la distancia.'])
   end

end