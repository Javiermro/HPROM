function m_CondBord = f_CBLineal(m_CondBord,m_ConecFront2,ndn)

   %Se asume que elementos lineales de la conectividades forman una frontera cerrada, por lo que
   %tomando una columna de las conectividades se obtiene todos los nodos que forman la frontera.
   %Esta forma de determinar la lista de nodos que forma la frontera puede fallar si se tiene áreas separadas
   %que coinciden en un nodo (en ese caso habría nodos repetidos y se estaría imponiendo la misma CB al mismo
   %nodo). Usando unique debería resolverse este caso.
   m_NodFront = unique(m_ConecFront2(:,1));
   %Hay que determinar que gdl ya tiene condiciones de borde para que no incorporar las lineales,
   %esto puede hacer que no se verifique la CB de MR en la frontera. Esto ocurre en el caso que la
   %frontera interna tenga nodos que coincidan con nodos de la externa, que ya tienen impuesto
   %condiciones de borde, y como por convención no se pisa esas CB (ver que se puede hace mejor para
   %resolver esto), puede quedar la frontera interna no verificando la condición de homogenización
   %de deformaciones (MR).
   %Se podría pasar la matriz doff o los grados de libertad que tiene restricciones, pero al nivel 
   %que se piensa llamar esta función no se la tiene definida, por lo que se la calcula con los
   %datos de m_CondBord. Ver que ahora se dispone de los doff cuando se llama esta función, ver si cambiarla.
   %Se transforma la indicación única de grados de libertad restringidos en columnas separadas.
   %Se ordena en la filas ndn y en las columnas los nodos con restricción.
   m_DirRestr = zeros(ndn,size(m_CondBord,1));
   m_DirRestr(ndn,:) = m_CondBord(:,2);
   for iGdl = 1:ndn-1
      decGdl = 10^(ndn-iGdl);
      m_DirRestrGdl = fix(m_DirRestr(ndn,:)/decGdl);
      m_DirRestr(iGdl,:) = m_DirRestrGdl;
      m_DirRestr(ndn,:) = m_DirRestr(ndn,:)-m_DirRestrGdl*decGdl;
   end
   m_NodRest = m_CondBord(:,1);
   %Como se puede ingresar líneas con nodos repetidos pero que tenga grados de libertad
   %restringidos, y se determina que nodos y gdl están restringidos teniendo en cuenta esto.
   %Hay que tomar los nodos pertenecientes a la frontera interna y de esos tomar los grados de
   %libertad libres. A esos nodos y gdl se le impone la condición lineal (se puede agregar
   %directamente a m_CondBord, aunque haya nodos repetidos, lo que importa es que no haya 
   %gdl repetidos y es lo que verifica el programa).
   m_CondNoResFront = true(length(m_NodFront),ndn);
   m_IndCol = 1:ndn;
   for iNodRes = 1:length(m_NodRest)
      m_CondNoResFront(m_NodFront==m_NodRest(iNodRes),m_IndCol(m_DirRestr(:,iNodRes)'~=0)) = false;
   end
   %
   %Acá se determina los nodos que que tienen una sola restricción o no tienen ninguna, y están sobre la
   %segunda frontera.
   m_NodNoComplRest = ~all(~m_CondNoResFront,2);
   m_NodFront = m_NodFront(m_NodNoComplRest);
   m_CondNoResFront = m_CondNoResFront(m_NodNoComplRest,:);
   %Como m_CondNoResFront es una matriz lógica y la condición lineal se define poniendo 1, se
   %puede utilizar directamente para armar los datos de las CB.
   m_CondBordLineal = zeros(size(m_NodFront,1),2+ndn);
   m_CondBordLineal(:,1) = m_NodFront;
   m_CondBordLineal(:,2) = m_CondNoResFront(:,1)*10+m_CondNoResFront(:,2);
   
   m_CondBord = [m_CondBord;m_CondBordLineal];

end