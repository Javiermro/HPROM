function m_GdLNod = f_DofElem(m_ListNod,nGdLNod)

   %m_ListNod es un vector fila conteniendo los nodos de los cuáles se quiere obtener los grados de
   %libertad.
   %Devuelve un vector que contiene todos los grados de libertad de cada nodo ordenados en forma
   %contigua.
   m_GdLNod = reshape(bsxfun(@plus,nGdLNod*m_ListNod,(-nGdLNod+1:0)'),[],1);
    
end