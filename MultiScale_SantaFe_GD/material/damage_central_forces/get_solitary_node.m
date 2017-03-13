function [n,grad_phi] = get_solitary_node(n,B,e_DatElem,e_VG)

   % Recupera variables globales
   ndn = e_VG.ndn;
   %npe = e_VG.npe;
   %dofpe = e_VG.dofpe;
   %ntens = e_VG.ntens;
   ndime = e_VG.ndime;
   npe = e_DatElem.npe;
   dofpe = e_DatElem.dofpe;

   % Matriz de derivadas de funciones de forma (se extraen de la matriz de deformaciones)
   %En el 3D habría que cambiar para que extraer la derivada en z también.
   dN_xy = zeros(ndime,npe);
   dN_xy(1,:) = B(1,1:ndn:dofpe);
   dN_xy(2,:) = B(2,2:ndn:dofpe);

   % Búsqueda del nodo solitario
   dDirNxy = n'*dN_xy;
   dDirNxy(1) = dDirNxy(1)/norm(dN_xy(:,1));
   dDirNxy(2) = dDirNxy(2)/norm(dN_xy(:,2));
   dDirNxy(3) = dDirNxy(3)/norm(dN_xy(:,3)); 

   [~,nsoli] = max(abs(dDirNxy));

   if dDirNxy(nsoli)<0
      n = -n;
   end

   grad_phi = dN_xy(:,nsoli);
   
end