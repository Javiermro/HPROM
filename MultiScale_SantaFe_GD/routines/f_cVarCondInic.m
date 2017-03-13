function c_GdlCond = f_cVarCondInic(e_DatSet,e_VG)

   nSet  = e_VG.nSet;
   ndime = e_VG.ndime;
   ntens = e_VG.ntens;
   %Cantidad máximas (para todos los tipos de set) de matrices que se van almancenar en la celda.
   nMax = 4;
   %En esta celda se coloca n matrices (b es el grado de libertad interno que se condensa):
   c_GdlCond = cell(nSet,nMax);
   for iSet = 1:nSet
      nElem = e_DatSet(iSet).nElem;
      eltype = e_DatSet(iSet).e_DatElem.eltype;
      switch eltype
         case 10
            %c_GdlCond{iSet,1} = m_Valores_b;
            %c_GdlCond{iSet,2} = m_Residuo_b;
            %c_GdlCond{iSet,3} = m_MatrizElemental_ub;
            %c_GdlCond{iSet,4} = (m_MatrizElemental_bb)^-1;   

            %Se inicializa la variable de salto.
            c_GdlCond{iSet,1} = zeros(ndime,nElem);
         case {21,22,23}
            %Se inicializa la variable de salto.
            c_GdlCond{iSet,1} = zeros(ndime,nElem);   %   beta_SDA 
            c_GdlCond{iSet,2} = zeros(ndime,nElem);   %   Dbeta_SDA
            c_GdlCond{iSet,3} = zeros(ntens,nElem);   %   beta_MSL
            c_GdlCond{iSet,4} = zeros(ntens,nElem);   %   Dbeta_MSL  
      end
   end

end