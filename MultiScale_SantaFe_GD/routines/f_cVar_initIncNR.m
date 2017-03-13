function c_GdlCond = f_cVar_initIncNR(e_DatSet,e_VG, c_GdlCond)

   nSet  = e_VG.nSet;
   ndime = e_VG.ndime;
   ntens = e_VG.ntens;
   for iSet = 1:nSet
      nElem = e_DatSet(iSet).nElem;
      eltype = e_DatSet(iSet).e_DatElem.eltype;
      switch eltype
         case 10
            %Se inicializa la variable de salto.
            c_GdlCond{iSet,1} = zeros(ndime,nElem);
         case {21,22,23}
            %Se inicializa la variable de salto.
            c_GdlCond{iSet,2} = zeros(ndime,nElem);   %   Dbeta_SDA
            c_GdlCond{iSet,4} = zeros(ntens,nElem);   %   Dbeta_MSL  
      end
   end

end