function norm_res = f_cVar_normRes(e_DatSet,e_VG, c_GdlCond)
   norm_res=0 ;
   nSet  = e_VG.nSet;
   for iSet = 1:nSet
      eltype = e_DatSet(iSet).e_DatElem.eltype;
      switch eltype
         case 10
            %Se inicializa la variable de salto.
            norm_res= norm_res+ norm(c_GdlCond{iSet,2});
         case 21
            norm_res= norm_res+ norm(c_GdlCond{iSet,7}) ... ;  %  m_Res_beta_SDA
                              + norm(c_GdlCond{iSet,10})  ;  % m_Res_beta_MSL
      end
   end
end