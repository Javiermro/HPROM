function [c_CT,KT] = f_TensorTangDescarga(e_VarEst,e_DatSet,e_VG)

   %Esta función devuelve todos los tensores tangentes constitutivos de los puntos de gauss de
   %la malla y la matriz de rigidez global si se fuerza la descarga en ese instante de tiempos.
   nSet = e_VG.nSet;
   ntens = e_VG.ntens;
   c_CT = cell(nSet,1);
   
   for iSet = 1:nSet
      
      e_DatMatSet = e_DatSet(iSet).e_DatMat;
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      m_BT = e_DatSet(iSet).m_BT;
      m_DetJT = e_DatSet(iSet).m_DetJT;
      nElem = e_DatSet(iSet).nElem;
      npg = e_DatMatSet.npg;
      sihvarpg = e_DatElemSet.sihvarpg;
      hvar = reshape(e_VarEst(iSet).hvar,sihvarpg,npg,nElem); 
      
      m_CT = zeros(ntens,ntens,npg,nElem);
      
      for iElem = 1:nElem
         
         switch eltype            
            case {2,4,8,5,7}
               %B = 
               for iPG = 1:nPG
                  switch conshyp
                     case {1,2}  %Elasticidad lineal y Elasto-Plasticidad J2: Hardening-Softening Isotrópico
                        ct = e_DatMatSet.ce;
                     case 4  %Daño isótropo y Daño isótropo regularizado
                        ct = e_DatMatSet.ce*hvar_old(2,iPG,iElem)/hvar_old(1,iPG,iElem);
                     case {50,51}
                        hvarMacro = hvar(:,iPG,iElem);
                        %Modelo multiescala clásico y multiescala
                        c_CTMicro = f_TensorTangDescarga(hvarMicro.e_VarEst,...
                          e_DatMatSet.e_DatSet,e_DatMatSet.e_VG);
                        %Me faltaría almcenar el KT como variable auxiliar
                        ct = f_ModTangHomog(KT,c_CTMicro,hvarMacro.m_LinCond,hvarMicro.dofl,...
                           hvarMicro.doff,e_DatMatSet.e_DatSet,e_DatMatSet.Omega_micro,...
                           e_DatMatSet.e_VG);
                     otherwise
                        error(['Tensor tangente constitutivo en descarga: Modelo constitutivo no ',...
                           'implementado.'])
                  end
                  m_CT(:,:,iPG,iElem) = ct;
               end
            %case 10
            otherwise
               %Ver que la mayoría de los elementos no influyen en el tensor, excepto por ejemplo el
               %SDA.
               error(['Tensor tangente constitutivo en descarga: Elemento finito no ',...
                           'implementado.'])
         end
         
      end
      
      c_CT{iSet} = m_CT;
      
   end     
   
end