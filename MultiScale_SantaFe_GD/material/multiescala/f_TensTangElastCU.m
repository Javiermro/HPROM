function [c_CT,KT] = f_TensTangElastCU(e_DatSet,e_VG)

   %Esta función devuelve todos los tensores tangentes constitutivos de los puntos de gauss de la malla y la
   %matriz de rigidez global si tiene un comportamiento elástico, sin considerar si es daño o plasticidad. Es
   %decir se considera todos los modelos constitutivos, pero se considera solo sus propiedades elásticas. Esto
   %permite utilizar un archivo de datos de la microcelda para elementos macro que se quiere que bifurquen y
   %para los elmentos que se quiere se comporten siempre elásticamente.
   
   nSet = e_VG.nSet;
   ntens = e_VG.ntens;
   ndn = e_VG.ndn;
   %
   c_CT = cell(nSet,1);
   c_Ke = cell(nSet,1);
   c_Fil = cell(nSet,1);
   c_Col = cell(nSet,1);
   
   for iSet = 1:nSet
      
      e_DatMatSet = e_DatSet(iSet).e_DatMat;
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      %conec = e_DatSet(iSet).conec;
      m_DofElem = e_DatSet(iSet).m_DofElem;
      m_BT = e_DatSet(iSet).m_BT;
      m_DetJT = e_DatSet(iSet).m_DetJT;
      nElem = e_DatSet(iSet).nElem;
      nPG = e_DatElemSet.npg;
      wg = e_DatElemSet.wg;
      dofpe = e_DatElemSet.dofpe;
      eltype = e_DatElemSet.eltype;
      conshyp = e_DatMatSet.conshyp;
      %
      m_Ke = zeros(dofpe,dofpe,nElem);      
      m_CT = zeros(ntens,ntens,nPG,nElem);
      % Grados de libertad y coordenadas de los nodos de los elementos del set
      %dofElemSet = f_DofElem(reshape(conec',1,[]),ndn);
      dofElemSet = m_DofElem(:);
      m_Fil = reshape(repmat(reshape(dofElemSet,dofpe,[]),dofpe,1),1,[]);
      m_Col = reshape(repmat(dofElemSet',dofpe,1),1,[]); 
      
      %parfor iElem = 1:nElem
      for iElem = 1:nElem
         
         switch eltype            
            case {2,4,8,31,32}
               
               m_kt = zeros(dofpe,dofpe);
               %En el peso de gauss ya viene multiplicado el espesor.
               m_pesoPG = m_DetJT(:,iElem).*wg;
               for iPG = 1:nPG
                  switch conshyp
                     case {1,2,4,5,10,11,12,52}  %Se considera todos como modelos constitutivos elásticos
                        ct = e_DatMatSet.ce;
                     case {50,51}  %Modelo multiescala clásico y cohesivo
                        e_VGMicro = e_DatMatSet.e_VG;
                        e_DatSetMicro = e_DatMatSet.e_DatSet;
                        nElemMicro = e_VGMicro.nElem;
                        [c_CTMicro,KTMicro] = f_TensTangElastCU(e_DatSetMicro,e_VGMicro);
                        [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro,...
                           e_DatMatSet.xx,e_DatSetMicro,e_VGMicro.m_ConecFront);
                        ct = f_ModTangHomog(KTMicro,c_CTMicro,m_LinCondMicro,doflMicro,doffMicro,...
                           e_DatSetMicro,e_DatMatSet.Omega_micro,true(nElemMicro,1),true(nElemMicro,1),...
                           e_VGMicro);
                     otherwise
                        error(['Tensores tangente constitutivos elásticos: Modelo constitutivo no ',...
                           'implementado.'])
                  end
                  m_CT(:,:,iPG,iElem) = ct;
                  B = m_BT(:,:,iPG,iElem);
                  m_kt = m_kt+B'*ct*B*m_pesoPG(iPG);
               end
               
            %case 10
               %Ver que la mayoría de los elementos no influyen en la forma de cálculo del tensor constitutivo
               %excepto por ejemplo el SDA.
               
            otherwise
               error(['Tensores tangentes constitutivos elásticos: Elemento finito no ',...
                           'implementado.'])
         end
         
         m_Ke(:,:,iElem) = m_kt;
         
      end
      
      c_CT{iSet} = m_CT;
      c_Ke{iSet} = m_Ke(:);
      c_Fil{iSet} = m_Fil;
      c_Col{iSet} = m_Col;
      
   end
   
   % Ensamble de matriz de rigidez global
   KT = sparse([c_Fil{:}],[c_Col{:}],cat(1,c_Ke{:})); 


end