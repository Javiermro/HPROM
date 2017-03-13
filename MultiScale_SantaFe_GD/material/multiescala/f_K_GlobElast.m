function [K_globElast,BTCe_Glob] = f_K_GlobElast(e_DatSet,e_VG)

   %Esta función devuelve todos los tensores tangentes constitutivos de los puntos de gauss de la malla y la
   %matriz de rigidez global si tiene un comportamiento elástico, sin considerar si es daño o plasticidad. Es
   %decir se considera todos los modelos constitutivos, pero se considera solo sus propiedades elásticas. Esto
   %permite utilizar un archivo de datos de la microcelda para elementos macro que se quiere que bifurquen y
   %para los elmentos que se quiere se comporten siempre elásticamente.
   
   nSet = e_VG.nSet;
   ntens   = e_VG.ntens;

   %
   c_Ke    = cell(nSet,1);
   c_BTCe1 = cell(nSet,1);
   c_BTCe2 = cell(nSet,1);
   c_BTCe3 = cell(nSet,1);
   c_BTCe4 = cell(nSet,1);
   
   
   c_Fil = cell(nSet,1);
   c_Col = cell(nSet,1);
   c_FilFza = cell(nSet,1);
   
   K_globElast = sparse([]);
   BTCe_Glob   = sparse([]);
   
   for iSet = 1:nSet
       e_DatMatSet = e_DatSet(iSet).e_DatMat;
       conshyp     = e_DatMatSet.conshyp;
       if conshyp ==1
           
          ct = e_DatMatSet.ce;
           
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
          %
          m_Ke  = zeros(dofpe,dofpe,nElem);
          m_BTC1 = zeros(dofpe,nElem);
          m_BTC2 = zeros(dofpe,nElem);
          m_BTC3 = zeros(dofpe,nElem);
          m_BTC4 = zeros(dofpe,nElem);
          % Grados de libertad y coordenadas de los nodos de los elementos del set
          %dofElemSet = f_DofElem(reshape(conec',1,[]),ndn);
          dofElemSet = m_DofElem(:);
          m_Fil = reshape(repmat(reshape(dofElemSet,dofpe,[]),dofpe,1),1,[]);
          m_Col = reshape(repmat(dofElemSet',dofpe,1),1,[]);
          m_FilFza = dofElemSet';

          %parfor iElem = 1:nElem
          for iElem = 1:nElem
              
              switch eltype
                  case {2,4,8,31,32}
                      
                      m_kt   = zeros(dofpe,dofpe);
                      m_BTCe = zeros(dofpe,ntens);
                      %En el peso de gauss ya viene multiplicado el espesor.
                      m_pesoPG = m_DetJT(:,iElem).*wg;
                      for iPG = 1:nPG
                          B      = m_BT(:,:,iPG,iElem);
                          m_kt   = m_kt+B'*ct*B*m_pesoPG(iPG);
                          m_BTCe = m_BTCe+B'*ct*m_pesoPG(iPG);
                      end
                      
                  otherwise
                      error(['Matriz Global Elástica: Elemento finito no ',...
                          'implementado.'])
              end
              
              m_Ke(:,:,iElem)  = m_kt;
              m_BTC1(:,iElem) = m_BTCe(:,1);
              m_BTC2(:,iElem) = m_BTCe(:,2);
              m_BTC3(:,iElem) = m_BTCe(:,3);
              m_BTC4(:,iElem) = m_BTCe(:,4);
              
          end
          
          c_Ke{iSet}     = m_Ke(:);
          c_BTCe1{iSet}   = m_BTC1(:);
          c_BTCe2{iSet}   = m_BTC2(:);
          c_BTCe3{iSet}   = m_BTC3(:);
          c_BTCe4{iSet}   = m_BTC4(:);
          c_Fil{iSet}    = m_Fil;
          c_Col{iSet}    = m_Col;
          c_FilFza{iSet} = m_FilFza;
          
      end
   end
   
   % Ensamble de matriz de rigidez global
   K_globElast = sparse([c_Fil{:}],[c_Col{:}],cat(1,c_Ke{:})); 
   BTCe_Glob(:,1)   = sparse([c_FilFza{:}],1,cat(1,c_BTCe1{:}));
   BTCe_Glob(:,2)   = sparse([c_FilFza{:}],1,cat(1,c_BTCe2{:}));
   BTCe_Glob(:,3)   = sparse([c_FilFza{:}],1,cat(1,c_BTCe3{:}));
   BTCe_Glob(:,4)   = sparse([c_FilFza{:}],1,cat(1,c_BTCe4{:}));
   
end