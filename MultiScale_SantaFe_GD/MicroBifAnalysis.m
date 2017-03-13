function [isBif,minQ] = MicroBifAnalysis(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,e_VGMicro)

omegaMicro = e_VGMicro.omegaMicro;
nElem      = e_VGMicro.nElem;
nSet       = e_VGMicro.nSet;

% Tensor Tangente Homogeneizado
%m_CTHomog = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,...
%    true(nElem,1),true(nElem,1),e_VGMicro);

% Tensor constitutivo homogenizado IMPLICITO
%IMPLEX: Para casos implex, es necesario tener en cuenta la organizacion de
%la estructura c_CT porque esta almacena tanto los tensores implex como los
%implicitos
for iSet = 1:nSet
    e_DatMatSet = e_DatSet(iSet).e_DatMat;
    nPG = e_DatSet(iSet).e_DatElem.npg;
    esImplex = e_DatMatSet.esImplex;
    if esImplex
        c_CT{iSet} = c_CT{iSet}(:,:,nPG+1:2*nPG,:);
    end
end
%Calculo de la matriz global implicita
KT = f_MatGlobal(c_CT,e_DatSet,e_VGMicro);

e_VGMicro.isMICRO.MICRO=0 ;  %JLM
%Homogenizaci�n del tensor implicito
m_CTHomogImpli = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,...
    true(nElem,1),true(nElem,1),e_VGMicro);

% Bifurcation analysis
[isBif,minQ] = BifurcMicro(m_CTHomogImpli,e_VGMicro) ;

end

function KT = f_MatGlobal(c_CT,e_DatSet,e_VG)

      %Esta funci�n devuelve la matriz de rigidez global si se conoce cu�les son los tensores constitutivos
      %tangentes de todos los puntos de gauss de la estructura.
      nSet = e_VG.nSet;
      %ndn = e_VG.ndn;
      %
      c_Ke = cell(nSet,1);
      c_Fil = cell(nSet,1);
      c_Col = cell(nSet,1);
      %
      for iSet = 1:nSet
         
         %Ver si no es peor usar esto, por ejemplo que haga alguna copia innecesaria, que usar siempre
         %e_DatSet(iSet).
         e_DatiSet = e_DatSet(iSet);
         e_DatElemSet = e_DatiSet.e_DatElem;
         nPG = e_DatElemSet.npg;
         wg = e_DatElemSet.wg;
         dofpe = e_DatElemSet.dofpe;
         eltype = e_DatElemSet.eltype;
         nElem = e_DatiSet.nElem;
         %conec = e_DatiSet.conec;
         m_DofElem = e_DatiSet.m_DofElem;
         m_BT = e_DatiSet.m_BT;
         m_DetJT = e_DatiSet.m_DetJT;
         m_CT = c_CT{iSet};
         %
         m_Ke = zeros(dofpe*dofpe,nElem);     
         %
         %dofElemSet = f_DofElem(reshape(conec',1,[]),ndn);
         dofElemSet = m_DofElem(:);
         m_Fil = reshape(repmat(reshape(dofElemSet,dofpe,[]),dofpe,1),1,[]);
         m_Col = reshape(repmat(dofElemSet',dofpe,1),1,[]); 
         %
         switch eltype
            case {2,4,8,20,31,32,108}
               m_pesoPG = bsxfun(@times,m_DetJT,wg);
               %Con las matrices de deformaci�n precalculadas la determinaci�n de las matrices elementales
               %son iguales para los elementos est�ndars.
               for iElem = 1:nElem
               %parfor iElem = 1:nElem
                  m_kt = zeros(dofpe,dofpe);
                  for iPG = 1:nPG
                     B = m_BT(:,:,iPG,iElem);
                     ct = m_CT(:,:,iPG,iElem);
                     m_kt = m_kt+B'*ct*B*m_pesoPG(iPG,iElem);
                  end
                  %Para acelerar un poco el c�lculo se guarda en forma de columna las matrices de rigidez
                  %elementales.
                  m_Ke(:,iElem) = m_kt(:);
               end   
            otherwise
               error(['Modelo Multiescala Cohesivo: Matriz de rigidez global impl�cita para el implex: ',...
                  'Tipo de elemento no definido'])
         end
         c_CT{iSet} = m_CT;
         c_Ke{iSet} = m_Ke(:);
         c_Fil{iSet} = m_Fil;
         c_Col{iSet} = m_Col;
         
      end
      
      % Ensamble de matriz de rigidez global
      KT = sparse([c_Fil{:}],[c_Col{:}],cat(1,c_Ke{:})); 
      
end