function m_Tracc = f_HomogTracc(e_DatSet,e_VarEst,e_VarAux,m_ElemLoc,facNormMicro,longFis,e_VG)

nDime = e_VG.ndime;
nTens = e_VG.ntens;
%Se asume que se está pasando los datos del PG (el 6 normalmente) sobre el que se quiere hacer el cálculo.
m_Tracc = zeros(nDime,1);
m_NormalMicro = zeros(nTens,nDime);


nSet = e_VG.nSet;
%
for iSet = 1:nSet
   
   e_DatElemSet = e_DatSet(iSet).e_DatElem;
   eltype = e_DatElemSet.eltype;
   nElem = e_DatSet(iSet).nElem;
   m_ElemLocSet = m_ElemLoc(e_DatSet(iSet).m_IndElemSet);
   %
   switch eltype
      case 32
         nPG = e_DatElemSet.npg;
         m_DetJT = e_DatSet(iSet).m_DetJT;
         wg = e_DatElemSet.wg;
         %
         n_Micro = e_DatElemSet.normal_micro;
         m_SentNorm = e_VarAux(iSet).VarAuxElem(1,:);
         ksb = e_DatElemSet.ksb;
         %Al ser lineal, y siempre se intregrase un solo PG, se podría usar volumen. 
         %VolElem = e_DatSet(iSet).m_VolElem;
         m_Tens = reshape(e_VarEst(iSet).sigma,nTens,nPG,nElem);
         %
         for iElem = 1:nElem
            if m_ElemLocSet(iElem)
               m_pesoPG = m_DetJT(:,iElem).*wg;
               m_nMicroi = m_SentNorm(iElem)*n_Micro(:,iElem);
               %m_NormalesMicro(1,1) = m_nMicroi(1);
               %m_NormalesMicro(2,2) = m_nMicroi(2);
               %m_NormalesMicro(4,1) = m_nMicroi(2);
               %m_NormalesMicro(4,2) = m_nMicroi(1);
               m_NormalMicro([1,8]) = m_nMicroi(1);
               m_NormalMicro([4,6]) = m_nMicroi(2);
               for iPG = 1:nPG
                  m_Tracc = m_Tracc+(m_NormalMicro'*m_Tens(:,iPG,iElem))/ksb(iElem)*m_pesoPG(iPG);                  
               end
            end
         end         
      otherwise
         if any(m_ElemLocSet)
            error('Tracción homogeneizada: Dominio Localizado: Tipo de elemento finito no definido.')
         end
   end
   
end

%m_SigmaHomog = f_HomogArea(c_Tens_new,ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VGMicro);
%m_Tracc = f_HomogArea(c_Tracc,nDime,facNormMicro*longFis,c_DetJTLoc,e_DatSet,e_VG);
m_Tracc = m_Tracc/facNormMicro/longFis;

end