function e_VarAux = f_eVarAuxInic(xx,e_DatSet,e_VG)

   %Las variables auxiliares son variables que se utilizan para pasar valores entre la funciones y
   %que pueden variar en cada iteración. Distinto a las variables históricas que se mantiene dos
   %copias de las mismas (old y new), en esta se mantiene el valor más nuevo de las variables,
   %pudiendo ser modificadas por funciones previas.
   %El mismo efecto se puede conseguir pasando como argumento las variables new, y no
   %inicializándolas, pero la desventaja es que si las variables no necesitan de la mecánica de la
   %variable histórica (new y old) se está duplicando resultados innecesariamente.
   
   nSet = e_VG.nSet;
   
   e_VarAux(nSet,1) = struct('VarAuxGP',[],'VarAuxElem',[]);
   for iSet = 1:nSet
      
      nElem = e_DatSet(iSet).nElem;
      siavare = e_DatSet(iSet).siavare;
      e_DatMatSet = e_DatSet(iSet).e_DatMat;
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      nVarAuxElem = e_DatElemSet.nVarAuxElem;
      conshyp = e_DatMatSet.conshyp;
      eltype = e_DatElemSet.eltype;
      
      %Variables auxiliares por punto de gauss
      switch conshyp
         case {1,2,3,4,5,6,7,8,9,10,50,52,53,100,110}
            e_VarAux(iSet).VarAuxGP = zeros(siavare,nElem);
         case {11,12,13}
            %siavarpg = e_DatMatSet.siavarpg;
            %esImplex = e_DatMatSet.esImplex;
            e_VarAux(iSet).VarAuxGP = zeros(siavare,nElem);
%             if esImplex
%                %El factor del incremento debe inicializarse en 1.
%                e_VarAux(iSet).VarAuxGP(3:siavarpg:siavare,:) = 1;
%             end
         case {51,54,55}
            %Como es el mismo set, se asume que las celdas unitarias de los todos los puntos de 
            %gauss de todos los elementos tiene las mismas condiciones de borde iniciales.
            %e_VGMicro = e_DatMatSet.e_VG;
            e_DatSetMicro = e_DatMatSet.e_DatSet;
%             [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro,...
%                e_DatMatSet.xx,e_DatSetMicro,e_VGMicro.m_ConecFront);
%             e_VarAux(iSet).VarAuxGP(1:siavare,1:nElem) = struct('m_LinCond',m_LinCondMicro,...
%                'doff',doffMicro,'dofl',doflMicro);
            e_VarAuxGP(1:siavare,1:nElem) = struct('m_ElemLocCalc',[],'c_ProyIncrEps',...
               {arrayfun(@(x)zeros(1,x.nElem),e_DatSetMicro,'UniformOutput',false)});
            e_VarAux(iSet).VarAuxGP = e_VarAuxGP;
         otherwise
            error('Inicialización variables auxiliares por PG: Modelo constitutivo no definido.')            
      end
      
      %Variables auxiliares por elemento
      switch eltype
         case {2,4,8,5,7,10,20,21,22,23,31}
            e_VarAux(iSet).VarAuxElem = zeros(nVarAuxElem,nElem);
         case 32
            e_VarAux(iSet).VarAuxElem = zeros(nVarAuxElem,nElem);
            %Se inicializa el sentido de las normales micro (la dirección se lleva en el e_DatSet del
            %elemento)
            e_VarAux(iSet).VarAuxElem(1,:) = 1;
      %
      %LARGE DEFORMATIONS
         case 108  %Cuadrángulo de 4 nodos FBar.
            e_VarAux(iSet).VarAuxElem = zeros(nVarAuxElem,nElem);
            %Se almacena la matriz de deformación del punto central para cálculo del gradiente de deformación
            %del punto central. Ver si no hacerlo de otra forma.
            m_Conec = e_DatSet(iSet).conec;
            xPG0 = [0,0];
            e_DatElemSetAux = e_DatElemSet;
            e_DatElemSetAux.xg = xPG0;
            e_DatElemSetAux.npg = 1;
            for iElem = 1:nElem
               m_CoordElem = f_CoordElem(xx,m_Conec(iElem,:));
               m_BTAux = squeeze(f_MatBe_quad_q1(m_CoordElem,e_DatElemSetAux,e_VG));
               e_VarAux(iSet).VarAuxElem(:,iElem) = m_BTAux(:);
            end
         otherwise
            error('Inicialización variables auxiliares por Elemento Finito: Tipo de elemento no definido.')
      end
      
   end

end