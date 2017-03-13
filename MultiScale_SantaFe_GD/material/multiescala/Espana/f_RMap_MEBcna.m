function [m_CTHomog,m_SigmaHomog,hvar_newMacro,vectVHElem_new] = f_RMap_MEBcna(...
    m_IDefMacroReg,hvar_oldMacro,...
    e_DatMatSetMacro,condBif,e_VGMacro,vectVHElem_old,kinf)

    %Se recupera variables micro
    xx            = e_DatMatSetMacro.xx;
    u             = hvar_oldMacro.u;
    omegaMicro    = e_DatMatSetMacro.omegaMicro;
    e_DatSet      = e_DatMatSetMacro.e_DatSet;
    e_VGMicro     = e_DatMatSetMacro.e_VG;
    esImplexMacro = e_DatMatSetMacro.esImplex;
    e_VarEst_old  = hvar_oldMacro.e_VarEst;
    c_GdlCond     = hvar_oldMacro.c_GdlCond;
    Fint          = hvar_oldMacro.Fint;
    m_LinCond     = hvar_oldMacro.m_LinCond;
    doff          = hvar_oldMacro.doff;
    dofl          = hvar_oldMacro.dofl;
    e_VarAux      = hvar_oldMacro.e_VarAux;
    c_DefMacro    = hvar_oldMacro.c_DefMacro;
    % VARIABLES GLOBALES
    nElem        = e_VGMicro.nElem;
    ntens        = e_VGMicro.ntens;
    ndoft        = e_VGMicro.ndoft;
    nSet         = e_VGMicro.nSet;
    nDime        = e_VGMicro.ndime;
    % INICIALIZACION DE VARIABLES
    Du_step_old             = zeros(ndoft,1);
    Fext                    = zeros(ndoft,1);
    e_VGMicro.elast         = e_VGMacro.elast;
    %Por si es necesario el paso tiempo a nivel micro.
    e_VGMicro.istep         = e_VGMacro.istep;
    e_VGMicro.iterMacro     = e_VGMacro.iter;
    %Para imprensión se guarda el número de elemento macro y número de PG macro.
    e_VGMicro.iElemNumMacro = e_VGMacro.iElemNum;
    e_VGMicro.iPGMacro      = e_VGMacro.iPG;
    %Nombre interno de la celda unitaria
    %(hacer esto en cada iteración puede ser medio lento, ver que hacer sino).
    %e_VGMicro.fileCompleto = [e_VG.fileCompleto,'_EM',int2str(e_VGMacro.iElemNum),'PGM',int2str(e_VGMacro.iPG)];
    %Se guarda los datos del matlabpool (se utiliza para imprimir de pasos y iteraciones a nivel micro que no
    %convergieron).
    e_VGMicro.nLab = e_VGMacro.nLab;
    e_VGMicro.tipoMPool = e_VGMacro.tipoMPool;

    % DEFORMACION MACRO TOTAL A APLICAR EN CADA PG
    c_DefMacro = arrayfun(@(x)repmat(m_IDefMacroReg,[1,x.e_DatElem.npg,x.nElem]),...
        e_DatSet,'UniformOutput',false);

    % ESQUEMA DE NEWTON-RAPHSON  #################################
    Du_step_new = zeros(ndoft,1);
    [u,c_GdlCond,Fint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT] = newton_raphson_MICRO(...
        xx,m_LinCond,dofl,doff,u,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,e_VarEst_old,e_VarAux,e_DatSet,...
        c_DefMacro,e_VGMicro);
    %Se guarda las fuerzas internas sólo para impresión de las gráficas X-Y
   %% Fint = Fint+DFint;

    %Para considerar que en el caso del cuadrángulo Q1 se
    %debe integrar las tensiones estabilizadas.
    c_Tens_old = cell(nSet,1);
    c_Tens_new = cell(nSet,1);
    for iSet = 1:nSet
        e_DatElemSet = e_DatSet(iSet).e_DatElem;
        eltype = e_DatElemSet.eltype;
        switch eltype
            case {2,4,8,31}
                c_Tens_old{iSet} = e_VarEst_old(iSet).sigma;
                c_Tens_new{iSet} = e_VarEst_new(iSet).sigma;
            case 20
                c_Tens_old{iSet} = e_VarEst_old(iSet).VarHistElem;
                c_Tens_new{iSet} = e_VarEst_new(iSet).VarHistElem;
            otherwise
                error(['Modelo Multiescala: Homogeneización: Tensiones de homogeneización: Modelo ',...
                    'constitutivo no definido.'])
        end
    end

m_SigmaHomog = f_HomogArea(c_Tens_new,ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VGMicro);
% Obtencion de los parametros de regularizacion y variables homogenizadas de la microescala

%% Solo si el PG == 6
if e_VGMacro.iPG==6
    [vectVHElem_new] = Regularization_parameters_MICRO(xx,e_VarEst_new,...
        e_VarEst_old,e_VarAux,Du_step_new,e_VGMicro,e_DatSet,vectVHElem_old,...
        kinf,omegaMicro);
else
    vectVHElem_new=[];
end

% Tensor Tangente Homogeneizado
m_CTHomog = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,...
    true(nElem,1),true(nElem,1),e_VGMicro);


%IMPLEX: Solo interesa que homogenice los tensores tangentes constitutivos implícitos hasta
%la bifurcación. Luego no se utiliza más, ya que la mismas se utilizan en el análisis de bifurcación.
if esImplexMacro 
    if e_VGMacro.iPG==6
        %Se convierte c_CT con las matrices implícitas micro únicamente. Ver si no poner condicional implex
        %dentro de la función f_ModTangHomog en su lugar (también había que considerarla dentro de la función
        %f_MatGlobal.
        for iSet = 1:nSet
            e_DatMatSet = e_DatSet(iSet).e_DatMat;
            nPG = e_DatSet(iSet).e_DatElem.npg;
            esImplex = e_DatMatSet.esImplex;
            if esImplex
                c_CT{iSet} = c_CT{iSet}(:,:,nPG+1:2*nPG,:);
            end
        end
        %Cálculo de la matriz global implícita
        KT = f_MatGlobal(c_CT,e_VarAux,e_DatSet,e_VGMicro);
        %Homogenización del tensor implícito
        m_CTHomogImpli = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,...
            true(nElem,1),true(nElem,1),e_VGMicro);
    else
        m_CTHomogImpli = zeros(e_VGMacro.ntens);
    end
    
     m_CTHomog = struct('Implex',m_CTHomog,'Impli',m_CTHomogImpli);
end

hvar_newMacro = struct('u',u,'c_GdlCond',{c_GdlCond},'Fint',Fint,'e_VarEst',e_VarEst_new,...
    'e_VarAux',e_VarAux,'m_LinCond',m_LinCond,'doff',doff,'dofl',dofl,...
    'c_DefMacro',{c_DefMacro});     %,'m_TensProy',m_TensProy

end

function KT = f_MatGlobal(c_CT,e_VarAux,e_DatSet,e_VG)

      %Esta función devuelve la matriz de rigidez global si se conoce cuáles son los tensores constitutivos
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
            case {2,4,8,20,31}
               m_pesoPG = bsxfun(@times,m_DetJT,wg);
               %Con las matrices de deformación precalculadas la determinación de las matrices elementales
               %son iguales para los elementos estándars.
               for iElem = 1:nElem
               %parfor iElem = 1:nElem
                  m_kt = zeros(dofpe,dofpe);
                  for iPG = 1:nPG
                     B = m_BT(:,:,iPG,iElem);
                     ct = m_CT(:,:,iPG,iElem);
                     m_kt = m_kt+B'*ct*B*m_pesoPG(iPG,iElem);
                  end
                  %Para acelerar un poco el cálculo se guarda en forma de columna las matrices de rigidez
                  %elementales.
                  m_Ke(:,iElem) = m_kt(:);
               end   
%             case 20
%                %El término constante de este elemento, proveniente del problema mixto, se integra con 1 solo
%                %punto de gauss, en el lugar de 4. Por ello al PG 5 se incorpora el peso 4, que corresponde al
%                %caso de un solo punto de gauss (wg viene con cero para realizar correctamente otras
%                %integraciones dentro del código).
%                %wg(5) = 4;
%                m_pesoPG = bsxfun(@times,m_DetJT,wg);
%                % Variable auxiliar del elemento
%                %m_VarAuxElem = e_VarAux(iSet).VarAuxElem;
%                for iElem = 1:nElem                  
%                %parfor iElem = 1:nElem
%                   %% VERSIÓN 1                  
%                   %Condición de bifurcación
%                   %Se determina mediante el análisis de bifurcación del PG central (5), una vez que este
%                   %bifurcó todo el elemento se considera bifurcado (y en consecuencia los restantes PGs).
%                   condBif = m_VarAuxElem(1,iElem);                  
%                   % Parámetros de estabilización por PG
%                   %Se multiplica los parámetros de estabilizacón por los pesos de los puntos de gauss directamente, así tener
%                   %en cuenta este parámetro en la integración. 
%                   if condBif>0
%                      %Si condBif es igual a 1, significa que se detectó bifurcación.
%                      %Ver si no vale la pena activar la subintegración unos pasos después que bifurcó (condBif=2), o eso se
%                      %deja para la SD.
%                      m_pesoPGElem(1:4) = e_DatElemSet.estabBif*m_pesoPG(1:4,iElem);
%                      %Se asume que el PG central, donde se subintegra, está en la posición 5 de la lista de puntos de gauss.
%                      m_pesoPGElem(5) = (1-e_DatElemSet.estabBif)*m_pesoPG(5,iElem); 
%                   elseif condBif==0
%                      %Si condBif es cero, significa que todavía no se detectó la bifurcación.
%                      m_pesoPGElem(1:4) = e_DatElemSet.estabNoBif*m_pesoPG(1:4,iElem);
%                      %Se asume que el PG central, donde se subintegra, está en la posición 5 de la lista de puntos de gauss.
%                      m_pesoPGElem(5) = (1-e_DatElemSet.estabNoBif)*m_pesoPG(5,iElem); 
%                   end  
%                   m_kt = zeros(dofpe,dofpe);
%                   for iPG = 1:nPG
%                      B = m_BT(:,:,iPG,iElem);
%                      ct = m_CT(:,:,iPG,iElem);
%                      m_kt = m_kt+B'*ct*B*m_pesoPGElem(iPG);
%                   end
%                   %Para acelerar un poco el cálculo se guarda en forma de columna las matrices de rigidez
%                   %elementales.
%                   m_Ke(:,iElem) = m_kt(:);
%                   %% VERSIÓN 2
%                   %Condición de bifurcación
%                   %Se determina mediante el análisis de bifurcación del PG central (5), una vez que este
%                   %bifurcó todo el elemento se considera bifurcado (y en consecuencia los restantes PGs).
%                   condBif = m_VarAuxElem(1,iElem);    
%                   % Parámetros de estabilización por PG
%                   %Se multiplica los parámetros de estabilizacón por los pesos de los puntos de gauss
%                   %directamente, así tener en cuenta este parámetro en la integración.
%                   if condBif>0
%                      %Si condBif es igual a 1, significa que se detectó bifurcación.
%                      %Ver si no vale la pena activar la subintegración unos pasos después que bifurcó (condBif=2), o eso se
%                      %deja para la SD.
%                      paramEstab = e_DatElemSet.estabBif;
%                   elseif condBif==0
%                      %Si condBif es cero, significa que todavía no se detectó la bifurcación, o que descargó.
%                      paramEstab = e_DatElemSet.estabNoBif;
%                   end
%                   paramEstabRest = 1-paramEstab;
%                   B5 = m_BT(:,:,5,iElem); 
%                   ct5Estab = paramEstabRest*m_CT(:,:,5,iElem);
%                   m_kt = zeros(dofpe,dofpe);
%                   for iPG = 1:nPG-1
%                      B = m_BT(:,:,iPG,iElem);
%                      ctEstab = paramEstab*m_CT(:,:,iPG,iElem);
%                      m_kt = m_kt+(B'*ctEstab*B+B5'*ct5Estab*B5)*m_pesoPG(iPG,iElem);
%                   end
%                   %Para acelerar un poco el cálculo se guarda en forma de columna las matrices de rigidez
%                   %elementales.
%                   m_Ke(:,iElem) = m_kt(:);
%                   %% VERSIÓN 3
%                   %Se asume que las matrices m_CT son las implícitas estabilizadas, por lo que se integra en
%                   %el elemento con cuatro punto de gauss (se podría integrar en los 5 puntos de gauss ya que
%                   %el punto de gauss 5 tiene peso nulo).
%                   m_kt = zeros(dofpe,dofpe);
%                   for iPG = 1:nPG-1
%                      B = m_BT(:,:,iPG,iElem);
%                      ct = m_CT(:,:,iPG,iElem);
%                      m_kt = m_kt+B'*ct*B*m_pesoPG(iPG,iElem);
%                   end
%                   %Para acelerar un poco el cálculo se guarda en forma de columna las matrices de rigidez
%                   %elementales.
%                   m_Ke(:,iElem) = m_kt(:);
%                end   
            otherwise
               error(['Modelo Multiescala Cohesivo: Matriz de rigidez global implícita para el implex: ',...
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