function [m_CTHomog,m_SigmaHomog,hvar_newMacro,vectVHElem_new,c_Tracc] =...
   f_RMap_MEStafe(...
   m_IDefMacroReg,m_IDefMacroLoc,hvar_oldMacro,...
   e_DatMatSetMacro,condBif,ind_ActState_old,...
   e_VGMacro,vectVHElem_old,ksd,kinf)

%Se recupera variables micro
xx            = e_DatMatSetMacro.xx;

omegaMicro    = e_DatMatSetMacro.omegaMicro;
e_DatSet      = e_DatMatSetMacro.e_DatSet;
e_VGMicro     = e_DatMatSetMacro.e_VG;
esImplexMacro = e_DatMatSetMacro.esImplex;
e_VarEst_old  = hvar_oldMacro.e_VarEst;
u             = hvar_oldMacro.u;
c_GdlCond     = hvar_oldMacro.c_GdlCond;
Fint          = hvar_oldMacro.Fint;
m_LinCond     = hvar_oldMacro.m_LinCond;
doff          = hvar_oldMacro.doff;
dofl          = hvar_oldMacro.dofl;
e_VarAux      = hvar_oldMacro.e_VarAux;
c_DefMacro    = hvar_oldMacro.c_DefMacro;

omegaMicroL    = hvar_oldMacro.omegaMicroL;
m_ElemLoc      = hvar_oldMacro.m_ElemLoc;
lMacro         = hvar_oldMacro.lMacro;
lMicro         = hvar_oldMacro.lMicro;
c_NormalesMicro = hvar_oldMacro.c_NormalesMicro;
longFis        = hvar_oldMacro.longFis;
facNormMicro   = hvar_oldMacro.facNormMicro;


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
e_VGMicro.SharedFS = e_VGMacro.SharedFS;

%Factor de regularización para la etapa intermedia (integración en un solo PG).
if ind_ActState_old==1
   %Se asume que en kinf viene el espesor de regularización que depende del espesor del elemento.
   hRVE = e_DatMatSetMacro.tamCaracMC;
   facReg = ksd/hRVE;
else
   facReg = 1;
end

% DEFORMACION MACRO TOTAL A APLICAR EN CADA PG
c_DefMacro = cellfun(@(x)bsxfun(@plus,x,facReg*m_IDefMacroReg),c_DefMacro,'UniformOutput',false);
if ind_ActState_old>1&&e_VGMacro.iPG==6
   %Se asume que IDefMacroLocPG es nula previo a la bifurcación y a la imposición de la SD.
   for iSet = 1:e_VGMicro.nSet
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      eltype = e_DatElemSet.eltype;
      switch eltype
         case 32
            nElem = e_DatSet(iSet).nElem;
            nPG = e_DatSet(iSet).e_DatElem.npg;
            n_Micro = e_DatSet(iSet).e_DatElem.normal_micro;
            m_SentNorm = e_VarAux(iSet).VarAuxElem(1,:);
            %                 m_lMicro = e_DatSet(iSet).e_DatElem.le_Elem;
            ksb = e_DatSet(iSet).e_DatElem.ksb;
            m_ElemLocSet = m_ElemLoc(e_DatSet(iSet).m_IndElemSet);
            m_DefMacro = c_DefMacro{iSet};
            for iElem = 1:nElem
               if m_ElemLocSet(iElem)
                  m_nMicroi = m_SentNorm(iElem)*n_Micro(:,iElem);
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  for iPG = 1:nPG
                     Def_sing = [m_nMicroi(1)*m_IDefMacroLoc(1);...
                                 m_nMicroi(2)*m_IDefMacroLoc(2);...
                                 0;...
                                 m_nMicroi(1)*m_IDefMacroLoc(2)+m_nMicroi(2)*m_IDefMacroLoc(1)];
                     %Se asume que las normales micro que no pertenecen al dominio localizado son nulas.
                     m_DefMacro(:,iPG,iElem) = m_DefMacro(:,iPG,iElem)+...
                        Def_sing/ksb(iElem)/facNormMicro;
                  end
               end
            end
            c_DefMacro{iSet} = m_DefMacro;
      end
   end
end

% ESQUEMA DE NEWTON-RAPHSON  #################################
%ticIDNewMicr = tic;
%fprintf('Tiempo del Newton micro: %f\n',toc(ticIDNewMicr));
Du_step_new = zeros(ndoft,1);
[u,c_GdlCond,DFint,e_VarEst_new,e_VarAux,Du_step_new,c_CT,KT] = newton_raphson_MICRO(...
   xx,m_LinCond,dofl,doff,u,Du_step_new,c_GdlCond,Du_step_old,Fint,Fext,e_VarEst_old,e_VarAux,e_DatSet,...
   c_DefMacro,e_VGMicro);
%Se guarda las fuerzas internas sólo para impresión de las gráficas X-Y
Fint = Fint+DFint;

%Para considerar que en el caso del cuadrángulo Q1 se debe integrar las tensiones estabilizadas.
c_Tens_old = cell(nSet,1);
c_Tens_new = cell(nSet,1);
for iSet = 1:nSet
   e_DatElemSet = e_DatSet(iSet).e_DatElem;
   eltype = e_DatElemSet.eltype;
   switch eltype
      case {2,4,8,31,32}
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

%% Solo si el PG == 6
if e_VGMacro.iPG==6
   vectVHElem_new = Regularization_parameters_MICRO(xx,e_VarEst_new,...
      e_VarEst_old,e_VarAux,Du_step_new,e_VGMicro,e_DatSet,vectVHElem_old,...
      kinf,omegaMicro,ind_ActState_old);
else
   vectVHElem_new = [];
end

% INCREMENTO DE TENSIONES HOMEGENEIZADAS
c_Tracc   = zeros(nDime,1);
if e_VGMacro.iPG==6&&ind_ActState_old>1
   %Esta matriz de los determinante de J con valores no nulos dentro del dominio omegaL podría
   %precalcularse si fuese necesario.
   %      c_DetJTLoc = arrayfun(@(x)bsxfun(@times,x.m_DetJT,m_ElemLoc(x.m_IndElemSet)),e_DatSet,...
   %         'UniformOutput',false);
   %Es necesario realizar esta resta en todos los puntos de gauss de los elementos micro porque justo en el
   %paso en que se empieza a homogeneizar en omegaL, la tensión previa homogenizada m_sigmaHomog_old, que
   %fue obtenida en una homogeneización en todo omegaMicro, puede ser distinta a la homogenización de las
   %tensiones sigma_old en omegaL (al ser dominios distintos). Se debe interpretar que todas las tensiones
   %previas de cada PG de la estructura micro son las tensiones previas del problema macro y no sus valores
   %homogeneizados.
   %c_Dsigma = cellfun(@(x,y)minus(x,y),c_Tens_new,c_Tens_old,'UniformOutput',false);
   
   m_VarFluc = f_VarFluct(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,m_ElemLoc,e_VGMicro);
   m_CTHomog      = zeros(ntens,ntens);
   m_CTHomog_bb   = zeros(nDime,nDime);
   m_CTHomog_ub   = zeros(nDime,ntens);
   m_IndTens = eye(ntens);
   %%%%%%%%%%%%%%%%%%%
   nSet=e_VGMicro.nSet;
   %OMEG=0; LF=0;
   for iSet = 1:nSet
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      eltype       = e_DatElemSet.eltype;
      
      switch eltype
         case 32
            nElem           = e_DatSet(iSet).nElem;
            nPG             = e_DatSet(iSet).e_DatElem.npg;
            m_Sigma         = reshape(c_Tens_new{iSet},ntens,nPG,nElem);
            %m_Sigma_old     = reshape(c_Tens_old{iSet},ntens,nPG,nElem);
            m_DofElem       = e_DatSet(iSet).m_DofElem;
            %m_DetJT         = e_DatSet(iSet).m_DetJT;
            m_BT            = e_DatSet(iSet).m_BT;
            m_CT            = c_CT{iSet};
            %wg              = e_DatElemSet.wg;
            m_ElemLocSet    = m_ElemLoc(e_DatSet(iSet).m_IndElemSet);
            n_Micro = e_DatElemSet.normal_micro;
            m_SentNorm = e_VarAux(iSet).VarAuxElem(1,:);
            %m_lMicro        = e_DatElemSet.le_Elem;
            ksb = e_DatSet(iSet).e_DatElem.ksb;
            VolElem         = e_DatSet(iSet).m_VolElem;
            
            for iElem = 1:nElem
               if m_ElemLocSet(iElem)
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  dofElem       = m_DofElem(:,iElem);
                  %m_pesoPG      = m_DetJT(:,iElem).*wg;
                  m_VarFlucElem = m_VarFluc(dofElem,:);
                  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  m_nMicroi = m_SentNorm(iElem)*n_Micro(:,iElem);
                  m_NormalesMicro = [m_nMicroi(1),0;...
                                     0,m_nMicroi(2);...
                                     0,0;...
                                     m_nMicroi(2),m_nMicroi(1)];
                  
                  %En forma estricta acá está mal, ya que al usar VolElem se está asumiendo que tiene un solo
                  %PG. Hay que usar m_pesoPG por si llegara a integrar con más PG el triángulo lineal.
                  for iPG = 1:nPG
                     
                     c_Tracc = c_Tracc+(m_NormalesMicro'*m_Sigma(:,iPG,iElem))/ksb(iElem)*VolElem(iElem);
                     
                     pp = m_NormalesMicro'*m_CT(:,:,iPG,iElem)*(m_IndTens+m_BT(:,:,iPG,iElem)*m_VarFlucElem);
                     
%                      m_CTHomog_bb = m_CTHomog_bb+...
%                         (pp*m_NormalesMicro/facNormMicro/ksb(iElem))/ksb(iElem)*VolElem(iElem);
                     m_CTHomog_bb = m_CTHomog_bb+(pp*m_NormalesMicro/ksb(iElem))/ksb(iElem)*VolElem(iElem);
                     %Usar el mismo m_VarFlucElem, calculado a partir de la pertubación de la deformación
                     %macro solo en el dominio localizado, para el término ub solo es válido si en este
                     %dominio no se produce fluctuaciones (es decir si está restringido con condiciones de
                     %borde tipo Taylor).
                     m_CTHomog_ub  = m_CTHomog_ub + pp/ksb(iElem)*VolElem(iElem) ;

                  end
               end
            end
      end
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %               c_Tracc{iSet} = m_Tracc;
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   end
   
   %%%%%%%%%%%%%%%%%%%%%%%%%
   %Para que se calcula el sigma homogeneizado para el PG6, para llevar su historia??
   m_SigmaHomog = f_HomogArea(c_Tens_new,ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VGMicro);
   
   c_Tracc   = c_Tracc /(facNormMicro*longFis);
   m_CTHomog_bb = m_CTHomog_bb /(facNormMicro*longFis);
   m_CTHomog_ub = m_CTHomog_ub /(facNormMicro*longFis);
   %%    m_CTHomog_ub = m_CTHomog_ub /(facNormMicro);
   m_CTHomog(1:nDime,1:nDime)   = m_CTHomog_bb;
   m_CTHomog(nDime+1:2*nDime,:) = m_CTHomog_ub;
   %%%%%%%%%%%%%%%%%%%%%%%%%
   %      m_SigmaHomog = f_HomogArea(c_Tracc,nDime,facNormMicro*longFis,c_DetJTLoc,e_DatSet,e_VG);
   %      m_CTHomog = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicroL,...
   %          m_ElemLoc,m_ElemLoc,e_VG);
   %%%%%%%%%%%%%%%%%%%%%%%%%
   
   
   %Cuando el modelo constitutivo multiescala bifurca ya no se calcula más tensor tangente implícito, pero
   %igual se debe devolver una estructura, con cualquier valor para que el ensamblaje sea correcto en las
   %matrices m_TensorTang (ver como plantear esto de otra forma). No se puede descartar la parte que guarda
   %los tensores implícitos ya que puede haber PG que ya bifurcaron y otros que no.
   if esImplexMacro
      m_CTHomog = struct('Implex',m_CTHomog,'Impli',zeros(ntens,ntens));
   end
else
   %Se asume que previo a la bifurcación y al imponer la discontinuidad fuerte, las tensiones old
   %homogeneizadas son iguales a las que se obtiene con la homogeneización de las sigma_old micro
   %de la celda unitaria (esto se realiza para ahorrarse la resta de las todas las tensiones
   %micro y homogeneizar el delta sigma).
   %       m_DsigmaHomog = m_sigmaHomog-m_sigmaHomog_old;
   %    c_Dsigma = cellfun(@(x,y)minus(x,y),c_Tens_new,c_Tens_old,'UniformOutput',false);
   m_SigmaHomog = f_HomogArea(c_Tens_new,ntens,omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VGMicro);
   % Tensor Tangente Homogeneizado
   %       m_CTHomog = f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,e_VG);
   m_CTHomog = facReg*f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,...
      true(nElem,1),true(nElem,1),e_VGMicro);
   %IMPLEX: Solo interesa que homogenice los tensores tangentes constitutivos implícitos hasta
   %la bifurcación. Luego no se utiliza más, ya que la mismas se utilizan en el análisis de bifurcación.
   if esImplexMacro
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
      %Cálculo de la matriz implícita
      KT = f_MatGlobal(c_CT,e_VarAux,e_DatSet,e_VGMicro);
      %Homogenización del tensor implícito
      m_CTHomogImpli = facReg*f_ModTangHomog(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,omegaMicro,...
         true(nElem,1),true(nElem,1),e_VGMicro);
      m_CTHomog = struct('Implex',m_CTHomog,'Impli',m_CTHomogImpli);
   end
end

hvar_newMacro = struct('u',u,'c_GdlCond',{c_GdlCond},'Fint',Fint,'e_VarEst',e_VarEst_new,...
   'e_VarAux',e_VarAux,'m_LinCond',m_LinCond,'doff',doff,'dofl',dofl,'m_ElemLoc',m_ElemLoc,...
   'c_DefMacro',{c_DefMacro},'omegaMicroL',omegaMicroL,'lMacro',lMacro,'lMicro',lMicro,...
   'c_NormalesMicro',{c_NormalesMicro},'longFis',longFis,'facNormMicro',facNormMicro);     %,'m_TensProy',m_TensProy

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
      case {2,4,8,20,31,32}
         m_pesoPG = bsxfun(@times,m_DetJT,wg);
         %Con las matrices de deformación precalculadas la determinación de las matrices elementales
         %son iguales para los elementos estándares.
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