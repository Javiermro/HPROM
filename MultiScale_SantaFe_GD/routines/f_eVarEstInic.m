function e_VarEst = f_eVarEstInic(c_Field,e_DatSet,nSet,MOD_TYPE,varargin)
   
%Inicializa la estructuras de variable de estado
   nField = length(c_Field);
   c_StructDef = cat(1,c_Field,repmat({[]},1,nField));
   %e_VarEst(1:nSet,1) = struct(c_StructDef{:});
   %Si son estructuras vac�as tambi�n se puede usar
   e_VarEst(nSet,1) = struct(c_StructDef{:});
   for iSet = 1:nSet
      sitvare = e_DatSet(iSet).sitvare;
      sihvare = e_DatSet(iSet).sihvare;
      nElem = e_DatSet(iSet).nElem;
      e_DatMatSet = e_DatSet(iSet).e_DatMat;
      e_DatElemSet = e_DatSet(iSet).e_DatElem;
      conshyp = e_DatMatSet.conshyp;
      nVarHistElem = e_DatElemSet.nVarHistElem;

      m_InicZerosSet = zeros(sitvare,nElem);
      for iField = 1:nField
         if strncmp(c_Field{iField},'hvar',4)
            switch conshyp
               case {1,2,3,4,5,6,7,8,9,10,11,12,13,52,100,110}
                  v_InicHvarSet = zeros(sihvare,nElem);
                  if conshyp==11
                     %Se inicializa la variables r_old(5), q_old(6), rInic_old(8) y qInic_old(9). Todas tienen
                     %que ser igual a r0.
                     sihvarpg = e_DatMatSet.sihvarpg;
                     if e_DatMatSet.tit==1
                        %En el caso de ser exponencial tambi�n se inicializa las rInic_old(end-1) y
                        %qInic_old(end).                        
                        m_Ind = [5;6;sihvarpg-1;sihvarpg];
                     else
                        m_Ind = [5;6];
                     end
                     v_InicHvarSet(bsxfun(@plus,m_Ind,0:sihvarpg:sihvare-1),:) = e_DatMatSet.r_0;                       
                  elseif conshyp==110
                      if isfield(e_DatMatSet,'sihvarpg') %JLM PD y HROM
                          sihvarpg = e_DatMatSet.sihvarpg;
                      else
                          sihvarpg = e_DatSet.sihvarpg;
                      end

                     m_Ind = [3;4;5];
                     v_InicHvarSet(bsxfun(@plus,m_Ind,0:sihvarpg:sihvare-1),:) = 1;                       
                  end 
               case {50,55}   %Modelo MultiEscala continuo
                   clear v_InicHvarSet
                   if nargin>4&&varargin{1}==0
                       %Se genera la estructura interna sin inicializar ninguna de las variables internas. Esto
                       %se utiliza en el newton para las variables new y evita el c�lculo de las matrices de
                       %condiciones de borde, que principalmente las peri�dicas son muy lentas, y no se
                       %utilizan.
                       %v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                       %   'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],'m_DefMacro',[]);
                       switch MOD_TYPE
                           case 1 % FIRST REDUCTION
                               v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'EPS_ROMI',[],'eta',[],'c_GdlCond',[],'Fint',[],...
                                   'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],...
                                   'c_DefMacro',[]);    %,'m_TensProy',[]
                           case 2 % SECOND REDUCTION
                               v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'EPS_ROMI',[],'eta',[],'c_GdlCond',[],'Fint',[],...
                                   'e_VarEst',[],'e_VarEst_ROMII',[],'e_VarAux',[],'e_VarAux_ROMII',[],'m_LinCond',[],'doff',[],'dofl',[],...
                                   'c_DefMacro',[]);    %,'m_TensProy',[]
                           otherwise % HF AND TRAINING CASES
                               v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                                   'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],...
                                   'c_DefMacro',[]);    %,'m_TensProy',[]
                       end
                   else
                       e_DatSetMicro = e_DatMatSet.e_DatSet;
                       e_VGMicro = e_DatMatSet.e_VG;
                       if e_VGMicro.nStage>1||e_VGMicro.NTLoad(1)>1
                           warning(['Inicializaci�n variables hist�ricas: En problemas multiescala las condiciones de borde de la microcelda ',...
                               'corresponde al Stage No 1 y Load No 1, los demas stages no se consideran durante el analisis.']);
                       end
                       
                       m_CondBord=e_VGMicro.m_CondBord.m_SCondBord{1};
                       e_VGMicro = rmfield(e_VGMicro,'m_CondBord');
                       e_VGMicro.m_CondBord = m_CondBord;
                       
                       %Como es el mismo set, se asume que las celdas unitarias de todos los puntos
                       %de gauss de todos los elementos tiene las mismas condiciones de borde iniciales.
                       %Como las condiciones de borde luego puede variar seg�n la celda (si bifurca o
                       %no), se las almacena como variables hist�ricas (en lugar en las e_DatSet, que se
                       %asume que son datos fijos).
                       %Se descarta todo las matrices vfix, m_InvCRR y doffCondCte que se refieren a los
                       %desplazamientos impuestos de grados de libertad fijos, porque se asumen que son
                       %nulos. Ver si hay alguna tipo de celda unitaria donde se necesario considerar
                       %algo distinto a esto.
                       % [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro.m_CondBord,...
                       %   e_DatMatSet.xx,e_VGMicro.ndoft,e_VGMicro.ndn,e_VGMicro.ndime,e_VGMicro.nSet,...
                       %   e_DatSetMicro,e_VGMicro.m_ConecFront);
                       
                       switch MOD_TYPE
                           case 1 % FIRST REDUCTION
                               
                               [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro.m_CondBord,...
                                   e_DatMatSet.xx,e_VGMicro.ndoft,e_VGMicro.ndn,e_VGMicro.ndime,e_VGMicro.nSet,...
                                   e_DatSetMicro,e_VGMicro.m_ConecFront,e_VGMicro.RotacionGlobal);
                               v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                                   'EPS_ROMI',zeros(e_VGMicro.nModesEPS_TOT,1),'eta',zeros(e_VGMicro.ntens-1,1),'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                                   'Fint',zeros(e_VGMicro.ndoft,1),... % incluye las tres componentes del multiplicador algoritmico
                                   'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro.nSet),...
                                   'e_VarAux',f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet,e_DatMatSet),...
                                   'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                                   'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                                   e_DatSetMicro,'UniformOutput',false)});   %'m_TensProy',[]
                               
                           case 2 % SECOND REDUCTION
                               % hrom parameters (redefining the hf parameters for allowing the hrom simulation)
                               % Definiciones propias del modelo reducido
                               nModesEPS_TOT = e_DatSet(iSet).e_DatMat.ROM_II.nModesEPS_TOT ;
                               e_DatSet_ROM = e_DatSet(iSet).e_DatMat.ROM_II.e_DatSet_ROM ;
                               
                               v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                                   'EPS_ROMI',zeros(nModesEPS_TOT,1),'eta',zeros(e_VGMicro.ntens-1,1),...
                                   'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                                   'Fint',zeros(e_VGMicro.ndoft,1),... % incluye las tres componentes del multiplicador algoritmico
                                   'e_VarEst',[],'e_VarEst_ROMII',f_eVarEstInic(c_Field,e_DatSet_ROM,1,MOD_TYPE), ...  %,,nSet,MOD_TYPE,varargin)
                                   'e_VarAux',[],'e_VarAux_ROMII',f_eVarAuxInic(e_DatSet_ROM,1,e_DatMatSet,e_DatMatSet.e_VG),...
                                   'm_LinCond',[],'doff',[],'dofl',[],...
                                   'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                                   e_DatSetMicro,'UniformOutput',false)});
                               
                           otherwise % HF AND TRAINING CASES
                               
                               % Se ha modificado la estructura m_CondBord para incluir la posibilidad de usar multiples stages, por
                               % ese motivo se debe cambiar la estructura en el dato de entrada
                               [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro.m_CondBord.m_SCondBord{1},...
                                   e_DatMatSet.xx,e_VGMicro.ndoft,e_VGMicro.ndn,e_VGMicro.ndime,e_VGMicro.nSet,...
                                   e_DatSetMicro,e_VGMicro.m_ConecFront);
                               %Las variables u (desplazamiento fluctuante micro), c_GdlCond (variable condensada
                               %micro) y e_VarAux (variables auxiliares micro) se hacen variables tipo hist�ricas
                               %ya que en la celda unitaria se aplica una deformaci�n macro variable en cada
                               %iteraci�n macro mientras que en la iteraci�n micro se consigue que siempre se
                               %parta de la misma condici�n inicial u (fluctuante micro).
                               %Otra posibilidad es que esta variables se vayan actualizando en cada iteraci�n
                               %macro, es decir pas�ndola como variable auxiliar macro (tambi�n se podr�a
                               %conseguir el mismo efecto si se pasa como argumento las variables new macro, y se
                               %utiliza los valores obtenidos de la misma). Esto har�a que una mala iteraci�n
                               %macro, que se aleje de la soluci�n, haga que en la pr�xima iteraci�n macro se
                               %parta en la iteraci�n micro de una condici�n inicial capaz muy alejada de la
                               %soluci�n (aunque en forma opuesta, si las iteraciones macro son correctas, si se
                               %utiliza condiciones iniciales actualizadas, se necesitar�a capaz menos
                               %iteraciones micro). Ahora almacen�ndola como variable hist�rica y usando los
                               %valores obtenidos del old macro, siempre se parte de una condici�n inicial micro
                               %que se sabe convergido en el paso previo, pero que puede exigir m�s iteraciones
                               %micro ya que puede estar m�s alejado de la soluci�n.
                               v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                                   'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                                   'Fint',zeros(e_VGMicro.ndoft,1),...
                                   'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro.nSet),...
                                   'e_VarAux',f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet,e_DatMatSet),...
                                   'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                                   'm_DefMacro',zeros(e_DatMatSet.e_VG.ntens,e_DatMatSet.e_VG.nElem));
                       end
                   end
                   
               case 51   %Modelo MultiEscala con discontinuidad fuerte (fisura cohesiva)
                  clear v_InicHvarSet
                  if nargin>4&&varargin{1}==0
                     %Se genera la estructura interna sin inicializar ninguna de las variables internas. Esto
                     %se utiliza en el newton para las variables new y evita el c�lculo de las matrices de
                     %condiciones de borde, que principalmente las peri�dicas son muy lentas, y no se
                     %utilizan.
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                        'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],'m_ElemLoc',[],...
                        'c_DefMacro',[],'omegaMicroL',[],'lMacro',[],'lMicro',[],...
                        'c_NormalesMicro',[],'longFis',[],'facNormMicro',[]);    %,'m_TensProy',[]
                  else
                     e_DatSetMicro = e_DatMatSet.e_DatSet;
                     e_VGMicro = e_DatMatSet.e_VG;
                     [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro.m_CondBord,...
                        e_DatMatSet.xx,e_VGMicro.ndoft,e_VGMicro.ndn,e_VGMicro.ndime,e_VGMicro.nSet,...
                        e_DatSetMicro,e_VGMicro.m_ConecFront,e_VGMicro.RotacionGlobal);
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                        'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                        'Fint',zeros(e_VGMicro.ndoft,1),...
                        'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro.nSet),...
                        'e_VarAux',f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet,e_DatMatSet),...
                        'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                        'm_ElemLoc',[],...
                        'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                           e_DatSetMicro,'UniformOutput',false)},...
                        'omegaMicroL',e_DatMatSet.omegaMicro,'lMacro',0,'lMicro',0,...
                        'c_NormalesMicro',{cell(e_VGMicro.nSet,2)},'longFis',0,'facNormMicro',1);   %'m_TensProy',[]
                  end
               case 53   %Modelo MultiEscala BCNA, SANTA FE
                  clear v_InicHvarSet
                  if nargin>4&&varargin{1}==0
                     %Se genera la estructura interna sin inicializar ninguna de las variables internas. Esto
                     %se utiliza en el newton para las variables new y evita el c�lculo de las matrices de
                     %condiciones de borde, que principalmente las peri�dicas son muy lentas, y no se
                     %utilizan.
                     switch MOD_TYPE
                         case 1 % FIRST REDUCTION
                             v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'EPS_ROMI',[],'eta',[],'c_GdlCond',[],'Fint',[],...
                                 'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],...
                                 'c_DefMacro',[]);    %,'m_TensProy',[]
                         case 2 % SECOND REDUCTION
                             %v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'EPS_ROMI',[],'eta',[],'Q_ROMI',[],'c_GdlCond',[],'Fint',[],...
                             %    'e_VarEst',[],'e_VarEst_ROMII',[],'e_VarAux',[],'e_VarAux_ROMII',[],'m_LinCond',[],'doff',[],'dofl',[],...
                             %    'c_DefMacro',[]);    %,'m_TensProy',[]
                             v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'EPS_ROMI',[],'eta',[],'c_GdlCond',[],'Fint',[],...
                                 'e_VarEst',[],'e_VarEst_ROMII',[],'e_VarAux',[],'e_VarAux_ROMII',[],'m_LinCond',[],'doff',[],'dofl',[],...
                                 'c_DefMacro',[]);    %,'m_TensProy',[]
                         otherwise % HF AND TRAINING CASES
                             v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                                 'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],...
                                 'c_DefMacro',[]);    %,'m_TensProy',[]
                     end
                  else
                     e_DatSetMicro = e_DatMatSet.e_DatSet;
                     e_VGMicro = e_DatMatSet.e_VG;
                     if e_VGMicro.nStage>1||e_VGMicro.NTLoad(1)>1
                         warning(['Inicializaci�n variables hist�ricas: En problemas multiescala las condiciones de borde de la microcelda ',...
                             'corresponde al Stage No 1 y Load No 1, los demas stages no se consideran durante el analisis.']);
                     end
                     
                     m_CondBord=e_VGMicro.m_CondBord.m_SCondBord{1};
                     e_VGMicro = rmfield(e_VGMicro,'m_CondBord');
                     e_VGMicro.m_CondBord = m_CondBord;
                     
                     switch MOD_TYPE
                         case 1 % FIRST REDUCTION
                      %       %[m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro.m_CondBord,...
                      %       %    e_DatMatSet.xx,e_VGMicro.ndoft,e_VGMicro.ndn,e_VGMicro.ndime,e_VGMicro.nSet,...
                      %       %    e_DatSetMicro,e_VGMicro.m_ConecFront);
                      %       v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                      %           'EPS_ROMI',zeros(e_VGMicro.nModesEPS_TOT,1),'eta',zeros(e_VGMicro.ntens-1,1),'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                      %           'Fint',zeros(e_VGMicro.ndoft,1),... % incluye las tres componentes del multiplicador algoritmico
                      %           'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro.nSet),...
                      %           'e_VarAux',f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet),...
                      %           'm_LinCond',[],'doff',[],'dofl',[],...
                      %           'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                      %           e_DatSetMicro,'UniformOutput',false)});   %'m_TensProy',[]

                             [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro.m_CondBord,...
                                 e_DatMatSet.xx,e_VGMicro.ndoft,e_VGMicro.ndn,e_VGMicro.ndime,e_VGMicro.nSet,...
                                 e_DatSetMicro,e_VGMicro.m_ConecFront,e_VGMicro.RotacionGlobal);
                             v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                                 'EPS_ROMI',zeros(e_VGMicro.nModesEPS_TOT,1),'eta',zeros(e_VGMicro.ntens-1,1),'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                                 'Fint',zeros(e_VGMicro.ndoft,1),... % incluye las tres componentes del multiplicador algoritmico
                                 'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro.nSet),...
                                 'e_VarAux',f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet,e_DatMatSet),...
                                 'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                                 'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                                 e_DatSetMicro,'UniformOutput',false)});   %'m_TensProy',[]
                             
                         case 2 % SECOND REDUCTION
                             % hrom parameters (redefining the hf parameters for allowing the hrom simulation)                             
                             % Definiciones propias del modelo reducido
                             nModesEPS_TOT = e_DatSet(iSet).e_DatMat.ROM_II.nModesEPS_TOT ;
                             e_DatSet_ROM = e_DatSet(iSet).e_DatMat.ROM_II.e_DatSet_ROM ;
                             
                             %doflMicro = true(nmodesU,1) ;
                             %doffMicro = false(nmodesU,1) ;
                             %m_LinCondMicro = eye(nmodesU);
                             %%%e_VarEstROMII = struct('u_ROMII',zeros(nmodesU,1)) ;
                             % Arreglo de variables auxiliares se conservara de tamanio
                             % e_VarAux = f_eVarAuxInic(e_DatSetMicro(iSet),1) ;
                             %%%e_VarAuxROMII = f_eVarAuxInic(e_DatSet_ROM,1) ;
                             % Variables historicas se redimensionaran en funcion del numero de puntos de gauss a integrar
                             %e_DatSetMicro(iSet).nElem = length(e_DatSet(iSet).e_DatMat.ROM_II.ROM_ELEM_SET) ;
                             %e_VarEst_ME = f_eVarEstInic(c_Field,e_DatSetMicro(iSet),e_VGMicro.nSet);
                             %%%e_VarEst_ME = f_eVarEstInic(c_Field,e_DatSet_ROM,1);
                             
                             % AN ADDITIONAL VARIABLE FOR CONVEXITY PURPOSES
                             %v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                             %    'EPS_ROMI',zeros(nModesEPS_TOT,1),'eta',zeros(e_VGMicro.ntens-1,1),'Q_ROMI',zeros(nModesEPS_TOT,1),...
                             %    'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                             %    'Fint',zeros(e_VGMicro.ndoft,1),... % incluye las tres componentes del multiplicador algoritmico
                             %    'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro.nSet),'e_VarEst_ROMII',f_eVarEstInic(c_Field,e_DatSet_ROM,1), ...
                             %    'e_VarAux',f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet),'e_VarAux_ROMII',f_eVarAuxInic(e_DatSet_ROM,1),...
                             %    'm_LinCond',[],'doff',[],'dofl',[],...
                             %    'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                             %    e_DatSetMicro,'UniformOutput',false)});
                             
                             %v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                             %    'EPS_ROMI',zeros(nModesEPS_TOT,1),'eta',zeros(e_VGMicro.ntens-1,1),...
                             %    'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                             %    'Fint',zeros(e_VGMicro.ndoft,1),... % incluye las tres componentes del multiplicador algoritmico
                             %    'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro.nSet),'e_VarEst_ROMII',f_eVarEstInic(c_Field,e_DatSet_ROM,1), ...
                             %    'e_VarAux',f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet),'e_VarAux_ROMII',f_eVarAuxInic(e_DatSet_ROM,1),...
                             %    'm_LinCond',[],'doff',[],'dofl',[],...
                             %    'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                             %    e_DatSetMicro,'UniformOutput',false)});     
                             
                             
                             v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                                 'EPS_ROMI',zeros(nModesEPS_TOT,1),'eta',zeros(e_VGMicro.ntens-1,1),...
                                 'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                                 'Fint',zeros(e_VGMicro.ndoft,1),... % incluye las tres componentes del multiplicador algoritmico
                                 'e_VarEst',[],'e_VarEst_ROMII',f_eVarEstInic(c_Field,e_DatSet_ROM,1), ...
                                 'e_VarAux',[],'e_VarAux_ROMII',f_eVarAuxInic(e_DatSet_ROM,1,e_DatMatSet),...
                                 'm_LinCond',[],'doff',[],'dofl',[],...
                                 'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                                 e_DatSetMicro,'UniformOutput',false)});                                   

                             
                         otherwise % HF AND TRAINING CASES
                             [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro.m_CondBord,...
                                 e_DatMatSet.xx,e_VGMicro.ndoft,e_VGMicro.ndn,e_VGMicro.ndime,e_VGMicro.nSet,...
                                 e_DatSetMicro,e_VGMicro.m_ConecFront);
                             v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                                 'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                                 'Fint',zeros(e_VGMicro.ndoft,1),...
                                 'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro.nSet), ...
                                 'e_VarAux',f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet,e_DatMatSet),...
                                 'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                                 'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                                 e_DatSetMicro,'UniformOutput',false)});   %'m_TensProy',[]
                     end

                  end
               case 54   %Modelo MultiEscala BCNA, SANTA FE
                  clear v_InicHvarSet
                  if nargin>4&&varargin{1}==0
                     %Se genera la estructura interna sin inicializar ninguna de las variables internas. Esto
                     %se utiliza en el newton para las variables new y evita el c�lculo de las matrices de
                     %condiciones de borde, que principalmente las peri�dicas son muy lentas, y no se
                     %utilizan.
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                        'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],...
                        'm_ElemLoc',[],...
                        'c_DefMacro',[],'omegaMicroL',[],'lMacro',[],'lMicro',[],...
                        'c_NormalesMicro',[]);    %,'m_TensProy',[]
                  else
                     e_DatSetMicro = e_DatMatSet.e_DatSet;
                     e_VGMicro = e_DatMatSet.e_VG;
                     [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro.m_CondBord,...
                        e_DatMatSet.xx,e_VGMicro.ndoft,e_VGMicro.ndn,e_VGMicro.ndime,e_VGMicro.nSet,...
                        e_DatSetMicro,e_VGMicro.m_ConecFront,e_VGMicro.RotacionGlobal);
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                        'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                        'Fint',zeros(e_VGMicro.ndoft,1),...
                        'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro.nSet),...
                        'e_VarAux',f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet,e_DatMatSet),...
                        'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                        'm_ElemLoc',[],...
                        'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                           e_DatSetMicro,'UniformOutput',false)}, ...   %'m_TensProy',[]
                        'omegaMicroL',e_DatMatSet.omegaMicro,'lMacro',0,'lMicro',0,...
                        'c_NormalesMicro',{cell(e_VGMicro.nSet,1)});   %'m_TensProy',[]

                  end
               otherwise
                  error('Inicializaci�n variables hist�ricas: Modelo constitutivo no definido.')
            end
            e_VarEst(iSet).(c_Field{iField}) = v_InicHvarSet;
         else
            e_VarEst(iSet).(c_Field{iField}) = m_InicZerosSet;
         end
      end
      %
      %Variables hist�rica del elemento
      e_VarEst(iSet).VarHistElem = zeros(nVarHistElem,nElem);
   end
   
end