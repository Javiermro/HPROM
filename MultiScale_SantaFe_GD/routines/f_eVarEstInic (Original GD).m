function e_VarEst = f_eVarEstInic(c_Field,e_DatSet,e_VG,varargin)
   
%Inicializa la estructuras de variable de estado
   nSet = e_VG.nSet;
   nField = length(c_Field);
   c_StructDef = cat(1,c_Field,repmat({[]},1,nField));
   %e_VarEst(1:nSet,1) = struct(c_StructDef{:});
   %Si son estructuras vacías también se puede usar
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
                  %Se inicializa la variable histórica para el modelo constitutivo 11
                  if conshyp==11
                     %Se inicializa la variables r_old(5), q_old(6), rInic_old(8) y qInic_old(9). Todas tienen
                     %que ser igual a r0.
                     sihvarpg = e_DatMatSet.sihvarpg;
                     if e_DatMatSet.tit==1
                        %En el caso de ser exponencial también se inicializa las rInic_old(end-1) y
                        %qInic_old(end).                        
                        m_Ind = [5;6;sihvarpg-1;sihvarpg];
                     else
                        m_Ind = [5;6];
                     end
                     v_InicHvarSet(bsxfun(@plus,m_Ind,0:sihvarpg:sihvare-1),:) = e_DatMatSet.r_0;                       
                  elseif conshyp==110
                     sihvarpg = e_DatMatSet.sihvarpg;
                     m_Ind = [3;4;5];
                     v_InicHvarSet(bsxfun(@plus,m_Ind,0:sihvarpg:sihvare-1),:) = 1;                       
                  end                      
               case {50,55}   %Modelo MultiEscala continuo
                  clear v_InicHvarSet
                  if nargin>3&&varargin{1}==0
                     %Se genera la estructura interna sin inicializar ninguna de las variables internas. Esto
                     %se utiliza en el newton para las variables new y evita el cálculo de las matrices de
                     %condiciones de borde, que principalmente las periódicas son muy lentas, y no se
                     %utilizan.
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                        'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[]);
                  else
                     e_DatSetMicro = e_DatMatSet.e_DatSet;
                     e_VGMicro = e_DatMatSet.e_VG;
                     %Como es el mismo set, se asume que las celdas unitarias de todos los puntos
                     %de gauss de todos los elementos tiene las mismas condiciones de borde iniciales.
                     %Como las condiciones de borde luego puede variar según la celda (si bifurca o
                     %no), se las almacena como variables históricas (en lugar en las e_DatSet, que se
                     %asume que son datos fijos).
                     %Se descarta todo las matrices vfix, m_InvCRR y doffCondCte que se refieren a los
                     %desplazamientos impuestos de grados de libertad fijos, porque se asumen que son
                     %nulos. Ver si hay alguna tipo de celda unitaria donde se necesario considerar
                     %algo distinto a esto.
                     [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro,...
                        e_DatMatSet.xx,e_DatSetMicro,e_VGMicro.m_ConecFront);
                     %Las variables u (desplazamiento fluctuante micro), c_GdlCond (variable condensada
                     %micro) y e_VarAux (variables auxiliares micro) se hacen variables tipo históricas
                     %ya que en la celda unitaria se aplica una deformación macro variable en cada
                     %iteración macro mientras que en la iteración micro se consigue que siempre se
                     %parta de la misma condición inicial u (fluctuante micro).
                     %Otra posibilidad es que esta variables se vayan actualizando en cada iteración
                     %macro, es decir pasándola como variable auxiliar macro (también se podría
                     %conseguir el mismo efecto si se pasa como argumento las variables new macro, y se
                     %utiliza los valores obtenidos de la misma). Esto haría que una mala iteración
                     %macro, que se aleje de la solución, haga que en la próxima iteración macro se
                     %parta en la iteración micro de una condición inicial capaz muy alejada de la
                     %solución (aunque en forma opuesta, si las iteraciones macro son correctas, si se
                     %utiliza condiciones iniciales actualizadas, se necesitaría capaz menos
                     %iteraciones micro). Ahora almacenándola como variable histórica y usando los
                     %valores obtenidos del old macro, siempre se parte de una condición inicial micro
                     %que se sabe convergido en el paso previo, pero que puede exigir más iteraciones
                     %micro ya que puede estar más alejado de la solución.
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                        'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                        'Fint',zeros(e_VGMicro.ndoft,1),...
                        'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro),...
                        'e_VarAux',f_eVarAuxInic(e_DatMatSet.xx,e_DatSetMicro,e_VGMicro),...
                        'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                        'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.nElem),e_DatSetMicro,...
                           'UniformOutput',false)});
                        %Notar que la deformación para grandes deformaciones (Gradiente de deformación F) no
                        %se inicializa en la identidad ya que este es pisado siempre en la función f_RMap_ME.
                        %{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),e_DatSetMicro,'UniformOutput',false)});
                  end
               case 51   %Modelo MultiEscala con discontinuidad fuerte (fisura cohesiva)
                  clear v_InicHvarSet
                  if nargin>3&&varargin{1}==0
                     %Se genera la estructura interna sin inicializar ninguna de las variables internas. Esto
                     %se utiliza en el newton para las variables new y evita el cálculo de las matrices de
                     %condiciones de borde, que principalmente las periódicas son muy lentas, y no se
                     %utilizan.
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                        'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],'m_ElemLoc',[],...
                        'c_DefMacro',[],'omegaMicroL',[],'lMacro',[],'lMicro',[],...
                        'c_NormalesMicro',[],'longFis',[],'facNormMicro',[]);    %,'m_TensProy',[]
                  else
                     e_DatSetMicro = e_DatMatSet.e_DatSet;
                     e_VGMicro = e_DatMatSet.e_VG;
                     [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro,...
                        e_DatMatSet.xx,...
                        e_DatSetMicro,e_VGMicro.m_ConecFront);
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                        'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                        'Fint',zeros(e_VGMicro.ndoft,1),...
                        'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro),...
                        'e_VarAux',f_eVarAuxInic(e_DatMatSet.xx,e_DatSetMicro,e_VGMicro),...
                        'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                        'm_ElemLoc',[],...
                        'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                           e_DatSetMicro,'UniformOutput',false)},...
                        'omegaMicroL',e_DatMatSet.omegaMicro,'lMacro',0,'lMicro',0,...
                        'c_NormalesMicro',{cell(e_VGMicro.nSet,2)},'longFis',0,'facNormMicro',1);   %'m_TensProy',[]
                  end
               case 53   %Modelo MultiEscala BCNA, SANTA FE
                  clear v_InicHvarSet
                  if nargin>3&&varargin{1}==0
                     %Se genera la estructura interna sin inicializar ninguna de las variables internas. Esto
                     %se utiliza en el newton para las variables new y evita el cálculo de las matrices de
                     %condiciones de borde, que principalmente las periódicas son muy lentas, y no se
                     %utilizan.
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                        'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],...
                        'c_DefMacro',[]);    %,'m_TensProy',[]
                  else
                     e_DatSetMicro = e_DatMatSet.e_DatSet;
                     e_VGMicro = e_DatMatSet.e_VG;
                     [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro,...
                        e_DatMatSet.xx,...
                        e_DatSetMicro,e_VGMicro.m_ConecFront);
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                        'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                        'Fint',zeros(e_VGMicro.ndoft,1),...
                        'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro),...
                        'e_VarAux',f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet),...
                        'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                        'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                           e_DatSetMicro,'UniformOutput',false)});   %'m_TensProy',[]
                  end
               case 54   %Modelo MultiEscala BCNA, SANTA FE
                  clear v_InicHvarSet
                  if nargin>3&&varargin{1}==0
                     %Se genera la estructura interna sin inicializar ninguna de las variables internas. Esto
                     %se utiliza en el newton para las variables new y evita el cálculo de las matrices de
                     %condiciones de borde, que principalmente las periódicas son muy lentas, y no se
                     %utilizan.
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',[],'c_GdlCond',[],'Fint',[],...
                        'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],...
                        'm_ElemLoc',[],...
                        'c_DefMacro',[],'omegaMicroL',[],'lMacro',[],'lMicro',[],...
                        'c_NormalesMicro',[] ,...
                        'longFis',[],'facNormMicro',[] );   %'m_TensProy',[]
                  else
                     e_DatSetMicro = e_DatMatSet.e_DatSet;
                     e_VGMicro = e_DatMatSet.e_VG;
                     [m_LinCondMicro,~,~,doffMicro,doflMicro] = f_CondBord(e_VGMicro,...
                        e_DatMatSet.xx,...
                        e_DatSetMicro,e_VGMicro.m_ConecFront);
                     v_InicHvarSet(1:sihvare,1:nElem) = struct('u',zeros(e_VGMicro.ndoft,1),...
                        'c_GdlCond',{f_cVarCondInic(e_DatSetMicro,e_VGMicro)},...
                        'Fint',zeros(e_VGMicro.ndoft,1),...
                        'e_VarEst',f_eVarEstInic(c_Field,e_DatSetMicro,e_VGMicro),...
                        'e_VarAux',f_eVarAuxInic(e_DatSetMicro,e_VGMicro.nSet),...
                        'm_LinCond',m_LinCondMicro,'doff',doffMicro,'dofl',doflMicro,...
                        'm_ElemLoc',[],...
                        'c_DefMacro',{arrayfun(@(x)zeros(e_VGMicro.ntens,x.e_DatElem.npg,x.nElem),...
                           e_DatSetMicro,'UniformOutput',false)}, ...   %'m_TensProy',[]
                        'omegaMicroL',e_DatMatSet.omegaMicro,'lMacro',0,'lMicro',0,...
                        'c_NormalesMicro',{cell(e_VGMicro.nSet,1)},...
                        'longFis',0,'facNormMicro',1 );
                  end
               otherwise
                  error('Inicialización variables históricas: Modelo constitutivo no definido.')
            end
            e_VarEst(iSet).(c_Field{iField}) = v_InicHvarSet;
         else
            e_VarEst(iSet).(c_Field{iField}) = m_InicZerosSet;
         end
      end
      %
      %Variables histórica del elemento
      e_VarEst(iSet).VarHistElem = zeros(nVarHistElem,nElem);
   end
   
end