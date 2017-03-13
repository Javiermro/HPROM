function hVarNew = f_DefinicionhVar(conshyp,sihvarpg,nPG)

   %Esta función inicializa las variables históricas según el modelo constitutivo, considerando que el caso
   %multiescala puede ser una array de estructura y para el estándar una matriz (es para llamarse dentro de la
   %función del elemento). Esta generación evita transferencia de datos a los nodos (habría que ver el tiempo
   %adicional agregado por el llamado de esta función).
   %sihvarpg = e_DatMatSet.sihvarpg;
   %nPG = e_DatElemSet.npg;
   %conshyp = e_DatMatSet.conshyp;   
   switch conshyp
      case {1,2,4,5,8,10,11,12,13,52,100,110}
         hVarNew = zeros(sihvarpg,nPG);
      case {50,55}
         hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
            'm_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[]);
      case 51
         hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
            'm_LinCond',[],'doff',[],'dofl',[],'m_ElemLoc',[],'c_DefMacro',[],'omegaMicroL',[],...
            'lMacro',[],'lMicro',[],'c_NormalesMicro',[],'longFis',[],'facNormMicro',[]);        %,'m_TensProy',[]
      case 53
         hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
            'm_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[]);        %,'m_TensProy',[]
      case 54
         hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],...
             'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],...
             'doff',[],'dofl',[],'m_ElemLoc',[],'c_DefMacro',[],...
             'omegaMicroL',[],'lMacro',[],'lMicro',[],'c_NormalesMicro',[],...
             'longFis',[],'facNormMicro',[] );        %,'m_TensProy',[]        
       otherwise
         error('Matrices Elementales: Variables Históricas: Inicialización: Modelo constitutivo no definido.')         
   end
   
end