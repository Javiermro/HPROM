function hVarNew = f_DefinicionhVar(conshyp,sihvarpg,nPG,varargin)

   %Esta funci�n inicializa las variables hist�ricas seg�n el modelo constitutivo, considerando que el caso
   %multiescala puede ser una array de estructura y para el est�ndar una matriz (es para llamarse dentro de la
   %funci�n del elemento). Esta generaci�n evita transferencia de datos a los nodos (habr�a que ver el tiempo
   %adicional agregado por el llamado de esta funci�n).
   %sihvarpg = e_DatMatSet.sihvarpg;
   %nPG = e_DatElemSet.npg;
   %conshyp = e_DatMatSet.conshyp;
   
   % Flag for hiperreduction purposes
   if nargin == 4
       MOD_TYPE = varargin{1};
   end
   
   switch conshyp
      case {1,2,4,5,8,10,11,12,13,52}
         hVarNew = zeros(sihvarpg,nPG);
      case 50
         %hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
         %   'm_LinCond',[],'doff',[],'dofl',[],'m_DefMacro',[]);
          switch MOD_TYPE
              case 1
                  hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'EPS_ROMI',[],'eta',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],...
                      'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[]);        %,'m_TensProy',[]
              case 2
                  hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'EPS_ROMI',[],'eta',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],...
                      'e_VarEst_ROMII',[],'e_VarAux',[],'e_VarAux_ROMII',[],'m_LinCond',[],...
                      'doff',[],'dofl',[],'c_DefMacro',[]);        %,'m_TensProy',[]              
              otherwise
                  hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
                      'm_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[]);        %,'m_TensProy',[]
          end         
      case 51
         hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
            'm_LinCond',[],'doff',[],'dofl',[],'m_ElemLoc',[],'c_DefMacro',[],'omegaMicroL',[],...
            'lMacro',[],'lMicro',[],'c_NormalesMicro',[],'longFis',[],'facNormMicro',[]);        %,'m_TensProy',[]
      case {53,55}
          switch MOD_TYPE
              case 1
                  hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'EPS_ROMI',[],'eta',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],...
                      'e_VarAux',[],'m_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[]);        %,'m_TensProy',[]
              case 2
                  hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'EPS_ROMI',[],'eta',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],...
                      'e_VarEst_ROMII',[],'e_VarAux',[],'e_VarAux_ROMII',[],'m_LinCond',[],...
                      'doff',[],'dofl',[],'c_DefMacro',[]);        %,'m_TensProy',[]
                  
                  % FOR CONVEXITY PURPOSES
                  %hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'EPS_ROMI',[],'eta',[],'Q_ROMI',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],...
                  %    'e_VarEst_ROMII',[],'e_VarAux',[],'e_VarAux_ROMII',[],'m_LinCond',[],...
                  %    'doff',[],'dofl',[],'c_DefMacro',[]);        %,'m_TensProy',[]                  
              otherwise
                  hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],'e_VarEst',[],'e_VarAux',[],...
                      'm_LinCond',[],'doff',[],'dofl',[],'c_DefMacro',[]);        %,'m_TensProy',[]
          end
      case 54
         hVarNew(1:sihvarpg,1:nPG) = struct('u',[],'c_GdlCond',[],'Fint',[],...
             'e_VarEst',[],'e_VarAux',[],'m_LinCond',[],...
             'doff',[],'dofl',[],'m_ElemLoc',[],'c_DefMacro',[],...
             'omegaMicroL',[],'lMacro',[],'lMicro',[],'c_NormalesMicro',[]);        %,'m_TensProy',[]        
       otherwise
         error('Matrices Elementales: Variables Hist�ricas: Inicializaci�n: Modelo constitutivo no definido.')         
   end
   
end