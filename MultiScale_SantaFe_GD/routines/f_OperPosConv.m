function [e_DatSet,e_VarEst_new,e_VarAux,m_DsigmaHomog,e_VG] = ...
       f_OperPosConv(u,xx,e_VarEst_new,e_VarEst_old,e_VarAux,e_DatSet,c_CT,KT,m_LinCond,dofl,doff,e_VG)
   
   %Cálculos que se realizan después de la convergencia del newton. Generalmente son cálculos
   %correspondiente a la parte explícita del modelo y que se utiliza en el paso siguiente.
   %Esta función se coloca después de la impresión de resultados, ya que todo lo que se calcule en
   %esta función se utiliza en el paso siguiente y por lo tanto no se utilizó en el cálculo del
   %paso actual.
   %Si son operaciones que se utiliza para postproceso, como por ejemplo la tracción del elemento
   %SDA, que debe ser impresa en el mismo paso, no debe colocarse acá, ya que sino la impresión
   %estaría un paso atrasado.
   
   %CRACK PATH FIELD EVALUATION
   if e_VG.exist_CrackPath
       [e_DatSet,e_VarEst_new,e_VarAux,e_VG] = f_OperPosConv_CPF...
          (u,xx,e_VarEst_new,e_VarAux,e_DatSet,e_VG);
   end
   
   if ~e_VG.isMICRO.MICRO
       m_DsigmaHomog = zeros(e_VG.ntens,1);
       
       %BIFURCATION EVALUATION AND NORMAL SELECTION
       [e_VarEst_new,e_VarAux] = f_OperPosConv_BNS(u,xx,e_VarEst_new,e_VarEst_old,e_VarAux,e_DatSet,...
           c_CT,e_VG);
       
   else
      
%        [isBif,minQ] = MicroBifAnalysis(KT,c_CT,m_LinCond,dofl,doff,e_DatSet,e_VG); %JLM lo comento por ahora, luego hay que hacerlo
%        isBif= 1; minQ=1;
%        if isBif
%            fprintf('------ BIFURCATION - SMALL SCALE ------ \n');
%        else
%            fprintf('MICRO: Minimum of the determinant of the acustic tensor Qhom: %f\n',minQ)
%        end
       
       % STRESS HOMOGENIZATION - MICROSCALE MODEL
       
       if e_VG.MOD_TYPE==2 % HROM
           
           %m_DsigmaHomog = e_DatMatSetMacro.ROM_II.Q_oper*(e_VarEst_new.sigma(:)) ;
           nCOMP_TENS = e_VG.ROM_II.nCOMP_TENS;
           %W_STRESS   = e_DatMatSetMacro.ROM_II.W_STRESS;
           sigma_GEN = e_VarEst_new.sigma(:) ;

           m_DsigmaHomog = e_VG.ROM_II.Q_oper{1}*(sigma_GEN(nCOMP_TENS{1})) + ...
               e_VG.ROM_II.Q_oper{2}*(sigma_GEN(nCOMP_TENS{2})) ;
           
       else % ROM & HF
           %Para considerar que en el caso del cuadrangulo Q1 se
           %debe integrar las tensiones estabilizadas.
           nSet       = e_VG.nSet;
           c_Tens_old = cell(nSet,1);
           c_Tens_new = cell(nSet,1);
           for iSet = 1:nSet
               e_DatElemSet = e_DatSet(iSet).e_DatElem;
               eltype = e_DatElemSet.eltype;
               switch eltype
                   case {2,4,8,31,32,108}
                       c_Tens_old{iSet} = e_VarEst_old(iSet).sigma;
                       c_Tens_new{iSet} = e_VarEst_new(iSet).sigma;
                   case 20
                       c_Tens_old{iSet} = e_VarEst_old(iSet).VarHistElem;
                       c_Tens_new{iSet} = e_VarEst_new(iSet).VarHistElem;
                   otherwise
                       error(['Modelo Multiescala: Homogeneizaciï¿½n: Tensiones de homogeneizaciï¿½n: Modelo ',...
                           'constitutivo no definido.'])
               end
           end

           %Homogenizacion de las tensiones micro
           m_DsigmaHomog = f_HomogArea(c_Tens_new,e_VG.ntens,e_VG.omegaMicro,{e_DatSet.m_DetJT},e_DatSet,e_VG);

           % Normal and Beta homogenized from small scale
          %[vectVHElem_new] = Normal_and_beta_MICRO(xx,e_VarEst_new,...
          %    e_VarEst_old,e_VarAux,e_VG,e_DatSet,m_DsigmaHomog);       
   
       end
   end

   if e_VG.exist_CrackPath
   % Gradient of Phi and selection of lonely points in each finite element
       e_VarAux = ...
         SDA_Properties(xx,e_VG,e_VarEst_new,e_DatSet,e_VarAux);
   end


end