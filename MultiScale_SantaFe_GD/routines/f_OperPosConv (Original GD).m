function [e_DatSet,e_VarEst_new,e_VarAux,e_VG] = ...
       f_OperPosConv(u,xx,e_VarEst_new,e_VarEst_old,e_VarAux,e_DatSet,c_CT,e_VG)
   
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
   
    %BIFURCATION EVALUATION AND NORMAL SELECTION 
    [e_VarEst_new,e_VarAux] = f_OperPosConv_BNS(u,xx,e_VarEst_new,e_VarEst_old,e_VarAux,e_DatSet,...
      c_CT,e_VG);


   if e_VG.exist_CrackPath
   % Gradient of Phi and selection of lonely points in each finite element
       e_VarAux = ...
         SDA_Properties(xx,e_VG,e_VarEst_new,e_DatSet,e_VarAux);
   end


end