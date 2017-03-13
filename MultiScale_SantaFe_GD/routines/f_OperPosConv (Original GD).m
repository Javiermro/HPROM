function [e_DatSet,e_VarEst_new,e_VarAux,e_VG] = ...
       f_OperPosConv(u,xx,e_VarEst_new,e_VarEst_old,e_VarAux,e_DatSet,c_CT,e_VG)
   
   %C�lculos que se realizan despu�s de la convergencia del newton. Generalmente son c�lculos
   %correspondiente a la parte expl�cita del modelo y que se utiliza en el paso siguiente.
   %Esta funci�n se coloca despu�s de la impresi�n de resultados, ya que todo lo que se calcule en
   %esta funci�n se utiliza en el paso siguiente y por lo tanto no se utiliz� en el c�lculo del
   %paso actual.
   %Si son operaciones que se utiliza para postproceso, como por ejemplo la tracci�n del elemento
   %SDA, que debe ser impresa en el mismo paso, no debe colocarse ac�, ya que sino la impresi�n
   %estar�a un paso atrasado.
   
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