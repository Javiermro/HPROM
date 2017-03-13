function psi_value = get_psi_value(funbc,time,e_VG)

   %Dtime = e_VG.Dtime;
   %max_time = e_VG.max_time;

   %Se está considerando que el time ingresado siempre es mayor que cero, por lo que no se hace
   %ninguna verificación si no se cumple con esto. También se asume que está ordenados los tiempos
   %ingresados, pero si no siempre se usa el primer número mayor (según estén ordenados).
   npfun = size(funbc,1);

   % Cuando time es menor a la primera abscisa (tiempo) de la función.
   if time<funbc(1,1)
      %Puede devolver valores de Psi menores que cero (ver que hacer con esto).
      warning(['Análisis no lineal: Determinación del factor Psi: El tiempo de cálculo es menor',...
         ' al primer punto de la función Psi definida: Se utiliza una extrapolación de',...
         ' la primer recta ingresada.']) %#ok<WNTAG>
   end
   % Cuando time está dentro de los valores de la función definida.
   for i = 2:npfun
      %Esta forma realizar permite que se pueda ingresar funciones con saltos en Psi, solo se debe
      %ingresar para el mismo tiempo dos valores y que los mismos ordenados de menor a mayor
      %(siempre se busca el primer tiempo mayor a time, y se interpola entre este y el anterior, y
      %luego se sale de la función).
      if time<=funbc(i,1)
         %Interpolación lineal entre el valor anterior (i-1) y el siguiente.
         psi_value = funbc(i-1,2)+(time-(funbc(i-1,1)))*(funbc(i,2)-funbc(i-1,2))/...
            (funbc(i,1)-funbc(i-1,1));
         %Cuando se encuentra el primer tiempo valor mayor o igual a time, y se sale (esto permite 
         %imponer funciones con salto y que no se sobreescriba el valor de Psi).
         return
      end
   end
   % Cuando time es mayor a los valores de la función definida.
   %Por el return anterior acá se pasa sólo si ocurre que el time es mayor.
   %if time>funbc(npfun,1)
   psi_value = funbc(npfun-1,2)+(time-(funbc(npfun-1,1)))*(funbc(npfun,2)-funbc(npfun-1,2))/...
      (funbc(npfun,1)-funbc(npfun-1,1));
   warning(['Análisis no lineal: Determinación del factor Psi: Se superó el tiempo hasta ',...
      'donde se definió la función: Se utiliza una extrapolación de la última recta ingresada.']) %#ok<WNTAG>
   %end
   
end
