function psi_value = get_psi_value_Stage(m_Loads,time,e_VG)

%Dtime = e_VG.Dtime;
%max_time = e_VG.max_time;

% Numero de cargas que se tienen en el mismo Stage
nLoads    = size(m_Loads,1);
psi_value = zeros(nLoads,1);

% Loop sobre las cargas del Stage
for iLoad=1:nLoads
    funbc = m_Loads{iLoad};
    psi_value(iLoad) = get_psi_value(funbc,time,e_VG);
end

end

function   psi_value = get_psi_value(funbc,time,e_VG)

%Se esta considerando que el time ingresado siempre es mayor que cero, por lo que no se hace
%ninguna verificaci�n si no se cumple con esto. Tambi�n se asume que est� ordenados los tiempos
%ingresados, pero si no siempre se usa el primer n�mero mayor (seg�n est�n ordenados).
npfun = size(funbc,1);

% Cuando time es menor a la primera abscisa (tiempo) de la funci�n.
if time<funbc(1,1)
    %Puede devolver valores de Psi menores que cero (ver que hacer con esto).
    warning(['An�lisis no lineal: Determinaci�n del factor Psi: El tiempo de c�lculo es menor',...
        ' al primer punto de la funci�n Psi definida: Se utiliza una extrapolaci�n de',...
        ' la primer recta ingresada.']) %#ok<WNTAG>
end
% Cuando time est� dentro de los valores de la funci�n definida.
for i = 2:npfun
    %Esta forma realizar permite que se pueda ingresar funciones con saltos en Psi, solo se debe
    %ingresar para el mismo tiempo dos valores y que los mismos ordenados de menor a mayor
    %(siempre se busca el primer tiempo mayor a time, y se interpola entre este y el anterior, y
    %luego se sale de la funci�n).
    if time<=funbc(i,1)
        %Interpolaci�n lineal entre el valor anterior (i-1) y el siguiente.
        psi_value = funbc(i-1,2)+(time-(funbc(i-1,1)))*(funbc(i,2)-funbc(i-1,2))/...
            (funbc(i,1)-funbc(i-1,1));
        %Cuando se encuentra el primer tiempo valor mayor o igual a time, y se sale (esto permite
        %imponer funciones con salto y que no se sobreescriba el valor de Psi).
        return
    end
end
% Cuando time es mayor a los valores de la funci�n definida.
%Por el return anterior ac� se pasa s�lo si ocurre que el time es mayor.
%if time>funbc(npfun,1)
psi_value = funbc(npfun-1,2)+(time-(funbc(npfun-1,1)))*(funbc(npfun,2)-funbc(npfun-1,2))/...
    (funbc(npfun,1)-funbc(npfun-1,1));
warning(['An�lisis no lineal: Determinaci�n del factor Psi: Se super� el tiempo hasta ',...
    'donde se defini� la funci�n: Se utiliza una extrapolaci�n de la �ltima recta ingresada.']) %#ok<WNTAG>
%end

end
