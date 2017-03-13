function[datos,var]=extraer_numero_de_un_string(datos)

%Busca de la cadena de texto hasta donde el caracter tiene el simbolo de
%"=" ya que guid manda en ese formato el string, y asi poder obtener el
%el valor numerico.

var = datos;
for i = 1:length(datos);
    if datos(i) == '=';
    var = datos(1:min(5,length(datos)));
    datos = str2num(datos((i+1):length(datos)));
    return;    
    end;   
end;