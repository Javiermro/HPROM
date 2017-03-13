function [datos,count,flag] = extraer_datos(fid,formato,tamano,flag)

    [datos,count] = fscanf(fid,formato,tamano);
    flag_guardar = flag;
    flag = 0;
    if count~=tamano & ~isempty(datos) 
        flag = flag_guardar;
        fclose(fid); 
        %return; 
    end

end


    