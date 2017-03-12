function [SSd,SS] =RetrieveSingV(NAME_SAVE_BASIS_DISP,NAME_SAVE_BasisR,INCLUDE_SVDISPL)

SSd = [] ;
if  INCLUDE_SVDISPL ==1
    try
        load(NAME_SAVE_BASIS_DISP,'Sp') ;
        SSd = Sp ;
    catch
        try
            load(NAME_SAVE_BASIS_DISP,'Sd') ;
            SSd = Sd ;
        catch
            % When  DATAGEN.CONSISTENT_INELASTIC_ELASTIC_DECOMP ==1
            load(NAME_SAVE_BASIS_DISP,'Se','Si') ;
            SSd = [Se ; Si] ;
        end
    end
    
end
try
    load(NAME_SAVE_BasisR,'Sp') ;
    SS = Sp ;
catch
    try
        load(NAME_SAVE_BasisR,'Sd') ;
        SS = Sd ;
    catch
        load(NAME_SAVE_BasisR,'Se','Si') ;
        SS = [Se ; Si] ;
    end
    
end