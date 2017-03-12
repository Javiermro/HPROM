function PHI = AlternateOrderBasis(U_B,U,imod_stress);
% See f_MatBT_ROM_II

 
 
 dbstop('7')
    PHI = zeros(size(U_B,1),imod_stress+size(U_B,2));
    ipppl = 1;
    iUUU = 1 ;jUUU =1; SSS_DONE = 0 ; UUU_DONE = 0 ;
    while ipppl<=imod_stress+size(U_B,2) & (SSS_DONE == 1 | UUU_DONE == 0)
        
        try
            PHI(:,ipppl) = U(:,iUUU); iUUU = iUUU + 1; ipppl = ipppl + 1;
        catch
            SSS_DONE = 1;
        end;
        try
            PHI(:,ipppl) = U_B(:,jUUU); ipppl = ipppl + 1;jUUU =jUUU +1   ;
        catch
            UUU_DONE = 1;
        end;
        
        
    end
 