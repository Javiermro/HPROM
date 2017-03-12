function PHI = ConsExpMatrix(ALTERNATE_BASIS_STRESS_B,imod_stress,U_B,U)
% See f_MatBT_ROM_III

%% JAHO,9-May-2012
%The basis matrix PHI is formed, on the one hand,  by the  basis vectors
%emanating from % the spectral analysis of the snapshot stress matrix, and on the other hand, by the basis vectors
% obtained  from the SVD of the B-strain matrix.
% Bear in mind that the points selected by greedy algorithms depends upon the order in which the
% basis vectors appear in the   basis matrix ---they are   hierarchical methods.

% The natural option is to group these vectors by simply putting together
% both matrices:
% --> PHI =[U(:,1:imod_stress)   U_B ];
% However, this strategy conflicts with the requirement of "consistenchy"
% that one would expect from a model reduction methods (the more bases one use, the better is the approximation), simply because the
% points/indices associated to the basis matrix placed in the second
% position is likely to change as the level(s) of truncation increases.
% To avoid this, use option ALTERNATE_BASIS_STRESS_B=1. By enabling this
% option, vector bases of both type appears intercalated in the basis
% matrix.
%dbstop('341')

% JAHO, 5-sept-2012. The above remarks are, in principle, irrelevant for
%                    BASIS_EXPANSION = 'SIMPLIFIED_METHOD'


if  ALTERNATE_BASIS_STRESS_B == 0 || (ALTERNATE_BASIS_STRESS_B > 1 &&  imod_stress == size(U_B,2) && imod_stress<ALTERNATE_BASIS_STRESS_B )
    PHI =[U(:,1:imod_stress)   U_B ];
    
elseif  ALTERNATE_BASIS_STRESS_B == 1 &&  imod_stress == size(U_B,2)
    % PHI = AlternateOrderBasis(U_B,U,imod_stress);
    % This option is only valid if imod_stress = size(U_B,2)
    iSSS = 1:2:2*imod_stress-1;
    iUUU = 2:2:2*size(U_B,2);
    PHI = zeros(size(U_B,1),imod_stress+size(U_B,2));
    PHI(:,iSSS) = U(:,1:imod_stress) ;
    PHI(:,iUUU) = U_B ;
elseif  ALTERNATE_BASIS_STRESS_B > 1 &&  imod_stress == size(U_B,2)
    % PHI = AlternateOrderBasis(U_B,U,imod_stress);
    % This option is only valid if imod_stress = size(U_B,2)
    % The first e_VG.ALTERNATE_BASIS_STRESS_B*2 basis are ordered as usual
    
    nbasis_mant = ALTERNATE_BASIS_STRESS_B  ;
    U_B_ini = U_B(:,1:nbasis_mant) ;
    U_ini = U(:,1:nbasis_mant);
    
    iSSS = 2*nbasis_mant+1:2:2*(imod_stress)-1;
    iUUU = 2*nbasis_mant+2:2:2*(size(U_B,2));
    PHI = zeros(size(U_B,1),imod_stress+size(U_B,2));
    PHI(:,1:2*nbasis_mant) = [U_ini U_B_ini] ;
    PHI(:,iSSS) = U(:,nbasis_mant+1:imod_stress) ;
    PHI(:,iUUU) = U_B(:,nbasis_mant+1:end) ;
else
    error('Option not implemented or not valid for this case')
end
