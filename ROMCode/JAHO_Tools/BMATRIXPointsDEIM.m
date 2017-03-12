function [NCOND_PHI_EXPANDED_FINAL,GAUSS_POINTS_K,INDICES_K,Mexpand] = ...
    BMATRIXPointsDEIM(INDICES_K,PHI_EXPANDED,NMODE_disp,ng,NMODES_stress,GAUSS_POINTS_K)
%---------------------------------------------------
% See  --> Greedy_Hierarchical_cycle.m  ------------
%---------------------------------------------------
%dbstop('7')
if nargin == 0
    load('/home/joaquin/USO_COMUN_MATLAB/multiescala_ALFREDO/BORTHO_DIRECT_modif/MultiEscala_Marzo_2012/jaho_routines/tmp.mat')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEIM HEURISTIC FOR DISCRIMINATING WHICH IS THE NEXT BEST SPATIAL POINTS
% IN THE SET OF "gappy"POINTS
% --- JAHO, 5, SEPT,2012
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% FIRST BASIS
% imodeB = 1 ;
% Uk = PHI_EXPANDED(:,NMODES_stress+1:NMODES_stress+imodeB) ;
%
%
% [rho ind] = max(abs(PHI(:,k)));
% INDICES(k) = ind ;

%dbstop('156')
%  hhhh = waitbar(0,'Computing norm(max(abs(r))) ...');

% We assume that Number of additional points = NMODE_disp

% In essence, this method attempts to write the "imodeB-th" mode as
% a linear combination of the preceeding NMODES_stress+(imodeB-1)-th modes

% B basis
% -------

PHI_EXPANDED = PHI_EXPANDED(:,NMODES_stress+1:end);
Uk = PHI_EXPANDED(:,1) ;
ngauss  = size(PHI_EXPANDED,1)/ng ;
fun_max = 0 ;
r = Uk ;
for igauss = 1:ngauss % This operation might be elegantly put in vectorial form ... Doesn't it ?
    ind_gauss = [(igauss-1)*ng+1: igauss*ng ];
    obj_fun = norm(r(ind_gauss));
    if obj_fun>fun_max
        fun_max = obj_fun ;
        ind = ind_gauss ;
        point_loc = igauss ;
    end
end
    INDICES_K= [ ind] ;
 GAUSS_POINTS_K = [ point_loc] ;


for imodeB = 2:NMODE_disp
    
    % The block of the Uk mode formed by the  INDICES_K rows
    U_red = Uk(INDICES_K,:);
    % Coefficients of the expansion
    M = U_red'*U_red ;
    c = M\U_red'*PHI_EXPANDED(INDICES_K,imodeB);
    % Residual
    r = PHI_EXPANDED(:,imodeB)-Uk*c ;
    
    % How to calculate the maximum ? The optimum value is now
    % a vector (and not a scalar, as in  the standard DEIM)
    % Now we check at which ``gauss point" the residual reaches its maximum
    % How to define the ``maximum" ? In which norm ?  --< euclidean norm
    % seems a reasonable choice
    fun_max = 0 ;
    for igauss = 1:ngauss % This operation might be elegantly put in vectorial form ... Doesn't it ?
        ind_gauss = [(igauss-1)*ng+1: igauss*ng ];
        obj_fun = norm(r(ind_gauss));
        if obj_fun>fun_max
            fun_max = obj_fun ;
            ind = ind_gauss ;
            point_loc = igauss ;
        end
    end
    
    
    
    Uk = PHI_EXPANDED(:,1:imodeB) ;
    INDICES_K= [INDICES_K ind] ;
    GAUSS_POINTS_K = [GAUSS_POINTS_K point_loc] ;
    
    
    
end

U_red = Uk(INDICES_K,:);
% Coefficients of the expansion
Mexpand= U_red'*U_red ;
NCOND_PHI_EXPANDED_FINAL = cond(Mexpand) ;

%[NCOND_PHI_EXPANDED_FINAL,GAUSS_POINTS_K,INDICES_K,Mexpand] = ...

