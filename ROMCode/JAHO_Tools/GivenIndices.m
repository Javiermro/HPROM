function [REPEATED_POINTS SELEC_POINTS B_SELECTED cond_NUMB_WITHOUT INDICES_EXPAND cond_NUMB_WITH] = GivenIndices(PHI,ng,CLASS_INTEGER,PARAM_STUDY,npg,ntens,INDICES)
% See greedy_algorithms.m

PHI_RED = PHI(INDICES,:) ;
M = PHI_RED'*PHI_RED;
cond_NUMB_WITHOUT = cond(M);

%% But it is instructive (and even necessary) to inquire about the
%% condition number of the reduced PHI matrix consnstructed by adding the files
%% corresponding to the associated stress components
INDICES_EXPAND = additional_rows(ng,INDICES,CLASS_INTEGER); % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg

PHI_RED = PHI(INDICES_EXPAND,:) ;
M = PHI_RED'*PHI_RED;
switch PARAM_STUDY
    case 'LEVEL_COND';
        cond_NUMB_WITH = cond(M);
    case 'NORMIDENT'
        MM = eye(size(M))-M;
        cond_NUMB_WITH = norm(MM,'fro');
    case 'NORMIDENT_JAHO'
        
        Pselect = SelectOperator(INDICES_EXPAND,size(PHI,1));
        %
        dPHI = PHI'-M\PHI_RED'*Pselect;
        condnum =  norm(dPHI,'fro'); %sqrt(sum(sum(MIDENT.*MIDENT)));
        
        
        cond_NUMB_WITH = norm(MM,'fro');
        
end
%dbstop('215')
INDICES_EXPAND = sort(INDICES_EXPAND);
%%%%%%%%%%%
% How to determine the "points indicies" from INDICES_EXPAND
%         INDICES = points_indices(INDICES_EXPAND,ng*npg);
%         INDICES = unique(INDICES)';
%         INDICES_B = repmat(npg*INDICES,npg,1)-repmat((npg-1:-1:0)',1,size(INDICES,2));

[SELEC_POINTS,B_SELECTED,REPEATED_POINTS] = selected_points(ng*npg,npg,ntens,INDICES,CLASS_INTEGER);
%%%%%%%%%%%
