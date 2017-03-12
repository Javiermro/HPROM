function [REPEATED_POINTS SELEC_POINTS B_SELECTED cond_NUMB_WITHOUT INDICES_EXPAND cond_NUMB_WITH INDICES_ORIGINAL,...
    INDICES_STRESS_RECONST,INDICES_GAUSS_RECONST,SELECT_ELEMENT_STRESS]...
    = Greedy_Hierarchical_cycle(PHI_EXPANDED,ng,npg,ntens, PARAMETER_GREEDY,NMODES_stress,NMODE_disp,...
    U_left,S_singular,BASIS_EXPANSION)


%dbstop('8')
if nargin == 0
    load('tmpp.mat')
    
    % PARAMETER_GREEDY.CRITERION_NUMBER_ADDPOINTS = 'NCOND';
    % This parameters indicate in routine Greedy_Hierarchical_cycle.m  how many
    %  points have to be introduced
    % --> ALL:  nadditional = Nbasis_U (default option)
    % ---> FIXED: Prescribed by the user  in (PARAMETER_GREEDY.NUMBER_ADDITIONAL)
    % ---> NCOND: --> Number of points required to render the condition number
    % below a given threshold (PARAMETER_GREEDY.NCOND_THRESHOLD)
    %PARAMETER_GREEDY.NUMBER_ADDITIONAL = 1;  % Only used if PARAMETER_GREEDY.CRITERION_NUMBER_ADDPOINTS = 'FIXED';
    %PARAMETER_GREEDY.NCOND_THRESHOLD = 30;  % Only used if PARAMETER_GREEDY.CRITERION_NUMBER_ADDPOINTS = 'FIXED';
end


% STEP 1: Selecting the first NMODES_stress points (according to optimal reconstruction criterion)
PHI = U_left(:,1:NMODES_stress) ;
K = NMODES_stress ;
PSI = U_left(:,NMODES_stress+1:end);
S = diag(S_singular(NMODES_stress+1:end));

PARAMETER_GREEDY.ng = ng ;
V = [] ;
%PARAMETERS.JAHO_OBJ_MULTnew.NAME_WS_IDENT = 'CoarseMeshNnew' ;  % --->

%dbstop('33')
[INDICES_FINAL MIN_TOTAL INDICES_GLO  min_comb INDICES_FINAL_GAUSS]= GreedyM_JAHOOBJ_rec_MULTI_NEW(PHI,K,PSI,S,...
    PARAMETER_GREEDY.LTRAIL,V,PARAMETER_GREEDY.NCOMB,PARAMETER_GREEDY);

% RELEVANT OUTPUTS:

% -->  INDICES_FINAL: Indexes of selected stress components (associated to gauss points INDICES_FINAL_GAUSS);
% -->   INDICES_FINAL_GAUSS: Indexes of selected gauss points  (the first NMODES_stress)



%
%dbstop('46')
INDICES_STRESS_RECONST = INDICES_FINAL ;      % For method --> BASIS_EXPANSION = 'SIMPLIFIED_METHOD'
INDICES_GAUSS_RECONST = INDICES_FINAL_GAUSS ;  % % JAHO, this strategy has proved unefficient; delete it when "cleaning" the code
%SELEC_POINTS(i) => [number of finite element   number of gauss point within the finite element]
IGGG =  mod(INDICES_FINAL_GAUSS,npg) ;
[ind1] = find(IGGG==0)  ;
if ~isempty(ind1)
    IGGG(ind1) = npg ;
end

elements = (INDICES_FINAL_GAUSS-IGGG)/npg+1;

SELECT_ELEMENT_STRESS = [elements', IGGG'];



%  The first was intended to select the integration points
%  that minimize the error between the HF stress snapshots  and the "reconstructed"
%  stress snapshots. In the second step, however, we have to discriminate
%  whether we are using BASIS_EXPANSION = 'NO_EXPANSION' or
%  BASIS_EXPANSION = 'SIMPLIFIED_METHOD'.



% In the case BASIS_EXPANSION = 'NO_EXPANSION', we work with the "expanded" basis matrix. Hence, selection of the
% "B-strain" points should also account for the stress bases and the
% corresponding sampling points. By contrast, in the case BASIS_EXPANSION =
% 'SIMPLIFIED_METHOD', reconstruction of stress tensor is entirely
% "decoupled" from the solution of the nonlinear, reduced-order equilibrium
% equations

switch  BASIS_EXPANSION
    case 'SIMPLIFIED_METHOD'
        % 5-Sept-2012 -------- **********
        % ........... -------- **********
        
        [REPEATED_POINTS B_SELECTED SELEC_POINTS INDICES_EXPAND INDICES_ORIGINAL cond_NUMB_WITH cond_NUMB_WITHOUT] =BMATRIXPointsSelection(INDICES_FINAL,INDICES_FINAL_GAUSS,PHI,NMODE_disp,PARAMETER_GREEDY,PHI_EXPANDED,ng,...
            NMODES_stress,npg);
        
        % IMPORTANT: Note that, in this case, B_SELECTED and SELEC_POINTS only include
        % the "B-matrix" points !!!!!
        
        % JAHO, this strategy has proved unefficient; delete it when
        % "cleaning" the code
    otherwise
        
        %
        
        
        
        
        [REPEATED_POINTS B_SELECTED SELEC_POINTS INDICES_EXPAND INDICES_ORIGINAL cond_NUMB_WITH cond_NUMB_WITHOUT] =AdditionalPointsSelection(INDICES_FINAL,INDICES_FINAL_GAUSS,PHI,NMODE_disp,PARAMETER_GREEDY,PHI_EXPANDED,ng,...
            NMODES_stress,npg);
        
        
        
        % STEP 2: Choosing the last NMODE_disp. Sampling criterion: the selection
        % must render the expanded  "Mgappy"5
        % matrix --- constructed from the block matrix formed by the selected rows of the EXPANDED, HYBRID MATRIX PHI_EXPANDED
        % invertible
        
        % How to guarantee this invertibility requirement ? We use a Greedy selection
        % procedure ---akin to Willcox's:
        
        %81    41
        %321   322   323   324   161   162   163   164
        
        
        
        
        
        %
        % PARAMETER_GREEDY.CRITERION_ADDPOINTS = 'COND_GREEDY'; % Criterion that specifies the procedure for
        % % selecting "additional sample points"
        % % (those corresponding to the B-matrix basis vectors)
        % % --> COND_GREEDY: Greedy
        % % algorithm that attemts to
        % % minimize the condition number
        % % of M = PHI_EXTED_gappy'*PHI_EXTED_gappy'
        % PARAMETER_GREEDY.CRITERION_NUMBER_ADDPOINTS = 'ALL';
        % % This parameters indicate in routine Greedy_Hierarchical_cycle.m  how many
        % %  points have to be introduced
        % % --> ALL:  nadditional = Nbasis_U (default option)
        % % ---> FIXED: Prescribed by the user  in (PARAMETER_GREEDY.NUMBER_ADDITIONAL)
        % % ---> NCOND: --> Number of points required to render the condition number
        % % below a given threshold (PARAMETER_GREEDY.NCOND_THRESHOLD)
        % PARAMETER_GREEDY.NUMBER_ADDITIONAL = 1;  % Only used if PARAMETER_GREEDY.CRITERION_NUMBER_ADDPOINTS = 'FIXED';
        % PARAMETER_GREEDY.NCOND_THRESHOLD = 100;  % Only used if
        
end

