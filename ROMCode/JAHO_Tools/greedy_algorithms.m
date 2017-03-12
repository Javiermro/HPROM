function [SELEC_POINTS,B_SELECTED, cond_NUMB_WITHOUT, INDICES_EXPAND,  cond_NUMB_WITH, WORST_COMB,INDICES_ORIGINAL,...
    INDICES_STRESS_RECONST,INDICES_GAUSS_RECONST,SELECT_ELEMENT_STRESS] = ...
    greedy_algorithms(PHI,GREEDY_ALGOR,ntens,npg,varargin)

% Greedy algoritm (Chaturantabut_2009_Model_reduc.pdf), pag 6
if nargin == 0
    load('DATA_JAHO/PHI.mat','U');
    NBASIS = 3 ;
    PHI = U(:,1:NBASIS); clear U ;
    ng = 4;
    
    GREEDY_ALGOR = 'WILLCOX' ;
    GREEDY_ALGOR = 'WILLCOX_INIGUESS' ;
    GREEDY_ALGOR = 'WILLCOX_EXHAUSTIVE' ;
    GREEDY_ALGOR = 'WILLCOX_SCR1' ;  % Screening criteria 1 (Astrid)
    GREEDY_ALGOR = 'NORMIDENT_JAHO';  % Minimizing norm(dPHI^T=  (PHI_T- inv(M)*PHI_G'*P))
    GREEDY_ALGOR = 'NORMIDENT_JAHO_COM';% Minimizing norm(dPHI^T). but accounting for the additional (ng-1) (stress components)
    GREEDY_ALGOR = 'CHATUR' ;
    %  GREEDY_ALGOR = 'BRUTE_FORCE1';  %% USELESS
    %  GREEDY_ALGOR = 'NORMIDENTRUE_JAHO_COM';% Similar to NORMIDENT_JAHO_COM, but it attempts to minimize
    %       (I- PHI*inv(M)*PHI_G'*P)  %%% VERY COSTLY
    %       !!!!!!!!!!!!!
    initialguess = 1:3 ;
    puntos_adicionales = 6;
end

if exist('SelectOperator') == 2
else
    addpath('jaho_routines');
end

% How to select a linearly independent set of NBASIS files in matrix PHI  ?
%

%B_JAHO
% JAHO Variables
%---------------
% DEFAULT INPUTS

CLASS_INTEGER = 'double'; % See doc uint16
ng=4;  %%% Compare with greedy algorithms
initialguess = 1:3 ;
name_SAVE_WEXH = 'DATA_JAHO/wexNEW.mat';
name_LOAD_WEXH = 'DATA_JAHO/wex.mat';
LOAD_WILLCOX_EXHAUSTIVE = 1;
nbins_WILLCOX = 50 ;
RANGE_CUMULATIVE_WILLCOX = [] ;
PARAM_STUDY = 'LEVEL_COND';
ncomb_rand_bf = 1e3;
nelemrep = 1;
INCLUDE_GAUSS_POINTS = 1 ;
LOAD_RESULT_GREEDY = 0 ;
NAME_RESULT_GREEDY = [] ;
NAME_SAVE_RESULT_GREEDY = 'DATA_JAHO/gggloc.mat';
% == CHATUR,  Chaturantabut_2009_Model_reduc.pdf
puntos_adicionales = 0 ;
CHATUR_NONREPEATED =0 ;   % New strategy, JAHO-10-may-2012
% It intends to avoid the nuisance of, later on, having to
% add new gauss points (for "refilling" the spots left by the
% repeated ones)
INDICES_GIVEN_LOC = [] ;
% for option GREEDY_ALGOR_GLO_CHOICE = 'HIERARCHICAL_cycle';---->
PARAMETER_GREEDY.LTRAIL = 1 ;   % Trailing basis vectors considered for computing the objective function
PARAMETER_GREEDY.NCOMB = 1 ;    % Number of combinations to be tested (neighbors in the topology defined by the
% objective function)
NMODES_stress =[]   ;  % Number of stress basis vectors
NMODE_disp =  []  ;   % Number of displacements basis vectors
U_left  = [] ; % Left singular vectors of the stress snapshot matrix
S_singular = [] ;% Singular values of the stress snapshot matrix
BASIS_EXPANSION = 'NO_EXPANSION';



% EXTRACTING INPUTS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varginOR = varargin ;
FdnamesInputs = {'CLASS_INTEGER','ng','initialguess','name_SAVE_WEXH',...
    'LOAD_WILLCOX_EXHAUSTIVE','nbins_WILLCOX','name_LOAD_WEXH','RANGE_CUMULATIVE_WILLCOX',...
    'ncomb_rand_bf','nelemrep','INCLUDE_GAUSS_POINTS','LOAD_RESULT_GREEDY','NAME_RESULT_GREEDY',...
    'NAME_SAVE_RESULT_GREEDY','puntos_adicionales','CHATUR_NONREPEATED','INDICES_GIVEN_LOC',...
    'PARAMETER_GREEDY','NMODES_stress','NMODE_disp','U_left','S_singular','BASIS_EXPANSION'};
AuxDATA=WriteAuxFdNamesNEW(FdnamesInputs,varginOR);
for id = 1:length(AuxDATA);
    eval(AuxDATA{id});
end
%E_JAHO



WORST_COMB = {};



%%%%%



%%%%%
%dbstop('99')
INDICES_STRESS_RECONST= [] ;  % Used only with BASIS_EXPANSION = 'SIMPLIFIED_METHOD' ;
INDICES_GAUSS_RECONST = [] ;  % (otherwise dummy)
SELECT_ELEMENT_STRESS = [] ;
%dbstop('101')


switch BASIS_EXPANSION
    case 'SIMPLIFIED_METHOD'  % JAHO, this strategy proves unefficient; delete it when "cleaning" the code
        switch GREEDY_ALGOR
            case 'HIERARCHICAL_cycle'
            otherwise
                error('BASIS_EXPANSION=SIMPLIFIED_METHOD only valid with GREEDY_ALGOR = HIERARCHICAL_cycle' )
        end
        
end


switch GREEDY_ALGOR
    case 'HIERARCHICAL_cycle'

        [REPEATED_POINTS SELEC_POINTS B_SELECTED cond_NUMB_WITHOUT INDICES_EXPAND cond_NUMB_WITH INDICES_ORIGINAL,...
            INDICES_STRESS_RECONST,INDICES_GAUSS_RECONST,SELECT_ELEMENT_STRESS] = ...
            Greedy_Hierarchical_cycle(PHI,ng,npg,ntens, PARAMETER_GREEDY,NMODES_stress,NMODE_disp,...
            U_left,S_singular,BASIS_EXPANSION);

    case 'BRUTE_FORCE1'
        PHIT = PHI' ;
        NPHIT =  sqrt(sum(PHIT.*PHIT));
        PHIT = PHIT./repmat(NPHIT,[size(PHIT,1),1]);
        COVT = abs(PHIT'*PHIT);
        [minVALUE indMIN]=min(COVT);
        %%%
        indMINSORT = sort(indMIN);
        [indUNIQ,first_element]= unique(indMINSORT,'first');
        
        repeated_element = [diff(first_element)  length(indMINSORT)-first_element(end)];
        
        figure(500)
        hold on
        xlabel('ROW')
        ylabel('Repetitions')
        title('Brute force algorithm: diagnosis')
        bar(indUNIQ,repeated_element)
        
        %%% LET US PICK UP THOSE ELEMENTS REPEATED MORE THAN ONCE
        IndRepeatedM1 = indUNIQ(repeated_element>=nelemrep);
        
        combinations = randi(length(IndRepeatedM1),[NBASIS,ncomb_rand_bf]);
        combinations = IndRepeatedM1(combinations);
        
        
        if INCLUDE_GAUSS_POINTS == 1
            % Combinations is to be "extended" to include additional stress
            % components
            %  combinations_old = combinations;
            disp(['APPENDING ROWS ...'])
            combinations = additional_rows(ng,combinations,CLASS_INTEGER);
        end
        
        disp(['COMPUTING norm(dPHIT)...'])
        h = waitbar(0,'COMPUTING norm(dPHIT)...');
        ncond = zeros(ncomb_rand_bf,1);
        for i = 1:ncomb_rand_bf
            waitbar(i/ncomb_rand_bf);
            ncond(i) = dPHIT(PHI,combinations(:,i));
            
        end
        close(h)
        
        figure(600)
        hold on
        xlabel('dPHIT')
        ylabel('Probability (%)')
        title('Result of BRUTE FORCE alg. based on PHI*PHI')
        % hist(cond_NUMB_WITHOUT_GLO,20);
        rANGEC = [0:5:200];
        bins_hist= histc(ncond,rANGEC);
        bins_hist =   bins_hist/length(ncond)*100;
        hhh=  plot(rANGEC,cumsum(bins_hist))
        
        legendLOC = ['nrep =',num2str(nelemrep),';','ncomb =',num2str(ncomb_rand_bf)];
        legend(hhh,legendLOC)

        medianCOVT = sum(COVT)/size(COVT,1);
        figure(501)
        hold on
        xlabel('ROW')
        ylabel('Average correlation factor')
        title('Brute force algorithm: diagnosis. Average correlation factor')
        plot(medianCOVT)
        
        
    case 'CHATUR'
        
        %dbstop('168')      
        [REPEATED_POINTS SELEC_POINTS B_SELECTED cond_NUMB_WITHOUT ...
            INDICES_EXPAND cond_NUMB_WITH INDICES_ORIGINAL] = ...
            ChaturGreedy(CHATUR_NONREPEATED,PHI,ng,CLASS_INTEGER,PARAM_STUDY,npg,ntens);

    case 'GIVEN_INDICES'
        INDICES_ORIGINAL  = INDICES_GIVEN_LOC ;
        % dbstop('174')
        [REPEATED_POINTS SELEC_POINTS B_SELECTED cond_NUMB_WITHOUT INDICES_EXPAND cond_NUMB_WITH] = ...
            GivenIndices(PHI,ng,CLASS_INTEGER,PARAM_STUDY,npg,ntens,INDICES_GIVEN_LOC');
        
    case 'NORMIDENT_JAHO'
        nrows = size(PHI,1);
        m = size(PHI,2);
        alias_error = zeros(nrows,1);
        for jrow = 1:nrows
            PHI_G = PHI(jrow,1);
            M = PHI_G'*PHI_G ;
            Pselect = zeros(1,nrows);
            Pselect(1,jrow) = 1 ;
            
            %     dbstop('226')
            dPHI = PHI(:,1)'-M\PHI_G'*Pselect;
            
            condnum =  norm(dPHI,'fro'); %sqrt(sum(sum(MIDENT.*MIDENT)));
            
            alias_error(jrow) =condnum;
        end
        [mindiff ind_INI] = min(alias_error);
        
        INDICES = [ind_INI];
        
        hhhh = waitbar(0,'Computing norm(dPHI^T) ...');
        
        for imod = 2:m
            waitbar(imod/m);
            
            for jrow = 1:nrows
                INDICES_loc = [INDICES jrow];
                PHI_RED = PHI(INDICES_loc,1:imod);
                PHI_C = PHI(:,1:imod);
                M = PHI_RED'*PHI_RED ;
                if rank(M) == length(INDICES_loc)
                    Pselect = SelectOperator(INDICES_loc',nrows);
                    %
                    dPHI = PHI_C'-M\PHI_RED'*Pselect;
                    condnum =  norm(dPHI,'fro'); %sqrt(sum(sum(MIDENT.*MIDENT)));
                else
                    condnum = 1e10 ;
                end
                
                
                ncond_glo(jrow) = condnum;
                
            end
            % Minimizing ncond_glo
            [cond_NUMB_WITHOUT indnext] = min(ncond_glo);
            INDICES = [INDICES indnext(1)];
            INDICES_GLO(imod,1) =indnext(1) ;
            cond_NUMB_WITHOUT_GLO(1) = cond_NUMB_WITHOUT;
        end
        
        close(hhhh);
        
        INDICES_EXPAND = additional_rows(ng,INDICES',CLASS_INTEGER); % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
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
        
        
        INDICES_EXPAND = sort(INDICES_EXPAND);
        %%%%%%%%%%%
        % How to determine the "points indicies" from INDICES_EXPAND
        %         INDICES = points_indices(INDICES_EXPAND,ng*npg);
        %         INDICES = unique(INDICES)';
        %         INDICES_B = repmat(npg*INDICES,npg,1)-repmat((npg-1:-1:0)',1,size(INDICES,2));
        [SELEC_POINTS,B_SELECTED,REPEATED_POINTS] = selected_points(ng*npg,npg,ntens,INDICES',CLASS_INTEGER);
        %%%%%%%%%%%
        
        
    case 'NORMIDENT_JAHO_COM'
        
        %
        % dbstop('292')
        nrows = size(PHI,1);
        m = size(PHI,2);
        npoints = nrows/(ng);                                % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
        alias_error = zeros(npoints,1);
        
        for ipoint = 1:npoints
            jrow_ini = (ipoint-1)*ng+1;                      % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
            jrow_fin = ipoint*ng;                            % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
            jrow = jrow_ini:jrow_fin;
            PHI_G = PHI(jrow,1);
            M = PHI_G'*PHI_G ;
            Pselect = zeros(ng,nrows);                       % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
            Pselect = SelectOperator(jrow',nrows);
            %     dbstop('226')
            dPHI = PHI(:,1)'-M\PHI_G'*Pselect;
            
            condnum =  norm(dPHI,'fro'); %sqrt(sum(sum(MIDENT.*MIDENT)));
            
            
            alias_error(ipoint) =condnum;
        end
        [mindiff ind_INI] = min(alias_error);
        
        jrow_ini = (ind_INI-1)*ng+1;                         % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
        jrow_fin = ind_INI*ng;                               % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
        jrow = jrow_ini:jrow_fin;
        INDICES = [jrow];
        
        hhhh = waitbar(0,'Computing norm(dPHI^T) ...');
        IND_TO_CHOOSE = 1:npoints ;
        %ind_ACUM = ind_INI ;
        IND_TO_CHOOSE(ind_INI) = [] ;
        
        %           ind_ACUM = zeros(m,1);
        %         Nind_ACUM = 1 ;
        % ind_ACUM(1) = ind_INI ;
        
        
        for imod = 2:m
            waitbar(imod/m);
            ncond_glo = 1e20*ones(npoints,1);
            %for jpoint = 1:npoints
            for jpoint_loc = 1:length(IND_TO_CHOOSE);
                jpoint = IND_TO_CHOOSE(jpoint_loc) ;
                jrow_ini = (jpoint-1)*ng+1;                  % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
                jrow_fin = jpoint*ng;                        % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
                jrow = jrow_ini:jrow_fin;
                INDICES_loc = [INDICES jrow];
                PHI_RED = PHI(INDICES_loc,1:imod);
                PHI_C = PHI(:,1:imod);
                M = PHI_RED'*PHI_RED ;
                if rank(M) == imod
                    Pselect = SelectOperator(INDICES_loc',nrows);
                    %
                    dPHI = PHI_C'-M\PHI_RED'*Pselect;
                    condnum =  norm(dPHI,'fro'); %sqrt(sum(sum(MIDENT.*MIDENT)));
                else
                    condnum = 1e10 ;
                end
                
                
                %                 [ndof_sigma_SC isort] = sort(ndof_sigma_SC);
                %                 npoints_sigma = points_indices(ndof_sigma_SC,ntens);
                %                 njjj = repmat([1:ntens],[1,length(npoints_sigma)/ntens]) ;
                %                 try
                %                     ndof_sigma_SC = (npoints_sigma-1)*ntens+njjj;
                %                 catch
                %                     ndof_sigma_SC = (npoints_sigma-1)*ntens+njjj';
                %
                %                 end
                ncond_glo(jpoint) = condnum;
                
            end
            % Minimizing ncond_glo
            [cond_NUMB_WITH indnext] = min(ncond_glo);
            jrow_ini = (indnext(1)-1)*ng+1;                  % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
            jrow_fin = indnext(1)*ng;                        % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
            jrow = jrow_ini:jrow_fin;
            INDICES = [INDICES jrow];
            INDICES_GLO(imod,1:ng) =jrow ;                   % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
            cond_NUMB_WITH(1) = cond_NUMB_WITH;
            ind_INI  = points_indicesCOMP(indnext(1),ng);    % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
            
            IND_TO_CHOOSE(ind_INI) = [] ;
            
            
        end
        
        cond_NUMB_WITHOUT = [] ;
        
        close(hhhh);
        
        
        INDICES_EXPAND = INDICES;
        PHI_RED = PHI(INDICES_EXPAND,:) ;
        M = PHI_RED'*PHI_RED;
        
        switch PARAM_STUDY
            case 'LEVEL_COND';
                cond_NUMB_WITH = cond(M);
            case 'NORMIDENT'
                MM = eye(size(M))-M;
                cond_NUMB_WITH = norm(MM,'fro');
            case 'NORMIDENT_JAHO'
                
                %                 Pselect = SelectOperator(INDICES_EXPAND,size(PHI,1));
                %                 %
                %                 dPHI = PHI'-M\PHI_RED'*Pselect;
                %                 condnum =  norm(dPHI,'fro'); %sqrt(sum(sum(MIDENT.*MIDENT)));
                %
                %
                %                 cond_NUMB_WITH = norm(MM,'fro');
                
        end
        
        
        %INDICES_EXPAND = sort(INDICES_EXPAND) ;
        %%%%%%%%%%%
        % How to determine the "points indicies" from INDICES_EXPAND
        %INDICES_EXPAND = sort(INDICES_EXPAND) ;
        %INDICES = points_indices(INDICES_EXPAND,ng);
        %INDICES = unique(INDICES);
        % WAIT !!!!
        % We have to redefine INDICES_EXPAND !!!!
        
        
        INDICES_EXPAND = sort(INDICES_EXPAND) ;
        %%%%%%%%%%%
        % How to determine the "points indicies" from INDICES_EXPAND
        %         INDICES = points_indices(INDICES_EXPAND,ng*npg);
        %         INDICES = unique(INDICES);
        %         INDICES_B = repmat(npg*INDICES,npg,1)-repmat((npg-1:-1:0)',1,size(INDICES,2));
        [SELEC_POINTS,B_SELECTED,REPEATED_POINTS] = selected_points(ng*npg,npg,ntens,INDICES,CLASS_INTEGER);
        % WAIT !!!!
        % We have to redefine INDICES_EXPAND !!!!
        
        %%%%%%%%%%%
        
    case 'NORMIDENTRUE_JAHO_COM'
        
        %
        nrows = size(PHI,1);
        m = size(PHI,2);
        npoints = nrows/ng;
        alias_error = zeros(npoints,1);
        IDENT = eye(nrows);
        for ipoint = 1:npoints
            jrow_ini = (ipoint-1)*ng+1;
            jrow_fin = ipoint*ng;
            jrow = jrow_ini:jrow_fin;
            PHI_G = PHI(jrow,1);
            M = PHI_G'*PHI_G ;
            
            Pselect = SelectOperator(jrow',nrows);
            
            dPHI = IDENT-PHI(:,1)*(M\PHI_G')*Pselect;
            
            condnum =  norm(dPHI,'fro'); %sqrt(sum(sum(MIDENT.*MIDENT)));
            
            
            alias_error(ipoint) =condnum;
        end
        [mindiff ind_INI] = min(alias_error);
        
        jrow_ini = (ind_INI-1)*ng+1;
        jrow_fin = ind_INI*ng;
        jrow = jrow_ini:jrow_fin;
        INDICES = [jrow];
        
        hhhh = waitbar(0,'Computing norm(dPHI_true^T) ...');
        
        for imod = 2:m
            waitbar(imod/m);
            
            for jpoint = 1:npoints
                jrow_ini = (jpoint-1)*ng+1;
                jrow_fin = jpoint*ng;
                jrow = jrow_ini:jrow_fin;
                INDICES_loc = [INDICES jrow];
                PHI_RED = PHI(INDICES_loc,1:imod);
                PHI_C = PHI(:,1:imod);
                M = PHI_RED'*PHI_RED ;
                if rank(M) == imod
                    Pselect = SelectOperator(INDICES_loc',nrows);
                    %
                    dPHI = IDENT-PHI_C*(M\PHI_RED')*Pselect;
                    condnum =  norm(dPHI,'fro'); %sqrt(sum(sum(MIDENT.*MIDENT)));
                else
                    condnum = 1e10 ;
                end
                
                
                ncond_glo(jpoint) = condnum;
                
            end
            % Minimizing ncond_glo
            [cond_NUMB_WITH indnext] = min(ncond_glo);
            jrow_ini = (indnext(1)-1)*ng+1;
            jrow_fin = indnext(1)*ng;
            jrow = jrow_ini:jrow_fin;
            INDICES = [INDICES jrow];
            INDICES_GLO(imod,1:ng) =jrow ;
            cond_NUMB_WITH(1) = cond_NUMB_WITH;
        end
        
        cond_NUMB_WITHOUT = [] ;
        
        close(hhhh);
        
        INDICES = sort(INDICES);
        INDICES_EXPAND = INDICES;
        PHI_RED = PHI(INDICES_EXPAND,:) ;
        M = PHI_RED'*PHI_RED;
        
        switch PARAM_STUDY
            case 'LEVEL_COND';
                cond_NUMB_WITH = cond(M);
            case 'NORMIDENT'
                MM = eye(size(M))-M;
                cond_NUMB_WITH = norm(MM,'fro');
            case 'NORMIDENT_JAHO'
                
                %                 Pselect = SelectOperator(INDICES_EXPAND,size(PHI,1));
                %                 %
                %                 dPHI = PHI'-M\PHI_RED'*Pselect;
                %                 condnum =  norm(dPHI,'fro'); %sqrt(sum(sum(MIDENT.*MIDENT)));
                %
                %
                %                 cond_NUMB_WITH = norm(MM,'fro');
                
        end
        
        
        %%%%%%%%%%%
        % How to determine the "points indicies" from INDICES_EXPAND
        %         INDICES = points_indices(INDICES_EXPAND,ng*npg);
        %         INDICES = unique(INDICES);
        %         INDICES_B = repmat(npg*INDICES,npg,1)-repmat((npg-1:-1:0)',1,size(INDICES,2));
        [SELEC_POINTS,B_SELECTED,REPEATED_POINTS] = selected_points(ng*npg,npg,ntens,INDICES,CLASS_INTEGER);
        %%%%%%%%%%%
        
        
    case {'WILLCOX','WILLCOX_EXHAUSTIVE','WILLCOX_INIGUESS','WILLCOX_SCR1'}
        
        % dbstop('555')
        if strcmp(GREEDY_ALGOR,'WILLCOX_EXHAUSTIVE') &&  LOAD_WILLCOX_EXHAUSTIVE == 1
            load(name_LOAD_WEXH)
            
        else
            
            
            nrows = size(PHI,1);
            switch GREEDY_ALGOR
                
                case 'WILLCOX'
                    ind_INI = 1 ;
                case 'WILLCOX_INIGUESS'
                    ind_INI = initialguess;%nrows;
                case 'WILLCOX_EXHAUSTIVE'
                    ind_INI = 1:nrows;
                    
                    
                case 'WILLCOX_SCR1'
                    % Screening criteria: the alias error
                    alias_error = zeros(nrows,1);
                    for jrow = 1:nrows
                        PHI_RED = PHI(jrow,1);
                        alias_error(jrow) = norm(1-PHI_RED'*PHI_RED);
                    end
                    [mindiff ind_INI] = min(alias_error);
                    
                    
            end
            
            m = size(PHI,2);
            INDICES = [] ;
            
            ncond_glo = zeros(length(ind_INI),1);
            % Loop over all possible placement points and evaluate M for each
            % point: choose the point that minimizes cond(M)
            
            
            
            INDICES_GLO  = zeros(m,length(ind_INI)) ;
            cond_NUMB_WITHOUT_GLO  = zeros(length(ind_INI),1);
            hhhh = waitbar(0,'Computing condition number for different initial guesses (WILLCOX greedy alg.) ...');
            
            for iini = 1:length(ind_INI)
                
                
                
                INDICES = [ind_INI(iini)];
                INDICES_GLO(1,iini)=ind_INI(iini);
                
                for imod = 2:m
                    
                    waitbar(imod/m);
                    
                    for jrow = 1:nrows
                        INDICES_loc = [INDICES jrow];
                        PHI_RED = PHI(INDICES_loc,1:imod);
                        % ncond_glo(jrow) = cond(PHI_RED'*PHI_RED);  % TOO
                        % COSTLY !!!!!!
                        ncond_glo(jrow) = 1/rcond(PHI_RED'*PHI_RED);
                    end
                    % Minimizing ncond_glo
                    [cond_NUMB_WITHOUT indnext] = min(abs(ncond_glo));
                    INDICES = [INDICES indnext(1)];
                    INDICES_GLO(imod,iini) =indnext(1) ;
                    cond_NUMB_WITHOUT_GLO(iini) = cond_NUMB_WITHOUT;
                end
                
                
                
                
                
                
            end
            
            
            close(hhhh);
            
        end
        
        % cond_NUMB_WITHOUT
        
        
        switch GREEDY_ALGOR
            
            
            case {'WILLCOX_EXHAUSTIVE','WILLCOX_INIGUESS'}
                figure(300)
                xlabel('Initial guess')
                ylabel('condition number')
                title('Result of WILLCOX greedy algorithm for different initial guesses ')
                plot(cond_NUMB_WITHOUT_GLO);
                
                figure(301)
                hold on
                xlabel('condition number')
                ylabel('Probability (%)')
                title('Result of WILLCOX greedy algorithm for different initial guesses. Histogram ')
                % hist(cond_NUMB_WITHOUT_GLO,20);
                
                
                RANGE_HIST_MAX = linspace(min(cond_NUMB_WITHOUT_GLO),max(cond_NUMB_WITHOUT_GLO),nbins_WILLCOX);
                
                
                bins_hist= histc(cond_NUMB_WITHOUT_GLO,RANGE_HIST_MAX);
                
                % Dividing by
                bins_hist =   bins_hist/length(cond_NUMB_WITHOUT_GLO)*100;
                bar(RANGE_HIST_MAX,bins_hist);
                save(name_SAVE_WEXH,'cond_NUMB_WITHOUT_GLO','INDICES_GLO')
                
                figure(302)
                hold on
                xlabel('condition number')
                ylabel('Probability (%)')
                title('Result of WILLCOX greedy algorithm for different initial guesses. CUMULATIVE histogram ')
                % hist(cond_NUMB_WITHOUT_GLO,20);
                
                bins_hist= histc(cond_NUMB_WITHOUT_GLO,RANGE_CUMULATIVE_WILLCOX);
                
                % Dividing by
                bins_hist =   bins_hist/length(cond_NUMB_WITHOUT_GLO)*100;
                bar(RANGE_CUMULATIVE_WILLCOX,cumsum(bins_hist));
                save(name_SAVE_WEXH,'cond_NUMB_WITHOUT_GLO','INDICES_GLO')
                
                
        end
        
        % Minimum of minimums
        [cond_NUMB_WITHOUT iii] = min(cond_NUMB_WITHOUT_GLO);
        INDICES=INDICES_GLO(:,iii);
        
        % Supremo
        [cond_NUMB_WITHOUT_supr iii] = max(cond_NUMB_WITHOUT_GLO);
        INDICES_supr=INDICES_GLO(:,iii);
        
        
        % if GREEDY_ADD ==1
        %% corresponding to the associated stress components
        % dbstop('625')
        INDICES_ORIGINAL = INDICES ;
        INDICES_EXPAND = additional_rows(ng,INDICES,CLASS_INTEGER);           % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
        
        PHI_RED = PHI(INDICES_EXPAND,:) ;
        M = PHI_RED'*PHI_RED;
        cond_NUMB_WITH = cond(M);
        
        switch PARAM_STUDY
            case 'LEVEL_COND';
                cond_NUMB_WITH = cond(M);
            otherwise
                MM = eye(size(M))-M;
                cond_NUMB_WITH = norm(MM,'fro');
                
                
                
                
        end
        
        INDICES_EXPAND_supr = additional_rows(ng,INDICES_supr,CLASS_INTEGER);  % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg
        
        PHI_RED = PHI(INDICES_EXPAND_supr,:) ;
        M = PHI_RED'*PHI_RED;
        cond_NUMB_WITH_supr = cond(M);
        
        
        INDICES_EXPAND = sort(INDICES_EXPAND);
        %%%%%%%%%%%
        % How to determine the "points indicies" from INDICES_EXPAND
        %         INDICES = points_indices(INDICES_EXPAND,ng*npg);
        %         INDICES = unique(INDICES)';
        %         INDICES_B = repmat(npg*INDICES,npg,1)-repmat((npg-1:-1:0)',1,size(INDICES,2));
        [SELEC_POINTS,B_SELECTED,REPEATED_POINTS] = selected_points(ng*npg,npg,ntens,INDICES,CLASS_INTEGER);
        %%%%%%%%%%%
        
        
        switch PARAM_STUDY
            case 'LEVEL_COND';
                cond_NUMB_WITHOUT_supr = cond(M);
            otherwise
                MM = eye(size(M))-M;
                cond_NUMB_WITHOUT_supr = norm(MM,'fro');
        end
        
        WORST_COMB = {cond_NUMB_WITHOUT_supr,INDICES_supr,cond_NUMB_WITH_supr,INDICES_EXPAND_supr} ;
        
        
        
        
        
end

%dbstop('675')


%%%%%
% AVOIDING REPEATED ELEMENTS: They are removed and replaces by randomly
% selected points

% JAHO: But, why the concern on removing "repeated elements" ?? why is it a
% nuisance to have repeated elements?


% % JAHO, UPGRADE 2-Mar-2012, BE CAREFUL !!!
% ntens =   ng;
%  [ndof_sigma_SC isort] = sort(INDICES);
%  npoints_sigma = points_indices(ndof_sigma_SC,ntens);
%  npoints_sigma_NU = unique(npoints_sigma);

% -------------------
% MODIFICACION PARA PONER LIMITE AL INCLUIR VARIOS MODOS DE TENSION Y NO HAY
% EL SUFICIENTE NUMERO DE PUNTOS DE GAUSS PARA INTEGRAR
% -------------------
SELEC_POINTS2 = [];
maxVal=max(SELEC_POINTS(:,1));
for iValue = 1:maxVal
   iLength = find(SELEC_POINTS(:,1)==iValue) ;
   
   iSET = SELEC_POINTS(iLength,:) ;
   PGS  = unique(iSET(:,2));   
   
   test=[repmat(iValue,length(PGS),1) PGS] ;
   SELEC_POINTS2 = [SELEC_POINTS2;test];   
end

SELEC_POINTS = SELEC_POINTS2 ;
B_SELECTED = unique(B_SELECTED) ;
REPEATED_POINTS = 0;
% -------------------


ADPOINT_TOT=[];
if REPEATED_POINTS>0
    
    nrows = size(PHI,1);
    m = size(PHI,2);
    npoints = nrows/(ng*npg);
    selec_test=SELEC_POINTS(:,1);
    INDADD=[];
    
    while (length(selec_test)<m)
        % while (length(INDADD)<REPEATED_POINTS)
        
        add_elements = m-length(selec_test);
        % add_elements = REPEATED_POINTS-length(INDADD);
        
        VECT_ELEMENTS = [1:npoints]';
        VECT_ELEMENTS(unique(selec_test)) = [];
        % VECT_ELEMENTS(selec_test) = [];
        
        combinations_adicional = randi(length(VECT_ELEMENTS),[add_elements,1]);
        INDADD = VECT_ELEMENTS(combinations_adicional);
        
        INDADD = unique(INDADD);
        
        selec_test =[selec_test;INDADD];
        ADPOINT_TOT = [ADPOINT_TOT;INDADD];
        
    end
    
    % pgauss = randi(npg,size(INDADD,1));
    % INDADD = [INDADD pgauss(:,1)];
    % SELEC_POINTS = [SELEC_POINTS;INDADD];
    pgauss = randi(npg,size(ADPOINT_TOT,1));
    ADPOINT_TOT = [ADPOINT_TOT pgauss(:,1)];
    SELEC_POINTS = [SELEC_POINTS;ADPOINT_TOT];
    SELEC_POINTS = sortrows(SELEC_POINTS,1);
    B_SELECTED = (SELEC_POINTS(:,1)-1)*npg + SELEC_POINTS(:,2);
    
    % GENERALIZED TENSOR STRESS COMPONENTES
    V_PREV=zeros(ntens,m);
    for i = 1:m
        V_PREV(:,i) = [(B_SELECTED(i)-1)*ntens+1:B_SELECTED(i)*ntens]';
    end
    INDICES_EXPAND = V_PREV(:);  %additional_rows(ng,B_SELECTED,CLASS_INTEGER);
    % INDICES_EXPAND = sort(INDICES_EXPAND);
    
end
%end
