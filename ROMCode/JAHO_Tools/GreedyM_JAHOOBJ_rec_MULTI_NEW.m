function [INDICES_FINAL MIN_TOTAL INDICES_GLO  min_comb INDICES_FINAL_GAUSS]= GreedyM_JAHOOBJ_rec_MULTI_NEW(PHI,K,PSI,S,L_trail,V,ncomb,....
    PARAMETERS)
% NOVEL PROPOSAL
% -------------
% Greedy minimization of function
% inv(M)*PHI_gappy^T*PSI_gappy_L*S_trailing_L
% where PSI_gappy_L = PSI_gappy(:,1:L)

% this is quite a convoluted strategy; it seems to believe it will work ...
%

%dbstop('13')
if nargin == 0
    load('tmp.mat')
   % ncomb = 3 ;
    % PARAMETERS.LOAD_PREVIOUS = 0 ;
    %L_trail = 30 ;
end

try
    WITH_VTRAIL = PARAMETERS.WITH_VTRAIL;
catch
    WITH_VTRAIL = 0 ;
end
%
% PARAMETERS.JAHO_OBJ_MULT.LOAD_PREVIOUS = 0 ;
% PARAMETERS.JAHO_OBJ_MULT.NAME_WS_IDENT = 'JAHOobjMULT_pr_' ;

ng = PARAMETERS.ng ;


NAME_LOAD = ['GREEDY_JAHO/',PARAMETERS.NAME_WS_IDENT,'L_',num2str(L_trail),'_M',num2str(size(PHI,2)),'_NC',num2str(ncomb),'.mat'];
STEP_BY_STEP = 0;
%dbstop('32')
if exist(NAME_LOAD)
    %dbstop('34')
    load(NAME_LOAD,'INDICES_FINAL')
    if ~exist('INDICES_FINAL') || isempty(INDICES_FINAL)
        STEP_BY_STEP = 1;
    else
        if   isempty(INDICES_FINAL) ; STEP_BY_STEP = 1; end
    end
    
    
end

INDICES_FINAL = [] ; MIN_TOTAL=[] ;INDICES_GLO = [] ; min_comb = [];

if (PARAMETERS.LOAD_PREVIOUS ~= 1 || exist(NAME_LOAD) == 0) || STEP_BY_STEP==1
    
    
    nrows = size(PHI,1);
    max_COND = 500000 ; % To ensure exisence of inv(V)
    PSI_before = PSI ;
    PHI_before = PHI ;
    
%dbstop('58')
    m = size(PHI,2);  % Number of leading modes
    L_trail =min(L_trail,size(PSI,2)+size(PHI,2));  % Number of trailing basis vectors
    
    %     dbstop('57')
    %
    % For computing the objective function, we presume that
    PHI = [PHI PSI(:,1:L_trail)];
    
    
    
    
    
    %
    
    INDICES = [] ;
    INDICES_GLO = zeros(ncomb,m);
    
    
    
    candidate_indices_next = [] ;
    
    % LOOP OVER NUMBER OF "MODES" (NS)
    for imod = 1:m
        disp(['-------------------------'])
        disp(['Mode m=',num2str(imod)])
        disp(['-------------------------'])
        % Loop over the number of modes
        % ------------------------------
        % waitbar(imod/m);
        
        NAME_LOAD_LOC = ['GREEDY_JAHO/',PARAMETERS.NAME_WS_IDENT,'L_',num2str(L_trail),'_M',num2str(imod),'_NC',num2str(ncomb),'LOC.mat'];
        
        %dbstop('90')
        if PARAMETERS.LOAD_PREVIOUS == 0 || exist(NAME_LOAD_LOC) == 0
            
            INDICES_GLO;
            
            % INDICES_GLO is a matrix of NCOMB rows and NMODES (m) columns
            % In cycle i, (1<i<=NMODES), the algorithm greedily optimizes
            % the objective function for each combination.  For instance,
            % for m=3, NCOMB = 10 and imod=3, the contents of matrix
            % INDICES_GLO at this point of the code might be something
            % like:
            
            %             INDICES_GLO =[
            %
            %          661       13183           0
            %          663        7946           0
            %          508       15632           0
            %          662       13183           0
            %          664       15630           0
            %          683       13183           0
            %          681       13565           0
            %          684       13183           0
            %          682        1838           0
            %         4749       13183           0]
            
            % It the sequel, the algorithm will fill the third column of  INDICES_GLO
            % by finding the "greedy" optimum for each combination. The "global"
            % minimizer will be computed at the end of the ensuing "combination" loop by determining
            % the combination that gives the minima value of the objective function.
            
            
            
            MIN_OBJ_FUN = 1e20 ;
            IND_WIN =[] ;
            FOUND_ITEM = 0 ;
            min_comb = zeros(ncomb,1);  % ---> EXIT
            fobj_total = zeros(1,nrows/ng);
            % size of the number of combinations
            
            if imod == 1
                size_comb = 1 ;
                %
            else
                size_comb = ncomb ;
            end
            
            
            MIN_PER_COMB =44444444444444444;
            for   icomb = 1:size_comb
                
                disp(['Alternative path icomb=',num2str(icomb)])
                % One has to perform size_comb cycles (except for the very first mode)
                
                INDICES = INDICES_GLO(icomb,1:imod-1);
                
                
                
                if ~any(INDICES>0)
                    INDICES=[] ;
                end
                
                % Single greedy step for the set of row indices "INDICES"
               % dbstop('152')
                [fobj_total FOUND_ITEM MIN_OBJ_FUN indices_min IND_GAUSS] = LOOP_loc_JAHO_OBJ_newMULTI(...
                    INDICES,nrows,imod,PHI,S,V,max_COND,FOUND_ITEM,fobj_total,ng);
                
                if isempty(indices_min)
                    % The search was not succesful (cond(M)>max_COND)
                    disp('Combination not found -----')
                    %             if icomb> 1
                    %                 jrow_min = INDICES_GLO( icomb-1,imod)
                    %             end
                    
                    % This combination will be rejected anyway (it yields a very
                    % large MIN_OBJ_FUN ), thus, it is irrelevant the corresponding index
                    %
                    indices_min = [87] ;
                    MIN_OBJ_FUN = fobj_total(1);
                end
                
                % Minimizing ncond_glo
                if FOUND_ITEM == 0
                    error('No index was found ....')
                    
                end
                
                % Which is the winnning index ??? The selected point.
                
                min_comb(icomb) = MIN_OBJ_FUN ;
                
                if imod>1
                    
                    % In the event of repeated indices ...
                    if icomb >1
                        
                        % [dummy indices_sort] = sort(fobj_total);
                        ibuscar = 1 ;
                        SALIR = 0 ;
                        while SALIR == 0 && ibuscar<=length(fobj_total)
                            IND_FIND = find(INDICES_GLO(1:icomb-1,imod) == IND_GAUSS);
                            
                            if isempty(IND_FIND)
                                INDICES_GLO(icomb,imod) = IND_GAUSS;
                                SALIR = 1 ;
                            else
                                % Repeated index. However, this does not
                                % necessarily imply that the index is to be
                                % rejected: min(min_comb(1:icomb-1)) -MIN_OBJ_FUN
                                
                                FUN_COMP = min(min_comb(1:icomb-1)) -MIN_OBJ_FUN;
                                if  FUN_COMP >1e-10
                                    %This means that the new index represent an
                                    %improvement
                                    INDICES_GLO(icomb,imod) = IND_GAUSS;
                                    SALIR = 1 ;
                                    
                                    
                                else
                                    disp('Repeated index (with no improvement)')
                                    ibuscar = ibuscar+1;  % What is this for ?== JAHO, Jul 6,2012
                                    if ibuscar == 2  % For selecting the next best point ??
                                        [dummy indices_sort] = sort(fobj_total);
                                    end
                                    IND_GAUSS = indices_sort(ibuscar);
                                    min_comb(icomb) = fobj_total(IND_GAUSS) ;
                                end
                                
                                
                                
                            end
                            
                            
                        end
                        
                        
                        %                 if  INDICES_GLO(icomb-1,imod) == jrow_min;
                        %                     % Repeated
                        %                     [dummy indices_sort] = sort(fobj_total);
                        %                     SALIR = 0 ;
                        %                     ibuscar = 1;
                        %                     while  indices_sort(ibuscar)== INDICES_GLO(icomb-1,imod) & ibuscar<=length(fobj_total)
                        %                         ibuscar = ibuscar + 1;
                        %                     end
                        %                 else
                        %                     INDICES_GLO(icomb,imod) = jrow_min;
                        %                 end
                    else
                        INDICES_GLO(icomb,imod) = IND_GAUSS;
                    end
                    
                else
                    % Which are the winning indices when imod = 1
                    % We do not yet
                    [dummy indices_sort] = sort(fobj_total);
                    
                    candidate_indices_next = indices_sort(1:ncomb);
                    INDICES_GLO(:,1) = candidate_indices_next';
                    min_comb = fobj_total(candidate_indices_next);
                end
                
                %%% STORING ....
                INDICES_GLO_MODE = INDICES_GLO(:,imod);
                min_comb_MODE = min_comb;
                save(NAME_LOAD_LOC,'INDICES_GLO_MODE','min_comb_MODE')
                %%%
                
                MIN_PER_COMB = min(MIN_PER_COMB,min_comb(icomb));
                 if MIN_PER_COMB == 0
                    dbstop('259')
                    warning('This might be an error')
                end
                disp(['Absolute minimum = ',num2str(MIN_PER_COMB)]) ;
                disp(['min_comb(icomb) = ',num2str(min_comb(icomb))]) ;
                
                
            end
            
        else
            
            
            load(NAME_LOAD_LOC,'INDICES_GLO_MODE','min_comb_MODE')
            INDICES_GLO(:,imod) = INDICES_GLO_MODE;
            min_comb = min_comb_MODE ;
        end
        
        
    end
    
    
    
    %INDICES_GLO
    %min_comb
    
    
    
    [MIN_TOTAL IND_COMB] = min(min_comb);
    INDICES_FINAL_GAUSS = INDICES_GLO(IND_COMB,:);
    INDICES_FINAL = zeros(length(INDICES_FINAL_GAUSS)*ng,1);
    %MIN_TOTAL
    
    for i = 1:length(INDICES_FINAL_GAUSS)
        igauss = INDICES_FINAL_GAUSS( i );
        iii_glo =  (igauss-1)*ng+1;
        fff_glo =   igauss*ng;
        
        igauss = i;
        iii =  (igauss-1)*ng+1;
        fff =   igauss*ng;
        
        INDICES_FINAL(iii:fff) = iii_glo:fff_glo ;
    end
    
    
    save(NAME_LOAD,'INDICES_FINAL','MIN_TOTAL','INDICES_GLO','min_comb')
    
    
    
    % close(hhhh);
    
    
else
    
    
    
    load(NAME_LOAD);
    
    
    
    
end

