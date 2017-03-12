function  [REPEATED_POINTS B_SELECTED SELEC_POINTS INDICES_EXPAND INDICES_ORIGINAL cond_NUMB_WITH cond_NUMB_WITHOUT] = AdditionalPointsSelection(INDICES_FINAL,INDICES_FINAL_GAUSS,PHI,NMODE_disp,PARAMETER_GREEDY,PHI_EXPANDED,...
    ng,NMODES_stress,npg)
% See Greedy_Hierarchical_cycle


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

INDICES_K = INDICES_FINAL' ;
GAUSS_POINTS_K = INDICES_FINAL_GAUSS ; % Initial set of gauss points

INDICES_OPT = INDICES_K ;
GAUSS_POINTS_OPT = GAUSS_POINTS_K;


PHI_K = PHI ;

disp('Determining "additional" sample points...')
switch PARAMETER_GREEDY.CRITERION_ADDPOINTS
    case 'COND_GREEDY'
        
        
        switch PARAMETER_GREEDY.CRITERION_NUMBER_ADDPOINTS
            case {'ALL','FIXED'}
                
                switch  PARAMETER_GREEDY.CRITERION_NUMBER_ADDPOINTS
                    case 'ALL'
                        NMODE_disp_LOC = NMODE_disp;
                    case 'FIXED'
                        NMODE_disp_LOC = PARAMETER_GREEDY.NUMBER_ADDITIONAL;
                end
                
                %NCOND_PHI_EXPANDED = [] ;
                %Mexpanded = PHI_EXPANDED(INDICES_K,:)'*PHI_EXPANDED(INDICES_K,:);
                %NCOND_PHI_EXPANDED(end+1) = cond(Mexpanded) ;
                
                
                ngauss  = size(PHI,1)/ng ;
                %  dbstop('107')
                for imodeB = 1:NMODE_disp_LOC
                    % Loop over "B-matrix" modes
                    
                    % New set of basis vectors
                    PHI_K = [PHI PHI_EXPANDED(:,NMODES_stress+1:NMODES_stress+imodeB)] ; %
                    %  ncond_glo = zeros(ngauss,1);
                    % Loop over gauss points
                    min_NCOND = 10000000000000000000000000000000 ;
                    for igauss = 1:ngauss
                        % Component indexes associated to the igauss-th gauss point
                        iii_glo =  (igauss-1)*ng+1;
                        fff_glo =   igauss*ng;
                        % Testing a New set of indexes
                        INDICES_test = [INDICES_K iii_glo:fff_glo];
                        % Constructing the gappy matrix PHI_G
                        PHI_G =    PHI_K(INDICES_test,:);
                        % Square Mgappy matrix
                        Mgappy = PHI_G'*PHI_G ;
                        % Condition number (estimated via rcond; it is cheaper than cond())
                        ncond_Mgappy =   1/rcond(Mgappy) ;
                        
                        if ncond_Mgappy<min_NCOND
                            min_NCOND = ncond_Mgappy ;
                            ind_g = igauss ;
                        end
                        
                        % ncond_glo(igauss) = ncond_Mgappy ;
                    end
                    
                    
                    
                    %  [mmm ind_g] =  min(ncond_glo) ;
                    GAUSS_POINTS_K = [GAUSS_POINTS_K ind_g] ;
                    iii_glo =  (ind_g-1)*ng+1;
                    fff_glo =   ind_g*ng;
                    INDICES_K = [INDICES_K iii_glo:fff_glo] ;
                    
                    % Mexpanded = PHI_EXPANDED(INDICES_K,:)'*PHI_EXPANDED(INDICES_K,:);
                    % NCOND_PHI_EXPANDED(end+1) = cond(Mexpanded) ;
                    
                    
                end
                
                Mexpanded = PHI_EXPANDED(INDICES_K,:)'*PHI_EXPANDED(INDICES_K,:);
                NCOND_PHI_EXPANDED_FINAL  = cond(Mexpanded) ;
                disp('*************************')
                disp('CONDITION NUMBER OF GAPPY EXPANDED MATRIX')
                disp(['cond(Mexpand)=',num2str(NCOND_PHI_EXPANDED_FINAL)])
                disp('*************************')
                
                %                 figure(10)
                %                 hold on
                %                 xlabel('Number of appended points')
                %                 ylabel('Condition number of Mgappy_expanded')
                %                 h = plot(0:NMODE_disp,NCOND_PHI_EXPANDED,'Color',coloresM(1,:),'Marker',MarkerM{1})
                %                 legend(h,['NU=NS=',num2str(NMODE_disp)])
                %
                %
                %                 INDICES_K_set = INDICES_FINAL
                %                 PHI_G = PHI(INDICES_K_set,:);
                %                 ReconsGAPPY = PHI_G*inv(PHI_G'*PHI_G)*PHI_G' ;
                %
                %                 UB = PHI_EXPANDED(:,NMODES_stress+1:end);
                %                 nind = length(INDICES_K_set);
                %                 % Does the inverse of UBgappy'*(I-ReconsGAPPY)*UBgappy exist  ???
                %                 UBgappy = UB(INDICES_K_set,:);
                %                 Ngappy =   UBgappy'*(eye(nind)-ReconsGAPPY)*UBgappy ;
                
            case 'NCOND'
                
                
                %  dbstop('178')
                Mexpanded = PHI_EXPANDED(INDICES_K,:)'*PHI_EXPANDED(INDICES_K,:);
                NCOND_PHI_EXPANDED_FINAL = cond(Mexpanded) ;
                imodeB = 1 ;
                while NCOND_PHI_EXPANDED_FINAL> PARAMETER_GREEDY.NCOND_THRESHOLD & imodeB<=NMODE_disp
                    
                    
                    ngauss  = size(PHI,1)/ng ;
                    %dbstop('42')
                    
                    % Loop over "B-matrix" modes
                    
                    % New set of basis vectors
                    PHI_K = [PHI PHI_EXPANDED(:,NMODES_stress+1:NMODES_stress+imodeB)] ; %
                    %  ncond_glo = zeros(ngauss,1);
                    % Loop over gauss points
                    min_NCOND = 10000000000000000000000000000000 ;
                    for igauss = 1:ngauss
                        
                        % Component indexes associated to the igauss-th gauss point
                        iii_glo =  (igauss-1)*ng+1;
                        fff_glo =   igauss*ng;
                        % Testing a New set of indexes
                        INDICES_test = [INDICES_K iii_glo:fff_glo];
                        % Constructing the gappy matrix PHI_G
                        PHI_G =    PHI_K(INDICES_test,:);
                        % Square Mgappy matrix
                        Mgappy = PHI_G'*PHI_G ;
                        % Condition number (estimated via rcond; it is cheaper than cond())
                        ncond_Mgappy =   1/rcond(Mgappy) ;
                        
                        if ncond_Mgappy<min_NCOND
                            min_NCOND = ncond_Mgappy ;
                            ind_g = igauss ;
                        end
                        
                        % ncond_glo(igauss) = ncond_Mgappy ;
                    end
                    
                    
                    
                    %  [mmm ind_g] =  min(ncond_glo) ;
                    GAUSS_POINTS_K = [GAUSS_POINTS_K ind_g] ;
                    iii_glo =  (ind_g-1)*ng+1;
                    fff_glo =   ind_g*ng;
                    INDICES_K = [INDICES_K iii_glo:fff_glo] ;
                    
                    % Mexpanded = PHI_EXPANDED(INDICES_K,:)'*PHI_EXPANDED(INDICES_K,:);
                    % NCOND_PHI_EXPANDED(end+1) = cond(Mexpanded) ;
                    
                    
                    Mexpanded = PHI_EXPANDED(INDICES_K,:)'*PHI_EXPANDED(INDICES_K,:);
                    NCOND_PHI_EXPANDED_FINAL = cond(Mexpanded) ;
                    
                    imodeB = imodeB + 1; % :NMODE_disp_LOC
                end
                
                disp('*************************')
                disp('CONDITION NUMBER OF GAPPY EXPANDED MATRIX')
                disp(['cond(Mexpand)=',num2str(NCOND_PHI_EXPANDED_FINAL)])
                disp('*************************')
                
                %                 figure(10)
                %                 hold on
                %                 xlabel('Number of appended points')
                %                 ylabel('Condition number of Mgappy_expanded')
                %                 h = plot(0:NMODE_disp,NCOND_PHI_EXPANDED,'Color',coloresM(1,:),'Marker',MarkerM{1})
                %                 legend(h,['NU=NS=',num2str(NMODE_disp)])
                %
                %
                %                 INDICES_K_set = INDICES_FINAL
                %                 PHI_G = PHI(INDICES_K_set,:);
                %                 ReconsGAPPY = PHI_G*inv(PHI_G'*PHI_G)*PHI_G' ;
                %
                %                 UB = PHI_EXPANDED(:,NMODES_stress+1:end);
                %                 nind = length(INDICES_K_set);
                %                 % Does the inverse of UBgappy'*(I-ReconsGAPPY)*UBgappy exist  ???
                %                 UBgappy = UB(INDICES_K_set,:);
                %                 Ngappy =
                %                 UBgappy'*(eye(nind)-ReconsGAPPY)*UBgappy
                %                 ;
                
                
                
                
            otherwise
                
                
                
                error('Option not implemented')
        end
        
    case 'DEIM'
        
        %  DEIM --> New points are selected following the same criterion
        % ----------
        % used in DEIM (CHATUR. algorithm) algorithm ....
        [NCOND_PHI_EXPANDED_FINAL,GAUSS_POINTS_K,INDICES_K,Mexpand] = ...
            AdditionalPointsDEIM(INDICES_K,PHI_EXPANDED,NMODE_disp,ng,NMODES_stress,GAUSS_POINTS_K) ;
        
    otherwise
        error('Option not implemented')
        
        
end



% OUTPUTS
REPEATED_POINTS = 0 ;  % In this refined implementation, this output variable is dummy
B_SELECTED = GAUSS_POINTS_K ; % B_SELECTED  -->    Certainly this is a misleading name; it is simply
 % the vector that contains the indexes of the selected quadrature points
%SELEC_POINTS(i) => [number of finite element   number of gauss point within the finite element]
IGGG =  mod(GAUSS_POINTS_K,npg) ;
[ind1] = find(IGGG==0)  ;
if ~isempty(ind1)
    IGGG(ind1) = npg ;
end

elements = (GAUSS_POINTS_K-IGGG)/npg+1;

SELEC_POINTS = [elements', IGGG'];


%cond_NUMB_WITHOUT INDICES_EXPAND cond_NUMB_WITH INDICES_ORIGINAL

INDICES_EXPAND = INDICES_K ;
INDICES_ORIGINAL = INDICES_K ;
cond_NUMB_WITH = NCOND_PHI_EXPANDED_FINAL ;
cond_NUMB_WITHOUT = cond_NUMB_WITH ;
