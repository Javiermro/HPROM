function [REPEATED_POINTS SELEC_POINTS B_SELECTED cond_NUMB_WITHOUT INDICES_EXPAND cond_NUMB_WITH INDICES_ORIGINAL] = ChaturGreedy(CHATUR_NONREPEATED,PHI,ng,CLASS_INTEGER,PARAM_STUDY,npg,ntens)
% See greedy_algorithms.m
if   CHATUR_NONREPEATED == 0
    % Default,classical strategy
    
    Uk = PHI(:,1) ;
    m = size(PHI,2);
    INDICES = zeros(m,1) ;
    k = 1 ;
    
    [rho ind] = max(abs(PHI(:,k)));
    INDICES(k) = ind ;
    
    %dbstop('156')
    %  hhhh = waitbar(0,'Computing norm(max(abs(r))) ...');
    
    for k = 2:m
        %     waitbar(k/m);
        ind_k = INDICES(1:(k-1));
        %U = PHI(:,1:(k-1));
        U_red = Uk(ind_k,:);
        c = U_red\PHI(ind_k,k);
        r = PHI(:,k)-Uk*c ;
        [  rho_max   ind_kp1   ] = max(abs(r));
        Uk = [Uk PHI(:,k)];
        INDICES(k) = ind_kp1 ;
    end
    
    
    
    
elseif CHATUR_NONREPEATED == 1
    
    % New strategy, JAHO-10-may-2012
    % It intends to avoid the nuisance of, later on, having to
    % add new gauss points (for "refilling" the spots left by the
    % repeated ones)
    
    % the total number of modes must be less than the total amount of 
    % gauss points, otherwise, it's not possible to assign an equal number of 
    % integation points that number of given modes.
    if size(PHI,1)>(size(PHI,2)/ng)
        warning('You must fix a lower number of total modes in this greedy strategy')
    end
    
    % dbstop('160')
    Uk = PHI(:,1) ;
    m = size(PHI,2);
    INDICES = zeros(m,1) ;
    k = 1 ;
    
    [rho ind] = max(abs(PHI(:,k)));
    INDICES(k) = ind ;
    
    %dbstop('156')
    %  hhhh = waitbar(0,'Computing norm(max(abs(r))) ...');
    GAUSS_POINTS = zeros(m,1) ;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Can one ensure that this "gauss point" is not repeated ???
    %dbstop('175')
    igauss =  mod(INDICES(k),ng) ;
    [ind1] = find(igauss==0)  ;
    if ~isempty(ind1)
        igauss(ind1) = ng ;
    end
    % disp('... calculating points ...')
    GAUSS_POINTS(k) = (INDICES(k)-igauss)/ng+1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    for k =2:m
        %     waitbar(k/m);
        ind_k = INDICES(1:(k-1));
        %U = PHI(:,1:(k-1));
        U_red = Uk(ind_k,:);
        c = U_red\PHI(ind_k,k);
        r = PHI(:,k)-Uk*c ;
        
        
        POINT_FOUND = 0;
        [  rho_max   ind_kp1   ] = max(abs(r));
        while  POINT_FOUND ==0
            %if k==153
            %    aaa=2; 
            %end
            
            %INDICES(k) = ind_kp1 ;
            
            % Can one ensure that this "gauss point" is not repeated ???
            igauss =  mod(ind_kp1,ng) ;
            [ind1] = find(igauss==0)  ;
            if ~isempty(ind1)
                igauss(ind1) = ng ;
            end
            %disp('... calculating points ...')
            NEW_POINT = (ind_kp1-igauss)/ng+1;
            if any(NEW_POINT == GAUSS_POINTS(1:k-1))
                %repeated
                % dbstop('213')
                r(ind_kp1) = 0 ;
                [  rho_max   ind_kp1   ] = max(abs(r));
            else
                POINT_FOUND = 1 ;
                
            end
        end
        
        INDICES(k) = ind_kp1 ;
        Uk = [Uk PHI(:,k)];
        GAUSS_POINTS(k)=NEW_POINT;
        
        % break the process when the total number of pg's have been assigned
%         if k>(size(PHI,1)/ng)
%             break
%         end
        
    end
    
end

INDICES_ORIGINAL = INDICES ;
PHI_RED = PHI(INDICES,:) ;
M = PHI_RED'*PHI_RED;
cond_NUMB_WITHOUT = cond(M);

%% But it is instructive (and even necessary) to inquire about the
%% condition number of the reduced PHI matrix consnstructed by adding the files
%% corresponding to the associated stress components
%%%%%INDICES_EXPAND = additional_rows(ng,INDICES,CLASS_INTEGER); % now is fixed the seleccion on Gauss Points, otherwise ng=ng*npg

%%%%%INDICES_EXPAND = unique(INDICES_EXPAND);
INDICES_EXPAND = unique(INDICES);

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

%%%%%[SELEC_POINTS,B_SELECTED,REPEATED_POINTS] = selected_points(ng*npg,npg,ntens,INDICES,CLASS_INTEGER);
[SELEC_POINTS,B_SELECTED,REPEATED_POINTS] = selected_points(npg,npg,1,INDICES,CLASS_INTEGER);
%%%%%%%%%%%


