function [fobj_total FOUND_ITEM MIN_OBJ_FUN IND_WIN IND_GAUSS] = LOOP_loc_JAHO_OBJ_newMULTI(INDICES_gauss,nrows,...
    imod,PHI,S,V,max_COND,FOUND_ITEM,fobj_total,ng)

% if nargin == 9
%     WITH_VTRAIL = 0 ;
% end

n_total = size(PHI,2);
try
S_trail = S(imod+1:n_total,imod+1:n_total) ;
catch
    error('Set NCOMB to  a lower value')
    
end
try
    V_trail = V(:,imod+1:n_total);
catch
    V_trail = [] ;
end

PSI_S = PHI(:,imod+1:end)*S_trail ;

MIN_OBJ_FUN = 444444666666666666666666666666666666666;
IND_WIN = [];



% INDICES_gauss  --->  INDICES
%INDICES_EXPAND = additional_rows(ng,INDICES_FINAL,'uint16');
%                        INDICES = reshape(INDICES_EXPAND,[size(INDICES_EXPAND,1)*size(INDICES_EXPAND,2) 1]);

if ~isempty(INDICES_gauss)
    INDICES = zeros(1,length(INDICES_gauss)*ng) ;
    % This can be optimized !!
    for i = 1:length(INDICES_gauss)
        igauss = INDICES_gauss( i );
        iii_glo =  (igauss-1)*ng+1;
        fff_glo =   igauss*ng;
        
        igauss = i;
        iii =  (igauss-1)*ng+1;
        fff =   igauss*ng;
        
        INDICES(iii:fff) = iii_glo:fff_glo ;
    end
    
else
    INDICES = [] ;
end


for igauss = 1:nrows/ng
    % dbstop('22')
    iii =  (igauss-1)*ng+1;
    fff =   igauss*ng;
    jrow = iii:fff; % Stored indexes
    
    INDICES_loc = [INDICES jrow];
    % Objective function
    %PHI_RED = PHI(INDICES_loc,1:imod);
    %ncond_glo(jrow) = cond(PHI_RED'*PHI_RED);
    %
    PHI_RED = PHI(INDICES_loc,1:imod); % PHI_GAPPY
    % IMPORTANT. The trailing singular vectors at stage "imod"
    % are the rest of vectors of PHI = [PHI_before PSI_before(:,1:Ltrail)] ....
    %PSI_RED = PHI(INDICES_loc,imod+1:end);  % PSI_GAPPY
    PSI_S_RED = PSI_S(INDICES_loc,:);
    
    % S_trail : trailing singular values
    
    
    M_RED = PHI_RED'*PHI_RED;
    % We compute the condition number of M_RED
    cond_M = cond(M_RED); % Too slow
    cond_M = 1/rcond(M_RED);
    
    % Is smaller than max_COND ?
    if cond_M < max_COND
        % Objective function
        %         if  jrow ==500
        %         disp('')
        %         end
        %dbstop('46')
        %if WITH_VTRAIL  == 1
        %     normOBJ_FUNCT = norm(M_RED\(PHI_RED'*PSI_S_RED)*V_trail','fro');
        %else
        
        normOBJ_FUNCT = norm(M_RED\(PHI_RED'*PSI_S_RED),'fro');
        %end
        
        
        
        fobj_total(igauss) = normOBJ_FUNCT ;
        %             if normOBJ_FUNCT<MIN_OBJ_FUN
        %                 MIN_OBJ_FUN = normOBJ_FUNCT ;
        %                 IND_WIN = jrow ;
        %             end
        FOUND_ITEM = 1 ;
        
        if normOBJ_FUNCT<MIN_OBJ_FUN
            MIN_OBJ_FUN = normOBJ_FUNCT ;
            IND_WIN = jrow ;
            IND_GAUSS = igauss ;
        end
    else
        fobj_total(igauss) = 10000000000000000000000000;
    end
    
end
