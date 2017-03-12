function [z,w,errorGLO,J,b] = OptimizedQuadraturePARTITIONED(SNAPfNW,ndom,nperm,W,DATA,TOL,DATA_CUBATURE)
%
% Partitioned version of the optimized quadrature algorithm
% Joaquín A. Hernández, January 8-th - 2016
% See PaperCubature.pdf
%---------------------------------------
% INPUTS:
% ------
% SNAPfNW --> Matrix of snapshots
% ndom  -->  Number of partitions
% nperm -->  Number of permutations for each partition
% W --> Vector of finite element weights
% TOL --> Tolerance for the optimization process
% DATA --> Miscellaneous data
%     Where to put the finite element weights in the J matrix
% %   DATA.vectWEIGHTS_inJ = 'NONE' ; %'NO'; %'SQUARE' ;  ;'ALL';   'SQUARE' ;   ;
%     Strategy upon encounter of negative weights
%     DATA.NNEGLSQ_ONLY_WHEN_NECESSARY = 0;
% END INPUTS
% OUTPUTS
% z --> Set of chosen points (their indices)
% w --> Corresponding weights
% ----------------------------------------------------
%dbstop('25')
if nargin == 0
    load('tmp.mat')
    DATA_CUBATURE.unitaryBvector = 0 ;
end

DATA = DefaultField(DATA,'npoints',[]);
TIMETOTAL = tic ;
DATA = DefaultField(DATA,'zINI',[]) ;
DATA = DefaultField(DATA,'DivideByNorm',0) ;
DATA.unitaryBvector = DATA_CUBATURE.unitaryBvector ;
SNAPfNW_glo = cell(ndom*nperm,1) ;

%nmodes_dom = size(SNAPfNW,2)/ndom ;
z_glo = cell(ndom*nperm,1) ; w_glo = cell(ndom*nperm,1) ;
iacum = 0;

if DATA.STANDARD_APPROACH == 1
    ndom = 1 ;
end
J = [] ; b= [] ;
if  ndom > 1
    switch  DATA.vectWEIGHTS_inJ
        case  'NO'
            error('Incompatible option')
    end
    
    %
    disp('Global J and b'  )
    disp('--------------------------------------------------------------')
    [J,b,DATA] = Jmin_bmin_calc(SNAPfNW,W,DATA) ;    %
    
    nmodes = size(J,1)-1 ;
    
    
    if mod(nmodes,ndom) == 0
        nmodes_dom = (nmodes/ndom)*ones(nmodes,1) ;
    else
        nmodes_dom = floor(nmodes/ndom)*ones(ndom,1) ; % Number of points per domain
        nmodes_ast = nmodes-sum(nmodes_dom) ;
        nmodes_dom(end) = nmodes_dom(end) + nmodes_ast;
    end
    iniACUM = 1 ;
    for idom = 1:ndom
        iacum = iacum + 1;
        disp('--------------------------------')
        disp(['Set of vectors = ',num2str(iacum)])
        disp('--------------------------------')
        if idom == 1
            ini = 1;
        else
            ini = sum(nmodes_dom(idom-1))   + 1;
        end
        
        ifin = sum(nmodes_dom(idom))   ;
        DATA.NOSHOWITERS = 1 ;
        m=1e10 ;
        
        [z_glo{iacum},w_glo{iacum},errorGLO,~,~]= OptimizedQuadratureMultiREST([],W,TOL,m,DATA,[],...
            [J(ini:ifin,:); J(end,:)],[b(ini:ifin); b(end) ],DATA_CUBATURE);
        
        for iPERM = 2:nperm
            iacum = iacum + 1;
            disp('--------------------------------')
            disp(['Set of vectors = ',num2str(iacum)])
            disp('--------------------------------')
            COL = randperm(nmodes,nmodes_dom(idom));
            [z_glo{iacum},w_glo{iacum},errorGLO,~,~]= OptimizedQuadratureMultiREST([],W,TOL,m,DATA,[],...
                [J(COL,:); J(end,:)],[b(COL); b(end) ]);
            
        end
        
    end
    znew = cell2mat(z_glo) ;
    znew = unique(znew) ;
    disp(['Number of nonrepeated candidates =',num2str(length(znew))])
else
    znew = [] ;
end


%SNAPfNWnew = SNAPfNW(znew,:) ;
if isempty(znew)
    znew =  DATA.zINI ;
    
end
%dbstop('107')
disp('GLOBAL INTEGRATION')
disp('---------------------------------')
disp('Starting Optimized Quadrature ....')
DATA.NOSHOWITERS = 0 ;
m=DATA.npoints ;
DATA = DefaultField(DATA,'NNEGLSQ_ONLY_WHEN_NECESSARY',1) ;
[z,w,errorGLO,b,J]= OptimizedQuadratureMultiREST(SNAPfNW,W,TOL,m,DATA,znew,J,b,DATA_CUBATURE);
disp('Done')

% Include boundary points
%[z,w] =  pointsaddener(DATA,z,J,b,W,w)


TOTALTIME = toc(TIMETOTAL) ;

disp(['TOTALTIME = ',num2str(TOTALTIME)])

PLOT_ERROR = 1;

if PLOT_ERROR == 1
    %  dbstop('128')
    figure(701)
    hold on
    xlabel('Number of Iterations')
    ylabel('log(ERROR)')
    plot(log10(errorGLO/100))
    
    
    figure(702)
    hold on
    xlabel('Number of Iterations')
    ylabel('(ERROR) %')
    plot((errorGLO))
end


