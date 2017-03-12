function [npoints_gappy PHI_AMP STUDY_GREEDY_LOC ndof_sigma_SC npoints_sigma ...
    numat_B] =DetermineNpoingGAUSS(e_VG,PHI,U_B,nGTOT,U,NGAUSS_GAPPY_loc)
    

% See f_MatBT_ROM_III.m

ndof_sigma_SC =[] ;
STUDY_GREEDY_LOC = 1 ;
npoints_sigma = [] ;
numat_B = [] ;
if isempty(NGAUSS_GAPPY_loc)
    npoints_gappy = size(PHI,2);
    %npoints_gappy = round(size(PHI,2)/4);
    PHI_AMP =PHI ;
else
    
    if isstr(NGAUSS_GAPPY_loc)
        if strcmp(NGAUSS_GAPPY_loc,'ALL')
            STUDY_GREEDY_LOC  = 0 ;
            ndof_sigma_SC = 1:e_VG.ntens*e_VG.npg*e_VG.nElem;
            npoints_gappy =   e_VG.npg*e_VG.nElem;
        end
        
    else
        
        
        if isnumeric(e_VG.NGAUSS_GAPPY_loc)
            NGAUSS_GAPPY_loc = e_VG.NGAUSS_GAPPY_loc ;
            
            if NGAUSS_GAPPY_loc == nGTOT
                % Full vector; no need of greedy (heuristic) selection
                %  dbstop('341')
                STUDY_GREEDY_LOC  = 0 ;
                ndof_sigma_SC = 1:e_VG.ntens*e_VG.npg*e_VG.nElem;
                npoints_gappy =   e_VG.npg*e_VG.nElem;
                
                % Definition of array "npoints_sigma"
                % = [nelement(1)  ngausspoint(1)
                %     ....          .......]
                % npoints_sigma = zeros()
                FUND_BLOQ = repmat([1:e_VG.nElem]',[1 e_VG.npg])';
                FUND_BLOQ = reshape(FUND_BLOQ,e_VG.npg*e_VG.nElem,1);
                NGAUSSTILE = repmat([1:e_VG.npg]',[e_VG.nElem,1]);
                npoints_sigma = [FUND_BLOQ NGAUSSTILE];
                
                %%%% Global index of chosen gauss points
                numat_B = 1:e_VG.npg*e_VG.nElem;
            end
            
        elseif iscell(e_VG.NGAUSS_GAPPY_loc)
            if strcmp(e_VG.NGAUSS_GAPPY_loc{1},'VARIABLE')
                nexp = e_VG.NGAUSS_GAPPY_loc{2};
                nU = size(U_B,2);
                nS = imod_stress ;
                rS = rankS;
                nGTOT = e_VG.npg*e_VG.nElem;
                cocLOC = (nGTOT-nU)/(rS-nU).^(nexp);
                NGAUSS_GAPPY_loc =  ceil((nS.^nexp)*cocLOC) ;
                
                
                %                  nGTOT = npg*nGAUSS;
                %                 cocLOC = (nGTOT-nU)/(rS-nU).^(nexp);
                %                 NGAUSS_GAPPY_loc =  ceil((nS.^nexp)*cocLOC) ;
                %                 NGAUSS_GAPPY_loc(NGAUSS_GAPPY_loc<nS) =  nS(NGAUSS_GAPPY_loc<nS);
                %                 NGAUSS_GAPPY_loc(NGAUSS_GAPPY_loc>nGTOT-nU) = nGTOT-nU;
                
                if NGAUSS_GAPPY_loc<nS
                    NGAUSS_GAPPY_loc=nS;
                elseif NGAUSS_GAPPY_loc>=nGTOT-nU
                    NGAUSS_GAPPY_loc = nGTOT-nU ;
                    STUDY_GREEDY_LOC  = 0 ;
                    ndof_sigma_SC = 1:e_VG.ntens*e_VG.npg*e_VG.nElem;
                    % npoints_gappy =   e_VG.npg*e_VG.nElem;
                end
            end
            
        end
        
        
        npoints_gappy = min(NGAUSS_GAPPY_loc + size(U_B,2),nGTOT) ;
        % PHI_AMP =[U(:,1:NGAUSS_GAPPY_loc)   U_B ];
        %dbstop('291')
        % elegir_aleatorios = 0 ;
        
        if  STUDY_GREEDY_LOC == 1
            STRATEGIA = 1 ; % IT SEEMS TO WORK BETTER
            if STRATEGIA == 1
                
                nampliac = min(NGAUSS_GAPPY_loc,size(U,2));
                PHI_AMP =[U(:,1:nampliac)   U_B ];
                
                % What happens if NGAUSS_GAPPY_loc>size(U,2)
                if NGAUSS_GAPPY_loc>size(U,2)
                    puntos_adicionales =NGAUSS_GAPPY_loc - size(U,2) ;
                end
            elseif STRATEGIA == 2
                mod_ad  = npoints_gappy-size(PHI,2);
                
                if mod_ad>0
                    % PHI_AMP =[PHI U(:,imod_stress+1:imod_stress+mod_ad)];
                    col_adic = ceil(mod_ad/size(PHI,2));
                    if col_adic>1
                        PHI_LOC_BOR = repmat(PHI,1,col_adic);
                    else
                        PHI_LOC_BOR = PHI ;
                    end
                    PHI_AMP =[PHI PHI_LOC_BOR(:,1:mod_ad)];
                end
            end
        end
        
        
        
        
    end
    
    
end