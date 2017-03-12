function [Hplot1LOC,legendLOC]=study_error1(PHI,SIGMA_SNAP_SC,ndof_sigma_SC,e_VG,TYPE_NORM)
% Compute the "snapshot-wise" error
% First of all, we plot the error contained in the "best
% approximation" (the projection)
% % DEFAULT INPUTS
% REDUCED_MODEL_II = 'NO' ; %
% %% EXTRACTING INPUTS
% %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% varginOR = varargin ;
% FdnamesInputs = {'REDUCED_MODEL_II'};
% AuxDATA=WriteAuxFdNamesNEW(FdnamesInputs,varginOR);
% for id = 1:length(AuxDATA);
%     eval(AuxDATA{id});
% end
% %E_JAHO




Hplot1LOC = [] ;
legendLOC = {} ;
coloresM =  ColoresMatrix(10,[],'ALEATORIO','YES');
MarkerM =  MarkerMatrix(10,[],'ALEATORIO','YES');
 SNAPSHOTS_ORIGIN = e_VG.SNAPSHOTS_ORIGIN  ;


% PROJECITON
%dbstop('29')
PROY_SIGMA = PHI*PHI'*SIGMA_SNAP_SC ;
%% GAPPY RECONSTRUCTION
PHI_G = PHI(ndof_sigma_SC,:);
M = PHI_G'*PHI_G ;
PROY_SIGMA_G = PHI*(M\PHI_G')*SIGMA_SNAP_SC(ndof_sigma_SC,:);

% switch e_VG.WHERE_WEIGHTS
%             case 'WinB'
%
%             case 'WinS'
%                 pondfactors = ones(size(e_VG.ponder_factors));
%             otherwise
%                 error('Not implemented ...')
% end

pondfactors = e_VG.ponder_factors ;

% If e_VG.WHERE_WEIGHTS == WinS, we must divide by the weighting factors
% !!!  (both PROY_SIGMA and PROY_SIGMA_G)


switch e_VG.WHERE_WEIGHTS
    case 'WinB'
        %              % Nothing is done
    case 'WinS'
        % Dividing by pondfactors
        %pondfactors
     %   dbstop('56')
        ppff_tot=  repmat(pondfactors,[1,e_VG.ntens]);
        ppff_tot = reshape(ppff_tot',[e_VG.ntens*e_VG.npg*e_VG.nElem,1]);
        PROY_SIGMA =  bsxfun(@rdivide,PROY_SIGMA,ppff_tot);
        PROY_SIGMA_G=  bsxfun(@rdivide,PROY_SIGMA_G,ppff_tot);
        SIGMA_SNAP_SC=  bsxfun(@rdivide,SIGMA_SNAP_SC,ppff_tot);
    otherwise
        error('Not implemented ...')
end


switch TYPE_NORM
    case 'NORM2'
        normSIGMA_SNAP = sqrt(sum((SIGMA_SNAP_SC).*(SIGMA_SNAP_SC)));
        errorLOC = sqrt(sum((SIGMA_SNAP_SC-PROY_SIGMA).*(SIGMA_SNAP_SC-PROY_SIGMA)));
        errorLOCG = sqrt(sum((SIGMA_SNAP_SC-PROY_SIGMA_G).*(SIGMA_SNAP_SC-PROY_SIGMA_G)));
    case 'MAX'
        normSIGMA_SNAP = max(abs(SIGMA_SNAP_SC));
        errorLOC = max(abs(SIGMA_SNAP_SC-PROY_SIGMA));
        errorLOCG = max(abs(SIGMA_SNAP_SC-PROY_SIGMA_G));
    case {'L2','AVERAGE','AVERAGE_RELGIV'}  % By components
        % e_VG.ponder_factors  is a nelem*npg vector containing the
        % quadrature weight and the jacobian volume transformation
        %IND_LOC = cell(4,1);
        
        %% Volume
        VOLUME = sum(e_VG.ponder_factors );
        
        
        
        ncompt = e_VG.ntens*e_VG.nElem*e_VG.npg ;
        errorLOC = SIGMA_SNAP_SC-PROY_SIGMA;
        errorLOCG = SIGMA_SNAP_SC-PROY_SIGMA_G;
        
        for icomp = 1:e_VG.ntens
            HP = [] ;
            HL = {} ;
            % Loop over stress components
            ind_loc = icomp:e_VG.ntens:ncompt;
            
            % bsxfun (SIGMA_SNAP_SC^2.*e_VG.ponder_factors)  (to avoid "repmat")
            switch TYPE_NORM
                case 'L2'
                    
                    stress_norm = sqrt(sum(bsxfun(@times, SIGMA_SNAP_SC(ind_loc,:).^2,pondfactors)));
                    error_comp = sqrt(sum(bsxfun(@times, errorLOC(ind_loc,:).^2,pondfactors)));
                    error_compG = sqrt(sum(bsxfun(@times, errorLOCG(ind_loc,:).^2,pondfactors)));
                    
                    leg_ylabel = ['Relative error (%) normL2(SNAPS.-APPOX.)/normL2(SNAPS.)'];
                    factor_prop =  100 ;
                    
                    
                    
                case {'AVERAGE','AVERAGE_RELGIV'}
                    
                    % dbstop('63')
                    factor_adimNORM = e_VG.factor_adimNORM;
                    switch TYPE_NORM
                        case 'AVERAGE'
                            stress_norm = abs(sum(bsxfun(@times, SIGMA_SNAP_SC(ind_loc,:),pondfactors)));
                            leg_ylabel = ['Relative error (%) \int(SNAPS.-APPOX.)/\int(SNAPS.)'];
                            factor_prop =  100 ;
                        case 'AVERAGE_RELGIV'
                            factor_prop =  1 ;
                            if isempty(factor_adimNORM)
                            elseif iscell(factor_adimNORM)
                                factor_loc = factor_adimNORM{1};
                                if factor_loc == 1
                                    leg_ylabel = ['Absolute error (MPa) \int(SNAPS.-APPOX.)dV/V '];
                                else
                                    leg_ylabel = ['Relative error (%) \int(SNAPS.-APPOX.)dV/V/',factor_adimNORM{2}];
                                    factor_prop = 100 ;
                                end
                            elseif factor_adimNORM == 1
                                factor_loc = 1 ;
                                leg_ylabel = ['Absolute error (MPa) \int(SNAPS.-APPOX.) '];
                            end
                            stress_norm = factor_loc*VOLUME*ones(1,size(SIGMA_SNAP_SC,2))    ;
                            
                            
                    end
                    
                    error_comp = abs(sum(bsxfun(@times, errorLOC(ind_loc,:),pondfactors)));
                    error_compG = abs(sum(bsxfun(@times, errorLOCG(ind_loc,:),pondfactors)));
                    
            end
            
            error_rel = error_comp./stress_norm ;
            error_relG = error_compG./stress_norm ;
            tolloc = 1e-4 ;
            [aaa] =  find(stress_norm<tolloc) ;
            if ~isempty(aaa)
                error_rel(aaa)  = error_comp(aaa);
                error_relG(aaa) = error_compG(aaa);
            end
            
            % Plots
            figure(e_VG.NFIG_ERROR_1+icomp)
            hold on
            xlabel('Column of the stress snapshot matrix')
            ylabel(leg_ylabel)
            title(['Error analysis for stress component: ',e_VG.STRESS_LABELS{icomp}])
            
            HP(end+1) = plot(error_rel*factor_prop,'Color',coloresM(icomp,:),'Marker',MarkerM{icomp});
            HL{end+1} = ['Best approximation, ','N_S','=',num2str(e_VG.NBASIS_STRESS_PREV),'+ N_U',';','N_U=',...
                num2str(e_VG.nmodesU),', ',noguion(TYPE_NORM),', ',noguion(e_VG.WHERE_WEIGHTS)];
            
            HP(end+1) = plot(error_relG*factor_prop,'Color',coloresM(e_VG.ntens+icomp,:),'Marker',MarkerM{e_VG.ntens+icomp});
            
            HL{end+1} = [strrep(e_VG.GREEDY_ALGOR,'_',' '),', ','N_S','=',num2str(e_VG.NBASIS_STRESS_PREV),'+ N_U',';','N_U=',...
                num2str(e_VG.nmodesU),';N_g= ',num2str(length(e_VG.npoints_sigma)),';',noguion(TYPE_NORM),', ',noguion(e_VG.WHERE_WEIGHTS),'; OSNAP =',SNAPSHOTS_ORIGIN];
            
            legend(HP,HL);
            
            
        end
        
        %
        %         dbstop('34')
        %         pff = repmat(e_VG.ponder_factors,[1,e_VG.ntens]);
        %         pff = reshape(pff',[e_VG.ntens*e_VG.nElem*e_VG.npg,1]);
        %         errorLOC = SIGMA_SNAP_SC-PROY_SIGMA;
        %         errorLOCG = SIGMA_SNAP_SC-PROY_SIGMA_G;
        %         normSIGMA_SNAP = SIGMA_SNAP_SC;
        % L2-norm of SIGMA_SNAP_SC
        
        
    otherwise
        error('Error analysis not implemented')
        
end





switch TYPE_NORM
    case {'MAX','NORM2'}
        figure(e_VG.NFIG_ERROR_1)
        hold on
        xlabel('Column of the stress snapshot matrix')
        ylabel('Relative error (%) norm(SNAPS.-APPOX.)/norm(SNAPS.)')
        errorLOC_rel = errorLOC./normSIGMA_SNAP ;
        tolloc = 1e-8 ;
        [aaa] =  find(normSIGMA_SNAP<tolloc) ;
        if ~isempty(aaa)
            errorLOC_rel(aaa) = errorLOC(aaa);
        end
        
        Hplot1LOC(end+1) = plot(errorLOC_rel*100,'Color',coloresM(1,:),'Marker',MarkerM{1});
        legendLOC{end+1} = ['Best approximation, ','N_S','=',num2str(e_VG.NBASIS_STRESS_PREV),'+ N_U',';','N_U=',...
            num2str(e_VG.nmodesU),', ',noguion(TYPE_NORM),', ',noguion(e_VG.WHERE_WEIGHTS),'; OSNAP =',SNAPSHOTS_ORIGIN];
        
        errorLOC_relG = errorLOCG./normSIGMA_SNAP ;
        
        
        if ~isempty(aaa)
            errorLOC_relG(aaa) = errorLOCG(aaa);
        end
        
        Hplot1LOC(end+1) = plot(errorLOC_relG*100,'Color',coloresM(2,:),'Marker',MarkerM{2});
        
        legendLOC{end+1} = [strrep(e_VG.GREEDY_ALGOR,'_',' '),', ','N_S','=',num2str(e_VG.NBASIS_STRESS_PREV),'+ N_U',';','N_U=',...
            num2str(e_VG.nmodesU),';N_g= ',num2str(length(e_VG.npoints_sigma)),';',noguion(TYPE_NORM),', ',noguion(e_VG.WHERE_WEIGHTS),'; OSNAP =',SNAPSHOTS_ORIGIN];
        
        legend(Hplot1LOC,legendLOC);
end

