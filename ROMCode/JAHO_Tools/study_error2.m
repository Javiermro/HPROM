function [Hplot1LOC,legendLOC,e_VG]=study_error2(e_VG)
% Compute the "snapshot-wise" error for the ROM-II results
% Similar to study_error1

Hplot1LOC = [] ;
legendLOC = {} ;
coloresM =  ColoresMatrix(10,[],'ALEATORIO','YES');
MarkerM =  MarkerMatrix(10,[],'ALEATORIO','YES');
factor_adimNORM = e_VG.factor_adimNORM;
%dbstop('11')
ndof_sigma_SC = e_VG.ndof_sigma_SC;
TYPE_NORM = e_VG.TYPE_NORM;
PROY_SIGMA_G =  e_VG.SIGMA_ROM_II ;   % ROM-II results
SNAPSHOTS_ORIGIN = e_VG.SNAPSHOTS_ORIGIN  ;


pondfactors = e_VG.ponder_factors ;
switch e_VG.COMPARE_WITH_SNAPSHOTS
    case 'YES'
        SIGMA_SNAP_SC= e_VG.SIGMA_SNAP_SC_NW; % ROM-I results
    case 'NO_JAHO'
        SIGMA_SNAP_SC= e_VG.SIGMA_HFM; % ROM-I results
end

ERROR_TO_PRINT = [] ;
e_VG.ERROR_TO_PRINT = ERROR_TO_PRINT ;

switch TYPE_NORM
    case 'NORM2'
        normSIGMA_SNAP = sqrt(sum((SIGMA_SNAP_SC).*(SIGMA_SNAP_SC)));
        errorLOCG = sqrt(sum((SIGMA_SNAP_SC-PROY_SIGMA_G).*(SIGMA_SNAP_SC-PROY_SIGMA_G)));
    case 'NORM2_relgiv'
        errorLOCG = sqrt(sum((SIGMA_SNAP_SC-PROY_SIGMA_G).*(SIGMA_SNAP_SC-PROY_SIGMA_G)));
        
        factor_prop =  1 ;
        if isempty(factor_adimNORM)
        elseif iscell(factor_adimNORM)
            factor_loc = factor_adimNORM{1};
            if factor_loc == 1
                leg_ylabel = ['Absolute error (MPa) norm(SNAPS.-APPOX.) '];
            else
                leg_ylabel = ['Relative error (%) \norm(SNAPS.-APPOX.)/sqrt(n_g)/',factor_adimNORM{2}];
                factor_prop = 100 ;
            end
        elseif factor_adimNORM == 1
            factor_loc = 1 ;
            leg_ylabel = ['Absolute error (MPa) norm(SNAPS.-APPOX.) '];
        end
        
        normSIGMA_SNAP = sqrt(length(errorLOCG))*factor_loc*ones(size(errorLOCG));
        
        
    case 'MAX'
        normSIGMA_SNAP = max(abs(SIGMA_SNAP_SC));
        errorLOCG = max(abs(SIGMA_SNAP_SC-PROY_SIGMA_G));
    case {'L2','AVERAGE','AVERAGE_RELGIV'}  % By components
        % e_VG.ponder_factors  is a nelem*npg vector containing the
        % quadrature weight and the jacobian volume transformation
        %IND_LOC = cell(4,1);
        
        %% Volume
        VOLUME = sum(e_VG.ponder_factors);
        
        
        
        ncompt = e_VG.ntens*e_VG.nElem*e_VG.npg ;
        errorLOCG = SIGMA_SNAP_SC-PROY_SIGMA_G;
        
        %%%%%% COMPLEMENTARY INFORMATION: 'MANUEL's' NORM
        % sum(abs(SIGMA_HOMONEG.))
        
        
        
        %%%%%%
        
        
        
        for icomp = 1:e_VG.ntens
            HP = [] ;
            HL = {} ;
            % Loop over stress components
            ind_loc = icomp:e_VG.ntens:ncompt;
            
            % bsxfun (SIGMA_SNAP_SC^2.*e_VG.ponder_factors)  (to avoid "repmat")
            switch TYPE_NORM
                case 'L2'
                    
                    stress_norm = sqrt(sum(bsxfun(@times, SIGMA_SNAP_SC(ind_loc,:).^2,pondfactors)));
                    error_compG = sqrt(sum(bsxfun(@times, errorLOCG(ind_loc,:).^2,pondfactors)));
                    
                    leg_ylabel = ['Relative error (%) normL2(SNAPS.-APPOX.)/normL2(SNAPS.)'];
                    factor_prop =  100 ;
                    
                    
                    
                case {'AVERAGE','AVERAGE_RELGIV'}
                    
                    % dbstop('63')
                    
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
                    
                    
                    error_compG = abs(sum(bsxfun(@times, errorLOCG(ind_loc,:),pondfactors)));
                    
                    
                    
            end
            
            %%%%%% COMPLEMENTARY INFORMATION: 'MANUEL's' NORM
            % sum(abs(SIGMA_HOMONEG.))
            % sumABS_HF = sum(abs(sum(bsxfun(@times, SIGMA_SNAP_SC(ind_loc,:),pondfactors))));
            sumABS_HF = sum(abs(sum(bsxfun(@times, SIGMA_SNAP_SC(ind_loc,:),pondfactors))));
            sumABS_ERROR =   sum(error_compG)   ;
            
            error_MANUEL = sumABS_ERROR/sumABS_HF*100 ;
            ERROR_TO_PRINT(end+1) = error_MANUEL ;
            
            error_relG = error_compG./stress_norm ;
            tolloc = 1e-4 ;
            [aaa] =  find(stress_norm<tolloc) ;
            if ~isempty(aaa)
                error_relG(aaa) = error_compG(aaa);
            end
            
            % Plots
            figure(e_VG.NFIG_ERROR_1+icomp)
            hold on
            xlabel('Column of the stress snapshot matrix')
            ylabel(leg_ylabel)
            title(['Error analysis for stress component: ',e_VG.STRESS_LABELS{icomp}])
            
            
            HP(end+1) = plot(error_relG*factor_prop,'Color',coloresM(e_VG.ntens+icomp,:),'Marker',MarkerM{e_VG.ntens+icomp});
            
            HL{end+1} = ['ROM-II,(',strrep(e_VG.GREEDY_ALGOR,'_',' '),'), ','N_S','=',num2str(e_VG.NBASIS_STRESS_PREV),'+ N_U',';','N_U=',...
                num2str(e_VG.nmodesU),';N_g= ',num2str(length(e_VG.npoints_sigma)),';',noguion(TYPE_NORM),', ',noguion(e_VG.WHERE_WEIGHTS),...
                '; e_m \approx ',num2str(error_MANUEL,2),' %','; OSNAP =',SNAPSHOTS_ORIGIN];
            
            legend(HP,HL);
            
            
        end
        
        e_VG.ERROR_TO_PRINT = ERROR_TO_PRINT ;
        
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
        
        tolloc = 1e-8 ;
        [aaa] =  find(normSIGMA_SNAP<tolloc) ;
        %         if ~isempty(aaa)
        %             errorLOC_rel(aaa) = errorLOC(aaa);
        %         end
        %
        %         Hplot1LOC(end+1) = plot(errorLOC_rel*100,'Color',coloresM(1,:),'Marker',MarkerM{1});
        %         legendLOC{end+1} = ['Best approximation, ','N_S','=',num2str(e_VG.NBASIS_STRESS_PREV),'+ N_U',';','N_U=',...
        %             num2str(e_VG.nmodesU),', ',noguion(TYPE_NORM),', ',noguion(e_VG.WHERE_WEIGHTS)];
        
        errorLOC_relG = errorLOCG./normSIGMA_SNAP ;
        
        
        if ~isempty(aaa)
            errorLOC_relG(aaa) = errorLOCG(aaa);
        end
        
        Hplot1LOC(end+1) = plot(errorLOC_relG*100,'Color',coloresM(2,:),'Marker',MarkerM{2});
        
        legendLOC{end+1} = ['ROM-II,(',strrep(e_VG.GREEDY_ALGOR,'_',' '),'), ','N_S','=',num2str(e_VG.NBASIS_STRESS_PREV),'+ N_U',';','N_U=',...
            num2str(e_VG.nmodesU),';N_g= ',num2str(length(e_VG.npoints_sigma)),';',noguion(TYPE_NORM),', ',noguion(e_VG.WHERE_WEIGHTS),...
            '; OSNAP =',SNAPSHOTS_ORIGIN];
        
        legend(Hplot1LOC,legendLOC);
    case 'NORM2_relgiv'
        
        %dbstop('222')
        figure(e_VG.NFIG_ERROR_1)
        hold on
        xlabel('Column of the stress snapshot matrix')
        ylabel(leg_ylabel)
        
        errorLOC_relG = errorLOCG./factor_adimNORM{1} ;
        
        
        
        Hplot1LOC(end+1) = plot(errorLOC_relG*100,'Color',coloresM(2,:),'Marker',MarkerM{2});
        
        legendLOC{end+1} = ['ROM-II,(',strrep(e_VG.GREEDY_ALGOR,'_',' '),'), ','N_S','=',num2str(e_VG.NBASIS_STRESS_PREV),'+ N_U',';','N_U=',...
            num2str(e_VG.nmodesU),';N_g= ',num2str(length(e_VG.npoints_sigma)),';',noguion(TYPE_NORM),', ',noguion(e_VG.WHERE_WEIGHTS)];
        
        legend(Hplot1LOC,legendLOC);
        
        
        
end

%%%%%%%%
% NORMS TO BE PLOTTED OUTSIDE (REGARDLESS TYPE_NORM variable)
% -----------------------------------------------------------

% EUCLIDEAN NORM (RELATIVE), maximum
% --------------
normSIGMA_SNAP = sqrt(sum((SIGMA_SNAP_SC).*(SIGMA_SNAP_SC)));
errorLOCG = sqrt(sum((SIGMA_SNAP_SC-PROY_SIGMA_G).*(SIGMA_SNAP_SC-PROY_SIGMA_G)));
tolloc = 1e-8 ;
[aaa] =  find(normSIGMA_SNAP<tolloc) ;
errorLOC_relG = errorLOCG./normSIGMA_SNAP*100 ;
if ~isempty(aaa)
    errorLOC_relG(aaa) = errorLOCG(aaa)*100;
end
ERRORS.EUCLIDRMAX  = max(errorLOC_relG) ;
%%%%%%%%
% MAXIMUM DIFFERENCE BETWEEN AVERAGED STRESSES



ncompt = e_VG.ntens*e_VG.nElem*e_VG.npg ;
errorLOCG = SIGMA_SNAP_SC-PROY_SIGMA_G;
HOMOG_ABS_COMP = [] ;
HOMOG_RELG_COMP = [] ;
%HOMOG_REL_COMP = [] ;
for icomp = 1:e_VG.ntens
    ind_loc = icomp:e_VG.ntens:ncompt;
    % Absolute error
    % --------------------------------------------------------------------
    % -
    error_compG = abs(sum(bsxfun(@times, errorLOCG(ind_loc,:),pondfactors)));
    HOMOG_ABS_COMP(end+1) = max(error_compG) ;
    %     % Error relative to norm STRESS
    %     % --------------------------------------------------------------------
    %     dbstop('273')
    %     stress_normREL = abs(sum(bsxfun(@times, SIGMA_SNAP_SC(ind_loc,:),pondfactors)));
    %     errorLOC_relG = error_compG./stress_normREL ;
    %     tolloc = 1e-8 ;
    %     [aaa] =  find(stress_normREL<tolloc) ;
    %     if ~isempty(aaa)
    %         errorLOC_relG(aaa) = [];
    %     end
    %     HOMOG_REL_COMP(end+1) = max(errorLOC_relG) ;
    % Error relative to a GIVEN PARAMETER
    VOLUME = sum(e_VG.ponder_factors);
    stress_normREL = factor_adimNORM{1}*VOLUME*ones(1,size(SIGMA_SNAP_SC,2))    ;
    errorLOC_relG = error_compG./stress_normREL*100 ;
    HOMOG_RELG_COMP(end+1) = max(errorLOC_relG) ;
end
HOMOG_ABS = max(HOMOG_ABS_COMP) ;

HOMOG_RELG =  max(HOMOG_RELG_COMP) ;


ERRORS.HOMOG_ABS_COMP = HOMOG_ABS_COMP ;
ERRORS.HOMOG_ABS = HOMOG_ABS ;
%ERRORS.HOMOG_REL_COMP =  HOMOG_REL_COMP ;
%ERRORS.HOMOG_REL =  HOMOG_REL ;
ERRORS.HOMOG_RELG_COMP =  HOMOG_RELG_COMP ;
ERRORS.HOMOG_RELG =  HOMOG_RELG ;

e_VG.ERRORS = ERRORS ;
