function parsave_INT_VAR(NAME_WS_STORE,SingleV_intern_GLO,PHI_intern_GLO,INV_VAR_PLOT_GLO,DATA_VARPLOT)

% parsave_SVD_U(NAME_WS_STORE,SNAP_DISP{igroup},SingleV_DISP,PHI_DISP,TRAJECTORIES,list_snapshots_show,NAME_SCRIPT);
% parsave_SVD_STRESS(NAME_WS_STORE,SNAP_STRESSES{igroup},SingleV_STRESS,PHI_STRESS,V_STRESS)

% save(LABEL_STORE(igroup).NAME_WS_STORE,'SingleV_intern_GLO','PHI_intern_GLO','INV_VAR_PLOT_GLO','DATA_VARPLOT','-append') ;
% save(NAME_WS_STORE,'SNAP_STRESSES','SingleV_STRESS','PHI_STRESS','V_STRESS','-append') ;

save(NAME_WS_STORE,'SingleV_intern_GLO','PHI_intern_GLO','INV_VAR_PLOT_GLO','DATA_VARPLOT','-append') ;

end