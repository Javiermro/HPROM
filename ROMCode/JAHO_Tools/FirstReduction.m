
% ##########################
% FIRST REDUCTION
% REDUCTION IN DISPLACEMENTS
% ##########################

DEFAULT_U = 1 ;  % Using nmodesU
if exist([path_file,'/DATAOUT/'])~=7; mkdir([path_file,'/DATAOUT/']); end
BASENN = [path_file,'/DATAOUT/','ROMI_1',IDENT] ;
% URINAME_WS_SAVE_RUN_GLO = [path_file,'/DATAOUT/','ROMI_1',IDENT,'.mat']; % For a single run ...
NOLOAD = 1 ; % By default, no workspace variable is loaded
try % It may happen that USE_NMODES_U is not defined
    if USE_NMODES_U == 1
        DEFAULT_U = 0 ; % USing nmodesU_GLO
    else
        nmodesU_GLO ={nmodesU} ;
        ROM_I_LOAD_WS_RUN_GLO{1} = 0;
    end
catch
    nmodesU_GLO ={nmodesU} ;
end

if ~exist('URINAME_WS_SAVE_RUN_GLO')  % To avoid conflict with previous versions ...
    URINAME_WS_SAVE_RUN_GLO = {} ;
end

for iproj = 1:length(nModesEPS_bar_GLO)
    nModesEPS_bar = nModesEPS_bar_GLO{iproj};
    nModesEPS_bar2 = nModesEPS_bar2_GLO{iproj};
    %nModesEPS_bar3 = nModesEPS_bar3_GLO{iproj};
    try
        nnn = URINAME_WS_SAVE_RUN_GLO{iproj} ;
    catch
        URINAME_WS_SAVE_RUN_GLO{iproj} = [BASENN,'_',num2str(nmodesU),'.mat'];
    end
      
    if ROM_I_LOAD_WS_RUN_GLO{iproj} == 0 || exist(URINAME_WS_SAVE_RUN_GLO{iproj})~=2
    nLabImp = 1; MOD_TYPE = 1;
         [tROMI,OTHEROUT] = analysis(file,path_file,nLabImp,MOD_TYPE,'SEL_SNPSH_STRESS',SEL_SNPSH_STRESS,'tol_disp',...
            tol_disp,'nModesEPS_bar',nModesEPS_bar,'nModesEPS_bar2',nModesEPS_bar2,...
            'PLOT_EIGSPEC',PLOT_EIGSPEC,'BASIS_DATA',BASIS_DATA);
        % As usual, DATAOUT is stored in an identifiable workspace so that
        % information can be retrieved without running again analysis_ROM
        DATAOUT.tROMI = tROMI;
        DATAOUT.OTHEROUT = OTHEROUT;
        
        save(URINAME_WS_SAVE_RUN_GLO{iproj},'DATAOUT');
        save([path_file,'tiempo_TOTAL_ROMI.mat'],'tROMI');
    else
        %  dbstop('83')
        load(URINAME_WS_SAVE_RUN_GLO{iproj},'DATAOUT');
        load([path_file,'tiempo_TOTAL_ROMI.mat'],'tROMI');
    end
end