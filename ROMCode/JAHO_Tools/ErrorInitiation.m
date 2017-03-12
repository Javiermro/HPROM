% see jaho_runGLOBAL_NEW




 
    
    % Initialization of sundry variables
    % -----------------------------------
    %dbstop('120')
 
    %NAME_WS_SAVE_RUN = NAME_WS_SAVE_RUN_GLO.(NNTTLL){iGLOBAL} ;
    ERROR.XX =  [];
    ERROR.YY =  [];
    ERROR.ZZ =  [];
    ERROR.XY =  [];
    
    ERROR_GLO.EUCLIDRMAX.VALUE = [];
    ERROR_GLO.EUCLIDRMAX.LEGY  = ['Max normEUC(\Sigma-\Sigma^{II})/normEUC(\Sigma), (%)'];
    ERROR_GLO.EUCLIDRMAX.TITLE  = {'Euclidean error'};
    ERROR_GLO.EUCLIDRMAX.NFIG  = 600;
    
    ERROR_GLO.HOMOG_ABS_COMP.VALUE = [];
    ERROR_GLO.HOMOG_ABS_COMP.LEGY  = ['Max  abs(\sigma^{homog}-\sigma^{homog,II}) (MPa)'];
    ERROR_GLO.HOMOG_ABS_COMP.TITLE  = {'Component: \sigma_{x}','Component: \sigma_{y}','Component:  \sigma_{z}','Component:  \tau_{xy}'};
    ERROR_GLO.HOMOG_ABS_COMP.NFIG  = [610 611 612 613];
    
    ERROR_GLO.HOMOG_ABS.VALUE = [];
    ERROR_GLO.HOMOG_ABS.LEGY  = ['Max  abs(\sigma^{homog}-\sigma^{homog,II}) (MPa)'];
    ERROR_GLO.HOMOG_ABS.TITLE  = {'Max(\sigma_{x},\sigma_{y},\sigma_{z}, \tau_{xy})'};
    ERROR_GLO.HOMOG_ABS.NFIG  = [620];
    
    % ERROR_GLO.HOMOG_REL_COMP.VALUE = [];
    % ERROR_GLO.HOMOG_REL_COMP.LEGY  = ['Max  abs((\sigma^{homog}-\sigma^{homog,II})/abs(\sigma^{homog})(%)'];
    % ERROR_GLO.HOMOG_REL_COMP.TITLE  = {'Component: \sigma_{x}','Component: \sigma_{y}','Component:  \sigma_{z}','Component:  \tau_{xy}'};
    % ERROR_GLO.HOMOG_REL_COMP.NFIG  = [630 631 632 633];
    
    % ERROR_GLO.HOMOG_REL.VALUE = [];
    % ERROR_GLO.HOMOG_REL.LEGY  = ['Max  abs((\sigma^{homog}-\sigma^{homog,II})/abs(\sigma^{homog})(%)'];
    % ERROR_GLO.HOMOG_REL.TITLE  = {'Max(\sigma_{x},\sigma_{y},\sigma_{z}, \tau_{xy})'};
    % ERROR_GLO.HOMOG_REL.NFIG  = [630];
    
    ERROR_GLO.HOMOG_RELG_COMP.VALUE = [];
    ERROR_GLO.HOMOG_RELG_COMP.LEGY  = ['Max  abs((\sigma^{homog}-\sigma^{homog,II})/',factor_adimNORM{2},'(%)'];
    ERROR_GLO.HOMOG_RELG_COMP.TITLE  = {'Component: \sigma_{x}','Component: \sigma_{y}','Component:  \sigma_{z}','Component:  \tau_{xy}'};
    ERROR_GLO.HOMOG_RELG_COMP.NFIG  = [630:634];
    
    ERROR_GLO.HOMOG_RELG.VALUE = [];
    ERROR_GLO.HOMOG_RELG.LEGY  = ['Max  abs((\sigma^{homog}-\sigma^{homog,II})/',factor_adimNORM{2},'(%)'];
    ERROR_GLO.HOMOG_RELG.TITLE  ={'Max(\sigma_{x},\sigma_{y},\sigma_{z}, \tau_{xy})'};
    ERROR_GLO.HOMOG_RELG.NFIG  = [640];
    
    
    
    CONVERGENCE_ERRORS.AT_NBASIS_STRESS = [] ;
    CONVERGENCE_ERRORS.AT_GAUSS = {} ;
    CONVERGENCE_ERRORS.FAIL_ALGORITHM = {} ;
    CONVERGENCE_ERRORS.FINAL_ALGORITHM =[] ;
    
    
    
 



% for itraj = 1:length(NAME_TEST_TRAJECT)  %% Extension to multiple trajectories
%     
%     % Initialization of sundry variables
%     % -----------------------------------
%     %dbstop('120')
%     NNTTLL = NAME_TEST_TRAJECT{itraj} ;
%     %NAME_WS_SAVE_RUN = NAME_WS_SAVE_RUN_GLO.(NNTTLL){iGLOBAL} ;
%     ERROR.(NNTTLL).XX =  [];
%     ERROR.(NNTTLL).YY =  [];
%     ERROR.(NNTTLL).ZZ =  [];
%     ERROR.(NNTTLL).XY =  [];
%     
%     ERROR_GLO.(NNTTLL).EUCLIDRMAX.VALUE = [];
%     ERROR_GLO.(NNTTLL).EUCLIDRMAX.LEGY  = ['Max normEUC(\Sigma-\Sigma^{II})/normEUC(\Sigma), (%)'];
%     ERROR_GLO.(NNTTLL).EUCLIDRMAX.TITLE  = {'Euclidean error'};
%     ERROR_GLO.(NNTTLL).EUCLIDRMAX.NFIG  = 600;
%     
%     ERROR_GLO.(NNTTLL).HOMOG_ABS_COMP.VALUE = [];
%     ERROR_GLO.(NNTTLL).HOMOG_ABS_COMP.LEGY  = ['Max  abs(\sigma^{homog}-\sigma^{homog,II}) (MPa)'];
%     ERROR_GLO.(NNTTLL).HOMOG_ABS_COMP.TITLE  = {'Component: \sigma_{x}','Component: \sigma_{y}','Component:  \sigma_{z}','Component:  \tau_{xy}'};
%     ERROR_GLO.(NNTTLL).HOMOG_ABS_COMP.NFIG  = [610 611 612 613];
%     
%     ERROR_GLO.(NNTTLL).HOMOG_ABS.VALUE = [];
%     ERROR_GLO.(NNTTLL).HOMOG_ABS.LEGY  = ['Max  abs(\sigma^{homog}-\sigma^{homog,II}) (MPa)'];
%     ERROR_GLO.(NNTTLL).HOMOG_ABS.TITLE  = {'Max(\sigma_{x},\sigma_{y},\sigma_{z}, \tau_{xy})'};
%     ERROR_GLO.(NNTTLL).HOMOG_ABS.NFIG  = [620];
%     
%     % ERROR_GLO.HOMOG_REL_COMP.VALUE = [];
%     % ERROR_GLO.HOMOG_REL_COMP.LEGY  = ['Max  abs((\sigma^{homog}-\sigma^{homog,II})/abs(\sigma^{homog})(%)'];
%     % ERROR_GLO.HOMOG_REL_COMP.TITLE  = {'Component: \sigma_{x}','Component: \sigma_{y}','Component:  \sigma_{z}','Component:  \tau_{xy}'};
%     % ERROR_GLO.HOMOG_REL_COMP.NFIG  = [630 631 632 633];
%     
%     % ERROR_GLO.HOMOG_REL.VALUE = [];
%     % ERROR_GLO.HOMOG_REL.LEGY  = ['Max  abs((\sigma^{homog}-\sigma^{homog,II})/abs(\sigma^{homog})(%)'];
%     % ERROR_GLO.HOMOG_REL.TITLE  = {'Max(\sigma_{x},\sigma_{y},\sigma_{z}, \tau_{xy})'};
%     % ERROR_GLO.HOMOG_REL.NFIG  = [630];
%     
%     ERROR_GLO.(NNTTLL).HOMOG_RELG_COMP.VALUE = [];
%     ERROR_GLO.(NNTTLL).HOMOG_RELG_COMP.LEGY  = ['Max  abs((\sigma^{homog}-\sigma^{homog,II})/',factor_adimNORM{2},'(%)'];
%     ERROR_GLO.(NNTTLL).HOMOG_RELG_COMP.TITLE  = {'Component: \sigma_{x}','Component: \sigma_{y}','Component:  \sigma_{z}','Component:  \tau_{xy}'};
%     ERROR_GLO.(NNTTLL).HOMOG_RELG_COMP.NFIG  = [630:634];
%     
%     ERROR_GLO.(NNTTLL).HOMOG_RELG.VALUE = [];
%     ERROR_GLO.(NNTTLL).HOMOG_RELG.LEGY  = ['Max  abs((\sigma^{homog}-\sigma^{homog,II})/',factor_adimNORM{2},'(%)'];
%     ERROR_GLO.(NNTTLL).HOMOG_RELG.TITLE  ={'Max(\sigma_{x},\sigma_{y},\sigma_{z}, \tau_{xy})'};
%     ERROR_GLO.(NNTTLL).HOMOG_RELG.NFIG  = [640];
%     
%     
%     
%     CONVERGENCE_ERRORS.(NNTTLL).AT_NBASIS_STRESS = [] ;
%     CONVERGENCE_ERRORS.(NNTTLL).AT_GAUSS = {} ;
%     CONVERGENCE_ERRORS.(NNTTLL).FAIL_ALGORITHM = {} ;
%     CONVERGENCE_ERRORS.(NNTTLL).FINAL_ALGORITHM =[] ;
%     
%     
%     
% end