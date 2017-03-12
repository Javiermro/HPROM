

% DEFINICION DEL OPERADOR PROYECCION
if MOD_TYPE==1 || MOD_TYPE==2
    if  ~isempty(BASIS_DATA)
        if MOD_TYPE==1
            INV_PROB   = varargin{7} ;
            TRU_POINT  = varargin{8} ;
        elseif MOD_TYPE==2
            EnergySETS = varargin{39} ; %char(varargin(39)) ; %JLM
            Set_MATRIX = varargin{40}{1} ; % JLM Set de EF del dominio continuo
            Set_BANDAS = varargin{40}{2} ; % JLM Set de EF del dominio discontinuo
            INV_PROB   = varargin{41} ;
            TRU_POINT  = varargin{42} ;            
            e_VG.nIntPoints = varargin{43} ; % npointsINPUT
        end
        
        %load(BASIS_DATA,'PHI_EPS_REG','PHI_EPS_DIS','PHI_EPS_ELAS2')
        load(BASIS_DATA,'PHI_EPS_REG','PHI_EPS_DIS')
        ModoPHI_REG = PHI_EPS_REG(:,1:nModesEPS_bar);
        ModoPHI_DIS = PHI_EPS_DIS(:,1:nModesEPS_bar2);
        %ModoPHI_ELAS2 = PHI_EPS_ELAS2(:,1:nModesEPS_bar3);
        %clear PHI_EPS_REG PHI_EPS_DIS PHI_EPS_ELAS2
        clear PHI_EPS_REG PHI_EPS_DIS
        
        % FIRST REDUCTION MODES
        e_VG.nModesEPS_bar = nModesEPS_bar;
        e_VG.nModesEPS_bar2 = nModesEPS_bar2;
        
        e_VG.nModesEPS_TOT = nModesEPS_bar+nModesEPS_bar2;
        
        if MOD_TYPE==2 % SOLO PARA SEGUNDA REDUCCION
            %    e_VG.nModesSTR_REG = NBASIS_STRESS_PREV_REG ;
            %    e_VG.nModesSTR_DIS = NBASIS_STRESS_PREV_DIS ;
            e_VG.nModesENER_REG = NBASIS_ENER_PREV_REG ;
            e_VG.nModesENER_DIS = NBASIS_ENER_PREV_DIS ;
        end
        
        %LOOP OVER THE SETS DEFINED IN THE FE TEST
        % ****************************************
        flagROM_I = 0; flagROM_II = 0;
        
        for iSet=1:e_VG.nSet
            
            if isfield(e_DatSet(iSet).e_DatMat,'e_VG') % MULTISCALE CASE
                e_DatSet(iSet).e_DatMat.e_VG.nIntPoints = varargin{43}; %JLM
                e_DatSet(iSet).e_DatMat.e_VG.nModesEPS_bar  = nModesEPS_bar;
                e_DatSet(iSet).e_DatMat.e_VG.nModesEPS_bar2 = nModesEPS_bar2;
                
                e_DatSet(iSet).e_DatMat.e_VG.nModesEPS_TOT  = nModesEPS_bar+nModesEPS_bar2;
                
                if MOD_TYPE==2 % SOLO PARA SEGUNDA REDUCCION
                    e_DatSet(iSet).e_DatMat.e_VG.nModesENER_REG = NBASIS_ENER_PREV_REG ;
                    e_DatSet(iSet).e_DatMat.e_VG.nModesENER_DIS = NBASIS_ENER_PREV_DIS ;
                end
                
                % Young modulus for consistency of internal force parameters
                e_DatSet(iSet).e_DatMat.e_VG.E_MATRIX = 10e3; %JLM
%                 e_DatSet(iSet).e_DatMat.e_VG.E_MATRIX = 1.0;
                
                % Amplification paramether for convexification procedure
                e_DatSet(iSet).e_DatMat.e_VG.K_PARAM = 10e5; %JLM
%                 e_DatSet(iSet).e_DatMat.e_VG.K_PARAM = 1.0;
                
                % SECOND REDUCTION
                if MOD_TYPE ==2
                    
                    % this part should be computed once, and stored in each multiscale set
                    if flagROM_II==0
                        
                        % HROM OFFLINE PROCEDURES
                        ROM_II = HROM_OFFLINE_PROC3(e_DatSet(iSet),ModoPHI_REG,ModoPHI_DIS,...
                            BASIS_DATA,BASIS_EXPANSION,NBASIS_ENER_PREV_REG,NBASIS_ENER_PREV_DIS,...
                            ALTERNATE_BASIS_STRESS_B,NGAUSS_GAPPY_loc,file,GREEDY_ALGOR,...
                            SNAPSHOTS_ORIGIN,CHATUR_NONREPEATED,INDICES_GIVEN_LOC,SET_GPS,...
                            PARAMETER_GREEDY,LOAD_RESULT_GREEDY,COMPARE_WITH_SNAPSHOTS,...
                            nModesEPS_bar,nModesEPS_bar2,PARTIAL_RECONSTR,WHERE_WEIGHTS,...
                            LEVEL_TRUN_INTERN,EnergySETS,Set_BANDAS,Set_MATRIX,INV_PROB,TRU_POINT,e_DatSet(iSet).e_DatMat.e_VG);
%                             LEVEL_TRUN_INTERN,EnergySETS,Set_BANDAS,Set_MATRIX,INV_PROB,TRU_POINT,e_VG);% JLM (asi era antes)                             
                        
                        flagROM_II = 1;
                    end
                    
                    % HROM_II structure for integrating the choosed Gauss Points
                    e_DatSet(iSet).e_DatMat.ROM_II = ROM_II;
                    
                    % Para impresion de valores de puntos seleccionados
                    e_DatSet(iSet).e_DatMat.e_VG.MOD_TYPE = MOD_TYPE;
                    e_DatSet(iSet).e_DatMat.e_VG.GaussFLAG = ROM_II.GaussFLAG;
                else
                    % this part should be computed once, and stored in each multiscale set
                    if flagROM_I==0
                        % HROM OFFLINE PROCEDURES - ROMI
                        [PHI_GEN,IntPhiTGen,MAT_WEAK,ValRM] = HROM_OFFLINE_PROC_ROMI(e_DatSet(iSet),ModoPHI_REG,ModoPHI_DIS,e_VG,isMICRO(1).MICRO);
                        flagROM_I = 1;
                    end
                    e_DatSet(iSet).e_DatMat.e_VG.PHI_GEN = PHI_GEN;
                    e_DatSet(iSet).e_DatMat.e_VG.IntPhiTGen = IntPhiTGen;
                    e_DatSet(iSet).e_DatMat.e_VG.MAT_WEAK = MAT_WEAK;
                    e_DatSet(iSet).e_DatMat.e_VG.ValRM = ValRM;
                    
                    e_DatSet(iSet).e_DatMat.e_VG.MOD_TYPE = MOD_TYPE;                    
                end
                
            else % MONOSCALE CASE
                
                e_VG.nModesEPS_bar  = nModesEPS_bar;
                e_VG.nModesEPS_bar2 = nModesEPS_bar2;
                e_VG.nModesEPS_TOT  = nModesEPS_bar+nModesEPS_bar2;
                
                if MOD_TYPE==2 % SOLO PARA SEGUNDA REDUCCION
                    e_VG.nModesENER_REG = NBASIS_ENER_PREV_REG ;
                    e_VG.nModesENER_DIS = NBASIS_ENER_PREV_DIS ;
                end
                
                % Young modulus for consistency of internal force parameters
                e_VG.E_MATRIX = 10e3; %JLM
%                 e_VG.E_MATRIX = 1 ;
                % Amplification paramether for convexification procedure
                e_VG.K_PARAM = 10e5; %JLM
%                 e_VG.K_PARAM = 1 ;
                
                % SECOND REDUCTION
                if MOD_TYPE ==2
                    
                    % this part should be computed once, and stored in each multiscale set
                    if flagROM_II==0
                        
                        % HROM OFFLINE PROCEDURES
                        ROM_II = HROM_OFFLINE_PROC3(e_DatSet,ModoPHI_REG,ModoPHI_DIS,...
                            BASIS_DATA,BASIS_EXPANSION,NBASIS_ENER_PREV_REG,NBASIS_ENER_PREV_DIS,...
                            ALTERNATE_BASIS_STRESS_B,NGAUSS_GAPPY_loc,file,GREEDY_ALGOR,...
                            SNAPSHOTS_ORIGIN,CHATUR_NONREPEATED,INDICES_GIVEN_LOC,SET_GPS,...
                            PARAMETER_GREEDY,LOAD_RESULT_GREEDY,COMPARE_WITH_SNAPSHOTS,...
                            nModesEPS_bar,nModesEPS_bar2,PARTIAL_RECONSTR,WHERE_WEIGHTS,...
                            LEVEL_TRUN_INTERN,EnergySETS,Set_BANDAS,Set_MATRIX,INV_PROB,TRU_POINT,e_VG);
                        
                        %clear e_DatSet
                        
                        % HROM_II structure for integrating the choosed Gauss Points
                        %e_DatSet.e_DatMat.ROM_II = ROM_II;
                        e_VG.ROM_II = ROM_II;
                        
                        % Para impresion de valores de puntos seleccionados
                        %e_DatSet.e_DatMat.e_VG.MOD_TYPE = MOD_TYPE;
                        e_VG.MOD_TYPE = MOD_TYPE;
                        %e_DatSet.e_DatMat.e_VG.GaussFLAG = ROM_II.GaussFLAG;
                        e_VG.GaussFLAG = ROM_II.GaussFLAG;
                        
                        %e_VG.ndoft  = ndoft;
                        e_VG.nElem  = ROM_II.e_DatSet_ROM.nElem ;
                        e_VG.dofpe  = ROM_II.e_DatSet_ROM.dofpe ;
                        
                        % Informacion necesaria para el proceso de calculo
                        e_VG.esImplex =e_DatSet(1).e_DatMat.esImplex;                        
                        
                        flagROM_II = 1;
                    end

                else
                    
                    % This part should be computed once, and stored in each multiscale set
                    if flagROM_I==0
                        % HROM OFFLINE PROCEDURES - ROMI
                        [PHI_GEN,IntPhiTGen,MAT_WEAK,ValRM] = HROM_OFFLINE_PROC_ROMI(e_DatSet,ModoPHI_REG,ModoPHI_DIS,e_VG,isMICRO(1).MICRO);
                        
                        e_VG.PHI_GEN = PHI_GEN;
                        e_VG.IntPhiTGen = IntPhiTGen;
                        e_VG.MAT_WEAK = MAT_WEAK;
                        e_VG.ValRM = ValRM;
                        
                        % Informacion necesaria para el proceso de calculo
                        e_VG.esImplex=e_DatSet(1).e_DatMat.esImplex;
                        
                        e_VG.MOD_TYPE = MOD_TYPE;
                        
                        flagROM_I = 1;
                    end
                    
                end
                
            end
        end
    else
        error('Snapshot basis not defined!');
    end
end

% MODEL FLAG
e_VG.MOD_TYPE = MOD_TYPE;