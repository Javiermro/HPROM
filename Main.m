
clear

% incluye todos los archivos .m de los subdirectorios /ROM...
addpath(genpath('/home/javiermro/Projects/HPROM'));
% exeptua los archivos .m del directorio /MEF-NL-MultiEscala_Totales_ROM/
% rmpath('C:/CIMNE/ROM_IntrinsicDD/ROMCode/MultiEscala/MEF-NL-MultiEscala_Totales_ROM/')
%Name of trajectory - user defined
NNTTLL = 'Generic_05_005_02';

path_file= '/home/javiermro/Projects/HPROM/Examples/Mi_GD_J2_StructHole/'; file = 'RVE_StructHole10.mfl';
% path_file= 'C:/CIMNE/Codes/ROM/SIMU/ROM/MS_Uni_02_RigidHole05/' ;file = 'Macro2.mfl';
% path_file= 'C:/CIMNE/Codes/ROM/SIMU/ROM/MS_GD_J2_Prob01/' ;file = 'Macro.mfl';
% path_file= 'C:/CIMNE/Codes/ROM/SIMU/ROM/Mi_GD_J2_TEST13/' ;file = 'RVE5X5_Periodico.mfl'; 
% path_file= 'C:/CIMNE/Codes/ROM/SIMU/ROM/Mi_GD_J2_Lam01/' ;
% path_file= 'C:/CIMNE/Codes/ROM/SIMU/ROM/Mi_GD_J2_RigidHole05/' ; file = 'RVE_Hole05.mfl';
% path_file= 'C:/CIMNE/Codes/ROM/SIMU/ROM/Mi_GD_J2_RigidHole04/' ; file = 'RVE_Hole04.mfl';
% path_file= 'C:/CIMNE/Codes/ROM/SIMU/ROM/Mi_GD_J2_RigidHole03/' ; file = 'RVE_Hole03.mfl';
% path_file= 'C:/CIMNE/Codes/ROM/SIMU/ROM/Mi_GD_J2_Hole03/' ; file = 'RVE_Hole03.mfl';
% path_file= '/home/javiermro/Projects/HPROM/SIMU/ROM/Mi_GD_J2_Hole05/' ; file = 'RVE_Hole05.mfl';
% path_file= 'C:/CIMNE/ROM_IntrinsicDD/SIMU/ROM/Mi_GD_J2_Struct_Hole01/' ; file = 'RVE_Struct_Hole01.mfl';
% path_file= 'C:/CIMNE/Codes/ROM/SIMU/ROM/Mi_GD_J2_Hole/' ;
% path_file= 'C:/CIMNE/Codes/ROM/SIMU/ROM/Mi_GD_J2_RigidHole/' ;
% path_file= 'C:/CIMNE/Codes/ROM/SIMU/ROM/Mi_GD_J2_TEST13/' ;
% file = 'RVE_Hole_Periodico.mfl' ; %'RVE2500_Periodico.mfl'; %'RVE11_Periodico.mfl'

% isMICRO(1).MICRO =0; % For macro & multiscale models
isMICRO(1).MICRO =1; FACT = 1;% For RVE analysis
isMICRO(1).epsilon_Macro0=FACT*[0.2; 0; 0; 0; 0]; 

%Strain and Energy Basis
% BASIS_DATA=[path_file,'/BASIS/allStrainEnergyModes4ROMI_TEST13_InSnap75_GDJ2.mat'];
%% RigidHole05                
% BASIS_DATA=[path_file,'/BASIS/allStrainEnergyModes4ROMI_RigidHOLE05_InSnap89_GDJ2_BIS3.mat'];
%% RigidHole04
% BASIS_DATA=[path_file,'/BASIS/allStrainEnergyModes4ROMI_RigidHOLE04_InSnap85_GDJ2_BIS3.mat']; %'];   %  archivo mat que contiena las bases de def. y energ-. JLM
%% RigidHole03
% BASIS_DATA=[path_file,'/BASIS/allStrainEnergyModes4ROMI_RigidHOLE03_InSnap85_GDJ2_BIS3.mat']; %'];   %  archivo mat que contiena las bases de def. y energ-. JLM
%% Hole03
% BASIS_DATA=[path_file,'/BASIS/allStrainEnergyModes4ROMI_HOLE03_nPos74_GDJ2.mat']; %'];   %  archivo mat que contiena las bases de def. y energ-. JLM
%% Hole05
% BASIS_DATA=[path_file,'/BASIS/allStrainEnergyModes4ROMI_HOLE05_nPos74_GDJ2.mat']; %'];   %  archivo mat que contiena las bases de def. y energ-. JLM
% BASIS_DATA=[path_file,'/BASIS/allStrainEnergyModes4ROMI_HOLE05_nPos74_GDJ2_PSI_E_.mat']; %'];   %  archivo mat que contiena las bases de def. y energ-. JLM
% BASIS_DATA=[path_file,'/BASIS/allStrainEnergyModes4ROMI_HOLE05_nPos74_GDJ2_PSI_EP.mat']; %'];   %  archivo mat que contiena las bases de def. y energ-. JLM
%% RVE_StructHole10
BASIS_DATA=[path_file,'/BASIS/allStrainEnergyModes_RVE_StructHole10_nPos83_Phi_ela.mat']; 


%% MODES %%%%%% ..._MODES = [REG_Domain(*)  DIS_Domain(*) ] (*) primero van los modos elasticos y los que siguen son los plasticos 
STRAIN_MODES = [4 20];%[3 23] ; %[137 172]  TOTAL: [434 1146]
ENERGY_MODES = [10 50]; %[106 149] TOTAL: [154 524]
%%
LEVELS_OF_TRUN_PHI_VAR  = STRAIN_MODES(1);  %6  5;       % 5 13;  % STRAINS -----> REGULAR DOMAIN SET=2 
LEVELS_OF_TRUN_PHI_VAR2 = STRAIN_MODES(2); %50; %50 5;       % 5 13;  % STRAINS -----> DISCONTINUOUS DOMAIN SET=1

LEVELS_OF_TRUN_ENER_REG = ENERGY_MODES(1); %10; %10 20;      % 16;  % 16 STRESS -----> REGULAR DOMAIN SET=2 
LEVELS_OF_TRUN_ENER_DIS = ENERGY_MODES(2); %12;%250;%280 20;      %9; %  12 STRESS -----> DISCONTINUOUS DOMAIN SET=1 

% file for save miscelaneous data 
% NAME_WE_TRAJE = [path_file,'/DATAOUT/TRAJECTORY_Generic_05_005_02.mat'] ;
NAME_WE_TRAJE = [path_file,'/DATAOUT/TRAJECTORY_Generic_01.mat'] ;

%name of files in which you want to save the SVD of matrix U,B and stresses
SVD_B_ws      = 'DATA_JAHO/SVD_B_OneElement_ROMHF.mat';
SVD_U_ws      = 'DATA_JAHO/SVD_U_OneElement_ROM_HF.mat' ;
SVD_STRESS_ws = 'DATA_JAHO/SVD_STRESS_OneElement_ROMHF.mat' ;

GREEDY_ALGOR_TRY = 'CHATUR' ; % CHATUR % HIERARCHICAL_cycle
CHATUR_NONREPEATED = 0;
LOAD_RESULT_GREEDY = 0;

WHERE_WEIGHTS = 'WinB';
SNAPSHOTS_ORIGIN = 'HF' ;
ALTERNATE_BASIS_STRESS_B = 0 ;
BASIS_EXPANSION = 'NO_NORMALIZATION';
ECONOMIC_STORAGE = 0;
SET_GPS = {[],[]} ; %ptos de gauss que son elegidos explicitamente {[Pto.Reg.][Pto.Dis]}
Set_BANDAS = {[1],[2]}; %Sets de EF del HF del dominio regular y discontinuo, puedes ser mas de uno, ej. {[1][2 3]}, {[REG][DIS]}
EnergySETS = 'DECOMP_01'; % Descomposicion de dominio: NO: NO_DECOMP; SI:DECOMP_01
nLabImp =1; 
MOD_TYPE = 2; % 1=ROM; 2=HROM; 0.5=Multiscale
INV_PROB = false;%true; %Inverse problem for displacement obtention, true = SI, falso = NO, for dense meshe could has large computational cost
TRU_POINT = 0 ; % Number of truncated points from the whole base
npointsINPUT = LEVELS_OF_TRUN_ENER_REG + LEVELS_OF_TRUN_ENER_DIS - 15 ; % Number of integration points used (no puede ser mayor a los modos de energia LEVELS_OF_TRUN_ENER_REG + LEVELS_OF_TRUN_ENER_DIS)

for iSTR = 1:length(LEVELS_OF_TRUN_ENER_DIS)
    
        switch MOD_TYPE
            case 0.5
                [TIEMPO_TOTAL_HF_CASES,OUTDATA] = analysis(file,path_file,nLabImp,isMICRO,MOD_TYPE);
            case 1 % FIRST REDUCTION
                [tROMI,OTHEROUT] = analysis(file,path_file,nLabImp,isMICRO,MOD_TYPE,...
                    'nModesEPS_bar',LEVELS_OF_TRUN_PHI_VAR,'nModesEPS_bar2',LEVELS_OF_TRUN_PHI_VAR2,...
                    'BASIS_DATA',BASIS_DATA,INV_PROB,TRU_POINT);
            case 2 % SECOND REDUCTION
                [TIEMPO_TOTAL_HF_CASES,OUTDATA] = analysis(file,path_file,nLabImp,isMICRO,MOD_TYPE,...
                    'NBASIS_ENER_PREV_REG',LEVELS_OF_TRUN_ENER_REG,'NBASIS_ENER_PREV_DIS',LEVELS_OF_TRUN_ENER_DIS(iSTR),...
                    'nModesEPS_bar',LEVELS_OF_TRUN_PHI_VAR,'nModesEPS_bar2',LEVELS_OF_TRUN_PHI_VAR2,'GREEDY_ALGOR',GREEDY_ALGOR_TRY,...
                    'LOAD_RESULT_GREEDY',LOAD_RESULT_GREEDY,'SET_GPS',SET_GPS,'WHERE_WEIGHTS',WHERE_WEIGHTS,'SVD_U_ws',SVD_U_ws,...
                    'SVD_STRESS_ws',SVD_STRESS_ws,'SVD_B_ws',SVD_B_ws,'SNAPSHOTS_ORIGIN',SNAPSHOTS_ORIGIN,'BASIS_DATA',BASIS_DATA,...
                    'NAME_TEST_TRAJECT_LOC',NNTTLL,'NAME_WE_TRAJE',NAME_WE_TRAJE,'ALTERNATE_BASIS_STRESS_B',ALTERNATE_BASIS_STRESS_B,...
                    'BASIS_EXPANSION',BASIS_EXPANSION,'CHATUR_NONREPEATED',CHATUR_NONREPEATED,...
                    'ECONOMIC_STORAGE',ECONOMIC_STORAGE,EnergySETS,Set_BANDAS,INV_PROB,TRU_POINT,npointsINPUT);
        end
        
end

