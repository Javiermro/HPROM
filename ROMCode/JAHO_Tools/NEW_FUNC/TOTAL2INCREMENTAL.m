
% CONVERT FROM TOTAL TO INCREMENTAL THE GIVEN SNAPSHOTS 
% IN ORDER TO CHECK THE PERFORMANCE OF THE SVD BASIS

clear

N_TRAINING = 10 ;

TEST_NAME = 'OneElement';
% SNAP_TRAIN_NAME = ['SNAPSHOTS_DISP_',TEST_NAME,'_SIMU_2_000.mat'] ;
SNAP_TRAIN_NAME{1} = ['SNAPSHOTS_DISP_',TEST_NAME,'_SIMU_2_00'] ;
SNAP_TRAIN_NAME{2} = ['SNAPSHOTS_DISP_STRESS_',TEST_NAME,'_SIMU_2_00'] ;
BASE_NAME = '/home/mcaicedo/Desktop/SVD_ALT/Examples/OneElement/HROM_BASIS/DATAOUT/' ;

INIT_ID_TRAIN = 0;
% IF INIT_ID _TRAIN IS 0, YOU MUST REDUCE THE TOTAL NUMBER OF TRAINING
% TRAJECTORIES IN ORDER TO FIX THE CORREC ONE.

for i_TRAINING = INIT_ID_TRAIN:N_TRAINING-1
   
    fprintf('CONVIRTIENDO TRAYECTORIA No: %i',i_TRAINING);
    
    for i_VAR = 1:length(SNAP_TRAIN_NAME)
        
        NAME_TRAIN = [BASE_NAME,SNAP_TRAIN_NAME{i_VAR},num2str(i_TRAINING),'.mat'] ;
        
        switch i_VAR
            case 1
                
                load(NAME_TRAIN,'USNAP')
                %INCR_USNAP = zeros(USNAP) ;
                N_SNAP = size(USNAP,2) ;
                
                ID_FRONT = 2:N_SNAP;
                ID_BACK = 1:N_SNAP-1;
                
                INCR_USNAP = USNAP(:,ID_FRONT) - USNAP(:,ID_BACK) ;
                
                %INCR_USNAP = [ INCR_USNAP zeros(size(USNAP,1),1) ] ;
                INCR_USNAP = [ USNAP(:,1) INCR_USNAP ] ;
                
                USNAP = INCR_USNAP ;
                
                save(NAME_TRAIN,'USNAP')
                
                clear USNAP INCR_USNAP
                
            case 2
                
                load(NAME_TRAIN,'SIGMA_SNAP')
                %INCR_USNAP = zeros(USNAP) ;
                N_SNAP = size(SIGMA_SNAP,2) ;
                
                ID_FRONT = 2:N_SNAP;
                ID_BACK = 1:N_SNAP-1;
                
                INCR_STRESS_SNAP = SIGMA_SNAP(:,ID_FRONT) - SIGMA_SNAP(:,ID_BACK) ;
                
                %INCR_STRESS_SNAP = [ INCR_STRESS_SNAP zeros(size(INCR_STRESS_SNAP,1),1) ] ;
                INCR_STRESS_SNAP = [ SIGMA_SNAP(:,1) INCR_STRESS_SNAP ] ;
                
                SIGMA_SNAP = INCR_STRESS_SNAP ;
                
                save(NAME_TRAIN,'SIGMA_SNAP')
                
                clear SIGMA_SNAP INCR_STRESS_SNAP
                
        end
    end
    
end