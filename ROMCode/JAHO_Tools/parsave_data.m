function parsave_data(NameT,NameWSTiempoTotal,NameMISC,INPUT_TRAJ,OTHEROUT)


%Trajectory info
save(NameT,'INPUT_TRAJ');

% HF time
tiempo_TOTAL_HF      = OTHEROUT.tiempo_TOTAL_HF;
save(NameWSTiempoTotal,'tiempo_TOTAL_HF');

% Miscelaneous data
INV_VAR_PLOT      = OTHEROUT.INV_VAR_PLOT;
ITERATION_COUNTER = OTHEROUT.ITERATION_COUNTER;
ponder_factors    = OTHEROUT.ponder_factors;
STATE_STEP        = OTHEROUT.STATE_STEP;
BIFURCATION_FLAG  = OTHEROUT.BIFURCATION_FLAG;
xx                = OTHEROUT.xx;
conec             = OTHEROUT.conec;

save(NameMISC, 'INV_VAR_PLOT', 'ITERATION_COUNTER','ponder_factors','STATE_STEP','BIFURCATION_FLAG','xx','conec');

end
