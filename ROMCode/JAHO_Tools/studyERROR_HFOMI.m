function studyERROR_HFOMI(path_file,file,nmodesU)
% Load worspace in which HF stresses are stored    
    nomArchsnapshot_sigma = [path_file,file(1:end-4),'_snaps_sigma_HF_TRAIN.mat'];
    load(nomArchsnapshot_sigma, 'SIGMA_SNAP') ;
    SIGMA_HF = SIGMA_SNAP ;
    nomArchsnapshot_sigmaROMI = [path_file,file(1:end-4),'_snaps_sigma.mat'];
    load(nomArchsnapshot_sigmaROMI, 'SIGMA_SNAP') ;
    SIGMA_ROMI = SIGMA_SNAP ;
    DIFFER = SIGMA_HF-SIGMA_ROMI ;
    error1 = sqrt(sum(DIFFER.*DIFFER));    
    normSIGMA_SNAP = sqrt(sum((SIGMA_HF).*(SIGMA_HF)));
    errorLOC_rel = error1./normSIGMA_SNAP ;
    tolloc = 1e-8 ;
    [aaa] =  find(normSIGMA_SNAP<tolloc) ;
    if ~isempty(aaa)
        errorLOC_rel(aaa) = error1(aaa);
    end
    coloresM =  ColoresMatrix(5,[],'ALEATORIO','YES');
    MarkerM =  MarkerMatrix(5,[],'ALEATORIO','YES');
    figure(5)
    hold on 
    hhh = plot(errorLOC_rel*100,'Color',coloresM(1,:),'Marker',MarkerM{2});
    xlabel('Snapshots matrix column'); 
    ylabel('Euclidear relative error between HF and ROMI stresses (%)');
    
    legend(hhh,['Nmod U = ',num2str(nmodesU)])
    
    %%% THE SAME FOR THE RESIDUAL 
    
    load([path_file,file(1:end-4),'_DATAOUT_ROMI.mat'],'DATAOUT');
    RESIDUO = DATAOUT.RESIDUO ;
    
   figure(6)
    hold on 
    hhh = plot(errorLOC_rel,'Color',coloresM(1,:),'Marker',MarkerM{2});
    xlabel('Time step'); 
    ylabel('Full-order residual');
    title('Full-order residual at each time step')
    hhh = plot(RESIDUO,'Color',coloresM(3,:),'Marker',MarkerM{4});
    legend(hhh,['Nmod U = ',num2str(nmodesU)])
    