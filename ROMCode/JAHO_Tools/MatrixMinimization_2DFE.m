function [Jmin INTF_all bMIN] = MatrixMinimization_2DFE(WEIGHTED_MODES,INCLUDE_SVDISPL,nmodes,NAME_SAVE_BASIS_DISP,...
    NAME_SAVE_BasisR,Jmin,nmodesFINTelas,TYPE_SVD_FINT,nmodesUelas,TYPE_SVD,nmodesU,...
    INTENSITY_ELASTIC_FACTOR_FORCE,INTENSITY_ELASTIC_FACTOR_DISP,bMIN,nmodesF,VOLUME_CONSTRAINT)


if nargin == 0
    load('tmp1.mat') ;
    %     INTENSITY_ELASTIC_FACTOR_FORCE = 10 ;
    %     INTENSITY_ELASTIC_FACTOR_DISP = 10 ;
end

%nmodesF =nmodes-nmodesU ;
%dbstop('14')
INTF_all = ones(nmodesF*nmodesU,1) ;
INTF_force = ones(nmodesF,1);
INTF_disp= ones(nmodesU,1);

%



if  WEIGHTED_MODES == 1
    if  INCLUDE_SVDISPL ==1
        try
            load(NAME_SAVE_BASIS_DISP,'Sp') ;
            SSd = Sp ;
        catch
            load(NAME_SAVE_BASIS_DISP,'Sd') ;
            SSd = Sd ;
        end
        
    end
    try
        load(NAME_SAVE_BasisR,'Sp') ;
        SS = Sp ;
    catch
        load(NAME_SAVE_BasisR,'Sd') ;
        SS = Sd ;
    end
    
    
    
    % Definition of     INTF_force
    switch TYPE_SVD_FINT
        case 'ELASTIC_INELASTIC_SVD'
            FACTOR_INI = 1;
            if ~isempty(SS)
                FACTOR_INI = SS(1) ;
            end
            INTF_force(1:nmodesFINTelas) = INTENSITY_ELASTIC_FACTOR_FORCE*FACTOR_INI ;
            
            if length(INTF_force)>nmodesFINTelas
                INTF_force(nmodesFINTelas+1:nmodesF) = SS(1:nmodesF-nmodesFINTelas) ;
            end
        case 'FULLSVD'
            INTF_force  = SS(1:nmodesF)  ;
    end
    
    % Definition of     INTF_disp
    if  INCLUDE_SVDISPL == 1
        switch TYPE_SVD
            case 'ELASTIC_INELASTIC_SVD'
                INTF_disp(1:nmodesUelas) = INTENSITY_ELASTIC_FACTOR_DISP;
                if length(INTF_disp)>nmodesUelas
                    INTF_disp(nmodesUelas+1:end) = SSd(1:nmodesU-nmodesUelas)/SSd(1) ;
                end
            case 'FULLSVD'
                INTF_disp  = SSd(1:nmodesU)/SSd(1)  ;
        end
    else
        INTF_disp = ones(nmodesU,1) ;
    end
    
    
    for imode =1:nmodesF
        %
        for jdof = 1:nmodesU
            iini = (imode-1)*nmodesU + jdof;
            Jmin(iini,:) =  INTF_force(imode)*INTF_disp(jdof)*Jmin(iini,:) ;
        end
    end
    
    
    if VOLUME_CONSTRAINT == 0
        for imode =1:nmodesF
            %
            for jdof = 1:nmodesU
                iini = (imode-1)*nmodesU + jdof;
                bMIN(iini,:) =  INTF_force(imode)*INTF_disp(jdof)*bMIN(iini) ;
            end
        end
    end
    
    
end
