
function SET_SEL = DEFINE_SET(STATE_STEP,STATE_BIF,ADV_SEL)

POS(1) = find(STATE_STEP==1,1);
POS(2) = find(STATE_BIF==1,1);
SET_SEL = zeros(1,size(STATE_STEP,2));
N_SEL = length(ADV_SEL) ;

SET_TOTAL = [] ;
for iSEL = 1:N_SEL
    i_SET = ADV_SEL{iSEL};
    iTEST=zeros(size(SET_SEL));
    switch i_SET{3}
        case 'BEGIN'
            SET_SEL(1,1:i_SET{1}) = 1;
            i_TEST(1,1:i_SET{1}) = 1;
        case 'END'
            SET_SEL(1,size(STATE_BIF,2)-i_SET{1}+1:size(STATE_BIF,2)) = 1;
            iTEST(1,size(STATE_BIF,2)-i_SET{1}+1:size(STATE_BIF,2)) = 1;
        case 'YIELD'
            switch i_SET{2}
                case 'BEFORE'
                    SET_SEL(1,(POS(1)-i_SET{1}):(POS(1)-1)) = 1;
                    iTEST(1,(POS(1)-i_SET{1}):(POS(1)-1)) = 1;
                case 'AFTER'
                    SET_SEL(1,POS(1):(POS(1)+i_SET{1}-1)) = 1;
                    iTEST(1,POS(1):(POS(1)+i_SET{1}-1)) = 1;
            end
        case 'BIF'
            switch i_SET{2}
                case 'BEFORE'
                    SET_SEL(1,(POS(2)-i_SET{1}):(POS(2)-1)) = 1;
                    iTEST(1,(POS(2)-i_SET{1}):(POS(2)-1)) = 1;
                case 'AFTER'
                    SET_SEL(1,POS(2):(POS(2)+i_SET{1}-1)) = 1;
                    iTEST(1,POS(2):(POS(2)+i_SET{1}-1)) = 1;
            end
        case 'UNIF'
            switch i_SET{2}
                case 'BEGTY'
                    VECT = 1:POS(1)-1;
                    SNAP_SHOW = linspace(1,size(VECT,2),i_SET{1});
                    SNAP_SHOW = floor(SNAP_SHOW) ;
                    SET_SEL(1,VECT(SNAP_SHOW)) = 1;
                    iTEST(1,VECT(SNAP_SHOW)) = 1;
                case 'YTB'
                    VECT = POS(1)+1:POS(2)-1;
                    SNAP_SHOW = linspace(1,size(VECT,2),i_SET{1});
                    SNAP_SHOW = floor(SNAP_SHOW) ;
                    SET_SEL(1,VECT(SNAP_SHOW)) = 1;
                    iTEST(1,VECT(SNAP_SHOW)) = 1;
                case 'BTEND'
                    VECT = POS(2)+1:size(SET_SEL,2);
                    SNAP_SHOW = linspace(1,size(VECT,2),i_SET{1});
                    SNAP_SHOW = floor(SNAP_SHOW) ;
                    SET_SEL(1,VECT(SNAP_SHOW)) = 1;
                    iTEST(1,VECT(SNAP_SHOW)) = 1;
            end
    end
    SET_TOTAL = [SET_TOTAL; iTEST] ;
end

CHECK_SETS = sum(SET_TOTAL) ;

if  ~isempty(find(CHECK_SETS>1,1))
    warning('There are repeated snapshots in the selection scheme!!') ;
    %SET_SEL(find(CHECK_SETS~=0)) = 1;
end

SET_SEL = find(CHECK_SETS~=0) ;

end