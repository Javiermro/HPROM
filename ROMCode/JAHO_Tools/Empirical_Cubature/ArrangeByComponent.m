function SNAPfNW = ArrangeByComponent(FINT,nmodesU) 

nsnap = size(FINT,2) ;
    M = size(FINT,1)/nmodesU ; 
    SNAPfint = zeros(M,nmodesU*nsnap) ;
    iacum = 1;
    for i=1:nsnap
        for  j=1:nmodesU
            SNAPfNW(:,iacum) = FINT(j:nmodesU:end,i);
            iacum = iacum+1 ;
        end
    end