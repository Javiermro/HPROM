function SNAPfNW = ArrangeByComponentSVD(FINT,nmodesU,SSd)
% Re-arrange basis matrix, and multiply it by singular values

%if isempty(SSd)
    SSd = ones(nmodesU,1) ;
%end

nsnap = size(FINT,2) ;
M = size(FINT,1)/nmodesU ;
SNAPfint = zeros(M,nmodesU*nsnap) ;
iacum = 1;
for i=1:nsnap
    for  j=1:nmodesU
        SNAPfNW(:,iacum) = SSd(j)*FINT(j:nmodesU:end,i);
        iacum = iacum+1 ;
    end
end