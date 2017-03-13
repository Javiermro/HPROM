function [xg,wg] = gauss_c(npg)

if npg > 1
    u = (1:npg-1)./sqrt((2*(1:npg-1)).^2-1);
    [vc,xg] = eig(diag(u,-1)+diag(u,1));
    [xg,k] = sort(diag(xg));
    wg = 2*vc(1,k)'.^2;
else
    xg(1)=0;
    xg(2)=0;
    wg = 2;
end
