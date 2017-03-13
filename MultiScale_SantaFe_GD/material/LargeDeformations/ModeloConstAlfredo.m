E= [ 1.2 0.2 ; 0.3 4];
[lam, V]=Eigs(E)=
Sig_ppal= ce*[log(lam(1)) 0 ; 0 log(lam(2))];
Sigma= Sig_ppal(1) * V(:,1)'*V(:,1) + Sig_ppal(2) * V(:,2)'*V(:,2) 