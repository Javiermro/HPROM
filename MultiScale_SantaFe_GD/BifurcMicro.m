function [Bif_flag,mindetQ1]=BifurcMicro(ct,e_VG)

SOLVE_METHOD = 1 ;

if SOLVE_METHOD == 0
    
    % ************************************
    % RAYLEIGH QUOTIENT METHOD ***********a
    % ************************************
    
    % Inicializaciones de las variables del procedimiento
    % ***************************************************
    sym=1;
    p=[ 0; 1];
    n_i1=[p(1) 0 ; 0 p(2) ; 0  0 ; p(2) p(1) ];
    Qpp=n_i1'*ct*n_i1;
    [V,lambda]=eig(Qpp);
    lambda_zero = max(diag(lambda));
    a=[lambda(1,1) lambda(2,2)]';
    x1= find(min(a)==a);
    lambda_i1=lambda(x1,x1);
    N_i1=V(:,x1);
    %n_i1=[N_i1(1) 0  ; 0  N_i1(2) ; 0  0 ; N_i1(2)  N_i1(1) ];
    lambda_i=100*lambda_i1;
    %detq_i1=1.e20;
    iter=0;
    tol=1.e-8;
    max_iter = 7000;
    N_i=[rand rand ]';
    N_i=N_i/norm(N_i);
    
    % Algoritmo Iterativo
    % *******************
    
    %while( abs(lambda_i1-lambda_i)>tol && iter < 3000 )
    while( abs((lambda_i1-lambda_i)/lambda_zero)>tol && iter < max_iter )
        iter=iter+1;
        lambda_i=lambda_i1;
        %  detq_i=detq_i1;
        %  N_i2=N_i;
        N_i=N_i1;
        %  n_i=n_i1;
        %  p=[N_i1(1) N_i1(2)];
        n_i1=[N_i1(1) 0 ; 0  N_i1(2) ; 0  0 ;N_i1(2) N_i1(1)];
        Qpp=n_i1'*ct*n_i1;
        
        %%% simetrizada %%%
        if  sym
            Qpp = (Qpp+Qpp')/2;
        end
        %%%%%%%%%%%%%%%%%%%
        
        [V,lambda]=eig(Qpp);
        a=[lambda(1,1) lambda(2,2) ]';
        x1= find(min(a)==a);
        % x2= find(max(a)==a);
        lambda_i1=lambda(x1(1),x1(1));
        % detq_i1= a(1)*a(2) ;
        
        N_i1=V(:,x1(1))/norm(V(:,x1(1)));
    end
    
    % Si es no simetrico y se simetriza para el analisis, en un apartado
    % posterior se debe revisar los valores propios del tensor acustico usando
    % el tensor constitutivo no simetrico
    %if sym
    % N_i1 = real(N_i1); N_i = real(N_i);
    
    n_i1=[N_i1(1) 0 ; 0  N_i1(2) ; 0  0 ;N_i1(2) N_i1(1)];
    Qpp1=n_i1'*ct*n_i1;
    
    det1=det(Qpp1);
    
    n_i2=[N_i(1) 0 ; 0  N_i(2) ; 0  0 ;N_i(2) N_i(1)];
    Qpp2=n_i2'*ct*n_i2;
    
    det2=det(Qpp2);
    
    if det2*det1 <= 0 || (det1<0 && det2<0)
        Bif_flag = 1;
    else
        Bif_flag = 0;
    end
    
    mindetQ1 = atan2(N_i1(2),N_i1(1)) ;
    %end
    
    if (iter == max_iter)
        warning('Bifucation analysis exceeds the maximum number of iterations - ELEMENT: %i: ',e_DatSet_iSet.m_NumElem(iElem_MACRO))
    end
    
else
    
    % METODO DEL BARRIDO
    [Bif_flag,m_thetaBif1,mindetQ1] = f_CondBifct(ct,e_VG) ;
    N_i1(1,1) = cos(m_thetaBif1(1));  N_i1(2,1) = sin(m_thetaBif1(1));
    N_i(1,1)  = cos(m_thetaBif1(2));  N_i(2,1)  = sin(m_thetaBif1(2));
    
end

end