function [rtrial,tension_efectiva2,theta] = rtrial_damage(MDtype,ce,eps_n1,E,nu,n)
%*************************************************************************
%* Defining damage criterion surface                                     *
%*                                                                       *
%*     MDtype=  1      : SYMMETRIC                                       *
%*     MDtype=  2      : ONLY TENSION                                    *
%*     MDtype=  3      : NON-SYMMETRIC                                   *
%*                                                                       *
%*     OUTPUT:                                                           *
%*                       rtrial                                          *
%*************************************************************************

if (MDtype==1)      %* Symmetric
    
    % rtrial= sqrt(eps_n1*ce*eps_n1');
    rtrial= sqrt(eps_n1'*ce*eps_n1);
    
    tension_efectiva2 = ce*eps_n1;
    theta = 0;
    
elseif (MDtype==2)  %* Only tension
    
%    [dir_p,sig_p]=calprinc(ce,eps_n1);
%    sign1=zeros(1,4);
%    tens_pos=(sig_p+abs(sig_p))/2;
%    tens=dir_p*tens_pos*dir_p';
    % tens=dir_p*tens_pos*inv(dir_p);
    % tens=dir_p*tens_pos\dir_p;
    % sign1(1)=tens(1,1);
    % sign1(2)=tens(2,2);
    % sign1(3)=tens(1,2);
    % sign1(4)=tens(3,3);
    
    tension_efectiva2 = ce*eps_n1;
    sign1 = f_AutoValPos2Dim(tension_efectiva2,4);


 %   sign1(1)=tens(1,1);
 %   sign1(2)=tens(2,2);
 %   sign1(3)=tens(3,3);
 %   sign1(4)=tens(1,2);
    
    rtrial= sqrt(sign1'*eps_n1);
    % rtrial= sqrt(sign1*eps_n1');
    
    tension_efectiva2 = sign1;
    theta = 0;
    
elseif (MDtype==3)  %*Non-symmetric

    theta = caltheta(ce,eps_n1);
    
    if (theta == 0) || (theta == 1)
        
        tension_efectiva2 = zeros(4,1);
        % rtrial=(theta+((1-theta)/n))*sqrt(eps_n1*ce*eps_n1');
        rtrial=(theta+((1-theta)/n))*sqrt(eps_n1'*ce*eps_n1);
        
    else
        
        dphi_dsigma = zeros(4,1);        
        rtrial=(theta+((1-theta)/n))*sqrt(eps_n1'*ce*eps_n1);
        
        tension_efectiva = ce*eps_n1;
        sigma_r = sqrt(((tension_efectiva(1,1)-tension_efectiva(2,1))/2)^2+tension_efectiva(4,1)^2);
        
        dphi_dsigma(1,1) = (1-1/n)*(1/(4*sigma_r))*(1-((tension_efectiva(1,1)-tension_efectiva(2,1))/(2*sigma_r))*...
                           ((tension_efectiva(1,1)+tension_efectiva(2,1))/(2*sigma_r)));
        dphi_dsigma(2,1) = (1-1/n)*(1/(4*sigma_r))*(1+((tension_efectiva(1,1)-tension_efectiva(2,1))/(2*sigma_r))*...
                           ((tension_efectiva(1,1)+tension_efectiva(2,1))/(2*sigma_r)));
        dphi_dsigma(3,1) = 0;
        dphi_dsigma(4,1) = -(1-1/n)*(tension_efectiva(4,1)/(4*sigma_r)^2)*((tension_efectiva(1,1)+tension_efectiva(2,1))/sigma_r);
        
        tension_efectiva2= E/((1+nu)*(1-2*nu))*(nu*(dphi_dsigma(1,1)+dphi_dsigma(2,1)).*[1 1 1 0]' + (1-nu)*dphi_dsigma);
        
    end
    
end
end

function m_SigmaPos = f_AutoValPos2Dim(m_Sigma,ntens)

   %La verificación de que si es tensor de tensiones es nulo y imponer que angle=0,sigmacri=0 no es
   %necesario ya que atan2 si es 0/0 (todas las tensiones nulas) devuelve angle=0.
   %Tampoco es necesario considerar el caso en que las dos tensiones principales son iguales
   %(TAU1==TAU2), ya que atan2 con 0/0 devuelve angle 0.
   %Es necesario usar atan2 en lugar de atan porque debe ser unívoca la solución entre los ángulos 
   %[-pi/2,pi/2] para abarcar todas las direcciones posibles de las normales (considerando que la 
   %normal en una dirección abarca ambos sentidos), ya que atan2 devuelve entre [-pi,pi], y al 
   %dividir por 2, se obtiene el rango [-pi/2,pi/2] buscado. 

   TAUXX = m_Sigma(1);
   TAUYY = m_Sigma(2);
   TAUZZ = m_Sigma(3);
   TAUXY = m_Sigma(4);

   ADDTAU = 0.5*(TAUXX+TAUYY);
   REDTAU = TAUXX-TAUYY;
   %Es más estable numéricamente usar hypot, pero no sé cuánto más lento es.
   %raiz = sqrt(0.25*REDTAU*REDTAU+TAUXY*TAUXY);
   raiz = hypot(REDTAU/2,TAUXY);
   %Tensiones principales Sigma 1 y Sigma 2 en el plano xy.
   %(en el caso que de tensión plana, SigmaZ es nula por definición, así que no influye, pero en deformación
   %plana sigmaZ es una tensión principal ya que los esfuerzos de corte sigmaXZ y sigmaYZ son nulos).
   %Lo que no se sabe es cuál es correctamente TAU1>TAU2>TAU3 para todos los modelos constitutivos.
   TAU1 = ADDTAU+raiz;
   TAU2 = ADDTAU-raiz;
   TAU3 = TAUZZ;
   %
   %Se utiliza para pasar de la base de los autovectores a la base global del problema.
   angle = 0.5*atan2(2*TAUXY,REDTAU);
   cosAng = cos(angle);
   senAng = sin(angle);
   m_SigmaPos = zeros(ntens,1);
   if TAU1>0
      m_SigmaPos(1) = TAU1*cosAng^2;
      m_SigmaPos(2) = TAU1*senAng^2;
      m_SigmaPos(4) = TAU1*(senAng*cosAng);
   end
   if TAU2>0
      m_SigmaPos(1) = m_SigmaPos(1)+TAU2*senAng^2;
      m_SigmaPos(2) = m_SigmaPos(2)+TAU2*cosAng^2;
      m_SigmaPos(4) = m_SigmaPos(4)-TAU2*(senAng*cosAng);
   end
   if TAU3>0
      m_SigmaPos(3) = TAU3;
   end
    
end