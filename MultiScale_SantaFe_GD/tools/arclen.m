function [du_iter,dlambda,iconv] = arclen(res,KT,Fext,dofl,Du_step,Du_prev,e_VG);

ds = e_VG.ds;
%dofl = e_VG.dofl;
istep = e_VG.istep;
iter = e_VG.iter;


%***********
 iconv=1; %*
%***********
%***********
%armadura  *
% ds=0.1; %*
%***********
iconv=1;
u=KT(dofl,dofl)\[-(res) Fext(dofl)];
du1=u(:,1);
du2=u(:,2);

Du_step_i= Du_step(dofl)+du1;

a=DOT(du2,du2);
b=2*DOT(Du_step_i,du2);
c=DOT(Du_step_i,Du_step_i)-ds^2 ;
%**********************************************
% Calcular el discirminante para ver si existen 
% soluciones reales
if b*b >=4*a*c
   d=sqrt(b*b-4*a*c);
   dlambda1=(-b+d)/(2*a);
   dlambda2=(-b-d)/(2*a);
else
   iconv=0;
   return
end
%**********************************************

%**********************************************
%Codigo para seleccionar cual de las dos 
%dlambda ocupamos "Crisfield"
dc1=Du_step_i+dlambda1*du2;
dc2=Du_step_i+dlambda2*du2;
teta1=(Dot(Du_prev(dofl),dc1))/norm(dc1);
teta2=(Dot(Du_prev(dofl),dc2))/norm(dc2);

if istep*iter==1 % Si estamos en la primera iter
   dlambda=dlambda1;
else
   dlambda=dlambda1;
   if teta2>=teta1   
      dlambda=dlambda2;  
   end
end
%**********************************************
du_iter=du1+dlambda*du2;




