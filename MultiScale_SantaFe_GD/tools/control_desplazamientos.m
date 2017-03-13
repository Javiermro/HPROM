function [du_iter,dlambda,iconv] = control_desplazamientos(res,KT,Fext_uimp,...
                       vfix,dofl,doff,Du_step,Du_prev,e_VG)


ds = e_VG.ds;
ndn = e_VG.ndn;
%***********
iconv=1; %*
%***********
node_controlled1      = e_VG.NODE1 ;
dof_local_controlled1 = e_VG.DOF1;
total_dof_controlled1 = (node_controlled1-1)*ndn+dof_local_controlled1;
%Se busca la posición (número de índice) en el vector (dofl x 1) del grado de libertad 1 a ser controlado.
free_dof_controlled1  = sum(dofl(1:total_dof_controlled1));

u=KT\[-(res) -Fext_uimp ] ;
du1 = u(:,1);
du2 = u(:,2);


%u=KT(dofl,dofl)\[-(res) Fext(dofl)];
%du1=u(:,1);
%du2=u(:,2);

%Armadura
%dlambda=(ds-du1(1))/du2(1);
Du_s=Du_step(total_dof_controlled1);
du_res= du1(free_dof_controlled1)  ;
du_fex= du2(free_dof_controlled1)  ;

if  isfield(e_VG, 'NODE2')
    node_controlled2      = e_VG.NODE2;
    dof_local_controlled2 = e_VG.DOF2;
    total_dof_controlled2 = (node_controlled2-1)*ndn+dof_local_controlled2;
    free_dof_controlled2  = sum(dofl(1:total_dof_controlled2));
    Du_s   = Du_s    - Du_step(total_dof_controlled2)   ;
    du_res = du_res  - du1(free_dof_controlled2);
    du_fex = du_fex  - du2(free_dof_controlled2);
end


%Cantiliver
%dlambda=(ds- (Du_s+du_res) )/du_fex;
dlambda=(ds- (Du_s + du_res) )/du_fex;

du_iter=du1+dlambda*du2;






