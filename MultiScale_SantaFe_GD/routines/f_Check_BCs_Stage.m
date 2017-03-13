function [doff_test,m_GDL] = f_Check_BCs_Stage(m_CondBord_iStage,doff_test,ndn)

LinCnd=find(m_CondBord_iStage(:,2)==10|m_CondBord_iStage(:,2)==01|m_CondBord_iStage(:,2)==11);
PerCnd=find(m_CondBord_iStage(:,2)==22);
HomCnd=find(m_CondBord_iStage(:,2)==03|m_CondBord_iStage(:,2)==30|m_CondBord_iStage(:,2)==33);

m_CondBord_LinCnd=m_CondBord_iStage(LinCnd,:); 
m_CondBord_PerCnd=m_CondBord_iStage(PerCnd,:);
m_CondBord_HomCnd=m_CondBord_iStage(HomCnd,:);

% Boundary conditions - Loads and Displacements
m_DirRestr = zeros(ndn,size(m_CondBord_LinCnd,1));
m_DirRestr(ndn,:) = m_CondBord_LinCnd(:,2);
for iGdl = 1:ndn-1
    decGdl = 10^(ndn-iGdl);
    m_DirRestrGdl = fix(m_DirRestr(ndn,:)/decGdl);
    m_DirRestr(iGdl,:) = m_DirRestrGdl;
    m_DirRestr(ndn,:) = m_DirRestr(ndn,:)-m_DirRestrGdl*decGdl;
end
for itest=1:size(m_DirRestr,2)
    m_Displtest=(m_CondBord_LinCnd(itest,1)*ndn-1:m_CondBord_LinCnd(itest,1)*ndn)';
    if doff_test(m_Displtest(m_DirRestr(:,itest)==1))~=true
        doff_test(m_Displtest(m_DirRestr(:,itest)==1))=true;
    else
        error('Condiciones de Borde: Asignacion de desplazamientos o cargas al mismo DOF!');
    end
end
%Grados de libertad de los nodos en que se le impuso la restricci�n
m_GdlNodRestr = reshape(f_DofElem(m_CondBord_LinCnd(:,1)',ndn),ndn,[]);
%Valores que se le pone a los nodos restringidos
m_ValNodRestr = m_CondBord_LinCnd(:,3:end)';

% Desplazamientos prescritos impuestos
condCte = 1;
m_IndGdlRestrCondCte = m_DirRestr==condCte;
m_GdlNodRestr=m_GdlNodRestr(:,any(m_IndGdlRestrCondCte));
m_ValNodRestr=m_ValNodRestr(:,any(m_IndGdlRestrCondCte));
m_IndGdlRestrCondCte=m_IndGdlRestrCondCte(:,any(m_IndGdlRestrCondCte));
% Restricciones debidas a desplazamientos prescriptos
m_GdlRestr = m_GdlNodRestr(m_IndGdlRestrCondCte);
m_RestrCte = m_ValNodRestr(m_IndGdlRestrCondCte);
m_GDL=m_GdlRestr(find(m_RestrCte~=0));

%doffCondCte = false(ndoft,1);
%doffCondCte(m_GdlNodRestr(m_IndGdlRestrCondCte)) = true;

%Periodic Boundary Conditions
for itest=1:size(m_CondBord_PerCnd,1)
    m_Displtest=(m_CondBord_PerCnd(itest,1)*ndn-1:m_CondBord_PerCnd(itest,1)*ndn)';
    if doff_test(m_Displtest)~=true
        doff_test(m_Displtest)=true;
    else
        error('Condiciones de Borde: Asignacion de desplazamientos o cargas al mismo DOF!');
    end
end

% Homogenization Restrictions
for itest=1:size(m_CondBord_HomCnd,1)
    m_Displtest=(m_CondBord_HomCnd(itest,1)*ndn-1:m_CondBord_HomCnd(itest,1)*ndn)';
    if doff_test(m_Displtest)~=true
        doff_test(m_Displtest)=true;
    else
        error('Condiciones de Borde: Asignacion de desplazamientos o cargas al mismo DOF!');
    end
end

end