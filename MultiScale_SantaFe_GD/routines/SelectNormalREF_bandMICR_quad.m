function [reference_vector,inject_factor] = SelectNormalREF_bandMICR_quad(inject_factor,data_micro,e_VG,VPG_new,VPG_old)

% seleccion de la normal con base en el promedio de las bandas activas en la microcelda
% Elementos cuadrilaterales

%delta_hvar = VPG_new.hvar_MICRO-VPG_old.hvar_MICRO;
ndn = e_VG.ndn;
hvar_MICRO = VPG_new.hvar_MICRO;
conec_MICRO = data_micro.conec_MICRO;
xx_MICRO = data_micro.xx_MICRO;
n_side_ElemSUM = zeros(1,2);
area_sum = 0;
nElem_MICRO = data_micro.e_VG_MICRO.nElem;
u_old_MICRO = VPG_old.u_MICRO;
u_new_MICRO = VPG_new.u_MICRO;

area_ELEM_MICRO = data_micro.area_elementos;

for iElem_MICRO = 1:nElem_MICRO
    
    n_side_Elem_selecc = zeros(1,2);
    ielem_material = data_micro.e_VG_MICRO.iel(iElem_MICRO);
    
    % material que dania y ablanda con deformacion - ielem_material==1 la segunda condicion revisa que en los 4 PGS se este presentando
    if (ielem_material==1) && (sum(hvar_MICRO(13:data_micro.e_VG_MICRO.sihvarpg:end,iElem_MICRO)>0)==4) % (sum(hvar_MICRO(11:data_micro.e_VG_MICRO.sihvarpg:end,iElem_MICRO)>0)==4)
        
        dof           = f_DofElem(conec_MICRO(iElem_MICRO,:),ndn);
        delta_u       = u_new_MICRO(dof)-u_old_MICRO(dof);
        conec_Elem    = conec_MICRO(iElem_MICRO,:);
        coord_n_MICRO = f_CoordElem(data_micro.xx_MICRO,conec_MICRO(iElem_MICRO,:));
        
        m_ConecFrontElem = [conec_Elem(1) conec_Elem(2);conec_Elem(2) conec_Elem(3);...
            conec_Elem(3) conec_Elem(4);conec_Elem(4) conec_Elem(1)];
        
        % VECTOR GRADIENTE DE LA NORMA DE DESPLAZAMIENTOS
        [m_n_vector,~] = matrix_grad_u_quad_q1_epd(coord_n_MICRO,delta_u);
        
        %Calculo de la normal al elemento
        %(esto es v�lido para elementos de frontera lineales, con n=cte a lo largo de la longitud de cada elemento)
        
        %Vector Tangente (tiene la longitud del elemento)
        tvector_Elem = xx_MICRO(m_ConecFrontElem(:,2),:)-xx_MICRO(m_ConecFrontElem(:,1),:);
        
        %Longitud del elemento lineal
        le_Elem = sqrt(sum(tvector_Elem.^2,2));
        
        %Vector Tangente normalizado
        tvector_Elem = bsxfun(@rdivide,tvector_Elem,le_Elem);
        
        %Vector Normal normalizado
        %Aca se define que los elementos de frontera se ingresa en sentido horario para que la normal
        %quede hacia afuera (no interesa en que direcci�n se ingrese los elementos mientras que sea la
        %misma en todos ellos).
        
        n_side_Elem = [tvector_Elem(:,2),-tvector_Elem(:,1)];
        dot_pr = abs(n_side_Elem*m_n_vector);
        dot_pr_Temp = n_side_Elem*[1;0];
        for i = 1:4
            %if dot_pr(i)<0.997
            %    le_Elem(i)=0;
            %end
            if dot_pr(i)>0.95 && dot_pr_Temp(i)>0
                n_side_Elem_selecc = n_side_Elem(i,:);
            end
        end
        
        n_side_ElemSUM = n_side_ElemSUM + n_side_Elem_selecc*area_ELEM_MICRO(iElem_MICRO);
        area_sum = area_sum + area_ELEM_MICRO(iElem_MICRO);
        
    end

end

n_side_ElemSUM = n_side_ElemSUM/area_sum;

norm_vecSUM = norm(n_side_ElemSUM);

%FACTOR DE AMPLIACION DE LA INTENSIDAD DE LA DEF ENRIQUECIDA
inject_factor=norm_vecSUM;

if norm_vecSUM == 0; norm_vecSUM = 1; end
reference_vector = n_side_ElemSUM'./norm_vecSUM;



