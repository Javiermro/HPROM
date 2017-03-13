function [reference_vector] = ...
    SelectNormalREF_betaMICR_quad(...
    e_VG_MICRO , e_DatMatMACRO, e_DatSet_MICRO,hvar_new,hvar_old )

% seleccion de la normal con base en el promedio de las bandas activas en la microcelda
% Elementos cuadrilaterales

%delta_hvar = VPG_new.hvar_MICRO-VPG_old.hvar_MICRO;

nSet = e_VG_MICRO.nSet;
beta_ElemSUM = zeros(2,1);

for iSet=1:nSet
    eltype = e_DatSet_MICRO(iSet).e_DatElem.eltype;
    if eltype==31

        nElem_MICRO =e_DatSet_MICRO(iSet).nElem;
        conec_MICRO = e_DatSet_MICRO(iSet).conec;
        xx_MICRO    = e_DatMatMACRO.xx;
        u_old_MICRO = hvar_old.u;
        u_new_MICRO = hvar_new.u;
        m_DofElem   = e_DatSet_MICRO(iSet).m_DofElem;
        dN_xy       = e_DatSet_MICRO(iSet).dN_xy  ;
        area_ELEM_MICRO =  e_DatSet_MICRO(iSet).m_VolElem;
        
        for iElem_MICRO = 1:nElem_MICRO
            
            n_side_Elem_selecc = zeros(1,2);
            dof                = m_DofElem(:,iElem_MICRO);
            delta_u            = u_new_MICRO(dof)  -u_old_MICRO(dof);
            conec_Elem         = conec_MICRO(iElem_MICRO,:);
            m_ConecFrontElem   = [conec_Elem(1) conec_Elem(2);...
                                  conec_Elem(2) conec_Elem(3);...
                                  conec_Elem(3) conec_Elem(4);...
                                  conec_Elem(4) conec_Elem(1)];
            
            [m_n_vector,m_Grad_u] = matrix_grad_u_quad_q1_epd(delta_u,dN_xy(:,:,iElem_MICRO));
            
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
                if dot_pr(i)>0.95 && dot_pr_Temp(i)>0
                    n_side_Elem_selecc = n_side_Elem(i,:);
                end
            end
            
%            GraduN_MICRO(iElem_MICRO,:) = (m_Grad_u'*n_side_Elem_selecc')';
            
            beta_ElemSUM = beta_ElemSUM + (m_Grad_u'*n_side_Elem_selecc')*area_ELEM_MICRO(iElem_MICRO);
            %n_side_ElemSUM = n_side_ElemSUM + n_side_Elem_selecc*area_ELEM_MICRO(iElem_MICRO);
            %area_sum = area_sum + area_ELEM_MICRO(iElem_MICRO);
            
        end
        
    end
end

%n_side_ElemSUM = n_side_ElemSUM/area_sum;

%norm_vecSUM = norm(n_side_ElemSUM);

%FACTOR DE AMPLIACION DE LA INTENSIDAD DE LA DEF ENRIQUECIDA
%inject_factor=norm_vecSUM;
% inject_factor=1;

%if norm_vecSUM == 0; norm_vecSUM = 1; end
%reference_vector = n_side_ElemSUM'./norm_vecSUM;

reference_vector = beta_ElemSUM;



