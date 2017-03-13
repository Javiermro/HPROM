function  [n_tensor,m_tensor,reference_vector,inject_factor] =...
 Normal_vector_selection...
          (xx,conec,uElem,n_bifurc_MACRO,inject_factor,e_VG,e_DatElemSet,...
          e_DatMat,dN_xy,EF_sidesCutCPF,hvar_new,hvar_old)

eltype = e_DatElemSet.eltype;
      
ndn                = e_VG.ndn;
%Cuidado acá si se ingresó las coordenadas en el archivo de datos, en dos columnas sin indicar las coordenadas
%z. Al realizar esto, la columna de z viene con NaN y como en la función SelectNormalREF_quad se realizan
%productos cruz, la normal de referencia (reference_vector) viene con NaNs.
coord_n_MACRO      = f_CoordElem(xx,conec);
smooth_dalpha_ELEM = e_VG.smooth_dalpha(conec);
%sihvarpg           = e_DatSet_MICRO.sihvarpg;
%hvar_old           = reshape(hvar_old,sihvarpg,[]);
%hvar_new           = reshape(hvar_new,sihvarpg,[]);
%hvar_old           = hvar_old(:,6);
%hvar_new           = hvar_new(:,6);


% GraduN_MICRO = VPG_new(6,iElem).GraduN_MICRO;

% Condicion de selecion de la normal de bifurcacion
% 1: --> gradiente de la norma del desplazamiento ||u||
% 2: --> gradiente de la variable de suavizado (mu)
% 3: --> seleccion de la normal con base en el promedio de las bandas activas en la microcelda
% 4: --> seleccion de la normal con base en la normal al tracking path
% 5: --> seleccion de la normal con base en la normal media de las bandas activas
% 6: --> seleccion de la normal con base en el salto beta de las bandas de la microcelda
%select_type = 6;

%select_type = e_VG.SELECT_NCRIT;
  if  ~isempty(e_VG.n_selec_mode)
      n_selec_mode=e_VG.n_selec_mode;
  else
      error('n_selec_mode no definido (rutina Normal_vector_selection).')
  end

switch n_selec_mode
    case 1
        % Calculo vector gradiente ||u|| normalizado
        switch eltype
            case {4,10,11,21} % Elementos cuadrilaterales
                %grad_norm_u_vector = matrix_grad_u_quad_q1_epd(coord_n_MACRO,delta_u);
                reference_vector = matrix_grad_u_quad_q1_epd(uElem,dN_xy);
            otherwise
                error('This case is not implemented yet!');
        end
        
        inject_factor=1;
        
    case 2
        % Calculo vector normal a d_mu = 0 (crack path)
        switch eltype
            case {4,10,11} % Elementos cuadrilaterales
                grad_mu = dN_xy*smooth_alpha(conec');
                reference_vector = grad_mu/norm(grad_mu);
            otherwise
                error('This case is not implemented yet!');
        end
        
        inject_factor=1;
        
    case 3
        % seleccion de la normal con base en el promedio de las bandas activas en la microcelda
        switch eltype
            case {4,10,11} % Elementos cuadrilaterales
                
                %delta_hvar = VPG_new.hvar_MICRO-VPG_old.hvar_MICRO;
                hvar_MICRO = VPG_new.hvar_MICRO;
                conec_MICRO = data_micro.conec_MICRO;
                xx_MICRO = data_micro.xx_MICRO;
                n_side_ElemSUM = zeros(1,2);
                nElem_MICRO = data_micro.e_VG_MICRO.nElem;
                u_old_MICRO = VPG_old.u_MICRO;
                u_new_MICRO = VPG_new.u_MICRO;
                
                for iElem_MICRO = 1:nElem_MICRO
                    
                    n_side_Elem_selecc = zeros(1,2);
                    ielem_material = data_micro.e_VG_MICRO.iel(iElem_MICRO);
                    
                    % material que dania y ablanda con deformacion - ielem_material==1 la segunda condicion revisa que en los 4 PGS se este presentando
                    if (ielem_material==1) && (sum(hvar_MICRO(11:data_micro.e_VG_MICRO.sihvarpg:end,iElem_MICRO)>0)==4)
                        
                        dof           = f_DofElem(conec_MICRO(iElem_MICRO,:),ndn);
                        delta_u       = u_new_MICRO(dof)-u_old_MICRO(dof);
                        conec_Elem    = conec_MICRO(iElem_MICRO,:);
                        coord_n_MICRO = f_CoordElem(data_micro.xx_MICRO,conec_MICRO(iElem_MICRO,:));
                        
                        m_ConecFrontElem = [conec_Elem(1) conec_Elem(2);conec_Elem(2) conec_Elem(3);...
                            conec_Elem(3) conec_Elem(4);conec_Elem(4) conec_Elem(1)];
                        
                        m_n_vector = Gradu_QuadQ1(coord_n_MICRO,delta_u);
                        
                        %Calculo de la normal al elemento
                        %(esto es vï¿½lido para elementos de frontera lineales, con n=cte a lo largo de la longitud de cada elemento)
                        %Vector Tangente (tiene la longitud del elemento)
                        tvector_Elem = xx_MICRO(m_ConecFrontElem(:,2),:)-xx_MICRO(m_ConecFrontElem(:,1),:);
                        %Longitud del elemento lineal
                        le_Elem = sqrt(sum(tvector_Elem.^2,2));
                        %Vector Tangente normalizado
                        tvector_Elem = bsxfun(@rdivide,tvector_Elem,le_Elem);
                        %Vector Normal normalizado
                        %Aca se define que los elementos de frontera se ingresa en sentido horario para que la normal
                        %quede hacia afuera (no interesa en que direcciï¿½n se ingrese los elementos mientras que sea la
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
                        
                    end
                    
                    n_side_ElemSUM = n_side_ElemSUM + n_side_Elem_selecc;
                    
                end
                
                norm_vecSUM = norm(n_side_ElemSUM);
                
                if norm_vecSUM == 0; norm_vecSUM = 1; end
                reference_vector = n_side_ElemSUM'./norm_vecSUM;
                
            otherwise
                error('This case is not implemented yet!');
        end
        
        inject_factor=1;
        
    case 4
        % seleccion de la normal con base en la normal al tracking path
        %Vector [x,y,z] normal al plano x-y, que sirve para obtener un vector orientado perpendicular a otro
        %en el plano x-y mediante el producto cruz con el vector REF_VECT.
        REF_VECT = [0 0 1]';
        switch eltype
            case {4,10,11,21,22,23} % Elementos cuadrilaterales
                reference_vector = SelectNormalREF_quad(EF_sidesCutCPF,REF_VECT,coord_n_MACRO,smooth_dalpha_ELEM);
            otherwise
                error('This case is not implemented yet!');
        end
        
        inject_factor=1;
        
    case 5 
        % seleccion de la normal con base en la normal media de las bandas activas de la microcelda
        switch eltype
            case {4,10,11} % Elementos cuadrilaterales
                [reference_vector,inject_factor] = SelectNormalREF_bandMICR_quad(inject_factor,data_micro,e_VG,VPG_new,VPG_old);
            otherwise
                error('This case is not implemented yet!');
        end
        
    case 6
        % seleccion de la normal con base en la normal media del salto beta en las bandas de la microcelda
        switch eltype
            case {4,10,11,21,22,23} % Elementos cuadrilaterales
                
                e_VG_MICRO         = e_DatMat.e_VG;
                e_DatSet_MICRO     = e_DatMat.e_DatSet;
                
                [reference_vector] = SelectNormalREF_betaMICR_quad(...
                    e_VG_MICRO,e_DatMat,e_DatSet_MICRO,hvar_new(6),hvar_old(6) );
                
                inject_factor=1;

%               [reference_vector,inject_factor,GraduN_MICRO] = SelectNormalREF_betaMICR_quad(inject_factor,...
 %                   data_micro,e_VG,VPG_new,VPG_old,...
 %                   e_VG_MICRO , e_DatSet_MICRO,hvar_new,hvar_old );              
  %%              VPG_new.GraduN_MICRO = GraduN_MICRO;
                
            otherwise
                error('This case is not implemented yet!');
        end
        
end

if n_selec_mode == 6 % ESTA ALTERNATIVA TOMA COMO NORMAL DE REFERENCIA EL VECTOR BETA PROMEDIO DE LA MICROESTRUCTURA
    
    if(abs(reference_vector'*n_bifurc_MACRO(:,1)) < abs(reference_vector'*n_bifurc_MACRO(:,2)))
        n_tensor(1,1) = n_bifurc_MACRO(1,1);
        n_tensor(2,2) = n_bifurc_MACRO(2,1);
        n_tensor(4,1) = n_bifurc_MACRO(2,1);
        n_tensor(4,2) = n_bifurc_MACRO(1,1);
        
        m_tensor(1,1) = n_bifurc_MACRO(1,2);
        m_tensor(2,2) = n_bifurc_MACRO(2,2);
        m_tensor(4,1) = n_bifurc_MACRO(2,2);
        m_tensor(4,2) = n_bifurc_MACRO(1,2);
        
    else
        n_tensor(1,1) = n_bifurc_MACRO(1,2);
        n_tensor(2,2) = n_bifurc_MACRO(2,2);
        n_tensor(4,1) = n_bifurc_MACRO(2,2);
        n_tensor(4,2) = n_bifurc_MACRO(1,2);
        
        m_tensor(1,1) = n_bifurc_MACRO(1,1);
        m_tensor(2,2) = n_bifurc_MACRO(2,1);
        m_tensor(4,1) = n_bifurc_MACRO(2,1);
        m_tensor(4,2) = n_bifurc_MACRO(1,1);
        
    end
    
else 

    if(abs(reference_vector'*n_bifurc_MACRO(:,1)) > abs(reference_vector'*n_bifurc_MACRO(:,2)))
        n_tensor(1,1) = n_bifurc_MACRO(1,1);
        n_tensor(2,2) = n_bifurc_MACRO(2,1);
        n_tensor(4,1) = n_bifurc_MACRO(2,1);
        n_tensor(4,2) = n_bifurc_MACRO(1,1);
        
        m_tensor(1,1) = n_bifurc_MACRO(1,2);
        m_tensor(2,2) = n_bifurc_MACRO(2,2);
        m_tensor(4,1) = n_bifurc_MACRO(2,2);
        m_tensor(4,2) = n_bifurc_MACRO(1,2);
        
    else
        n_tensor(1,1) = n_bifurc_MACRO(1,2);
        n_tensor(2,2) = n_bifurc_MACRO(2,2);
        n_tensor(4,1) = n_bifurc_MACRO(2,2);
        n_tensor(4,2) = n_bifurc_MACRO(1,2);
        
        m_tensor(1,1) = n_bifurc_MACRO(1,1);
        m_tensor(2,2) = n_bifurc_MACRO(2,1);
        m_tensor(4,1) = n_bifurc_MACRO(2,1);
        m_tensor(4,2) = n_bifurc_MACRO(1,1);
        
    end

end
