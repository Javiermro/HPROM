function [Phi_Grad,fii_new,n_tensor] = node_selection_quad_q1_SDA2...
    (coord_n,npe, n_tensor,Phi_Grad,dN_MACRO,...
    ind_state_MACRO_new,fii_old, EF_sidesCutCPF)

%******************************************************************
%*                PARA ELEMENTOS CUADRILATEROS                    *
%*     Seleccion lados intersectados por la discontinuidad        *
%******************************************************************
%*  coord_n : Coordenadas de los nodos que componen el elemento.  *
%*            coord_n = [x1  x2  x3 x4                            *
%*                       y1  y2  y3 y4]                           *
%*  numeracion de nodos: 4-------------3                          *
%*                       |             |                          *
%*                       |             |                          *
%*                       |             |                          *
%*                       1-------------2                          *
%******************************************************************

fii_new = zeros(npe,1);

%if ind_state_MACRO_new ~= 2
        
        % Coordenadas del nodo central.
        xcg = sum(coord_n(1,:))/npe;
        ycg = sum(coord_n(2,:))/npe;
        elem_centroid = [xcg; ycg; 0];
        
        % Seleccion de los nodos solitarios por elemento finito
        switch EF_sidesCutCPF
            case 12
                dir_vec = coord_n(:,2)-elem_centroid;
                dotpr = [n_tensor(1,1) n_tensor(2,2)]*dir_vec(1:2);
                if dotpr > 0
                    fii_new(2)=1;
                else
                    fii_new(1)=1;
                    fii_new(3)=1;
                    fii_new(4)=1;
                end
            case 13
                dir_vec = (coord_n(:,2)+coord_n(:,3))./2;
                dir_vec = dir_vec - elem_centroid;
                dotpr = [n_tensor(1,1) n_tensor(2,2)]*dir_vec(1:2);
                if dotpr > 0
                    fii_new(2)=1;
                    fii_new(3)=1;
                else
                    fii_new(1)=1;
                    fii_new(4)=1;
                end
            case 14
                dir_vec = coord_n(:,1)-elem_centroid;
                dotpr = [n_tensor(1,1) n_tensor(2,2)]*dir_vec(1:2);
                if dotpr > 0
                    fii_new(1)=1;
                else
                    fii_new(2)=1;
                    fii_new(3)=1;
                    fii_new(4)=1;
                end
            case 23
                dir_vec = coord_n(:,3)-elem_centroid;
                dotpr = [n_tensor(1,1) n_tensor(2,2)]*dir_vec(1:2);
                if dotpr > 0
                    fii_new(3)=1;
                else
                    fii_new(1)=1;
                    fii_new(2)=1;
                    fii_new(4)=1;
                end
            case 24
                dir_vec = (coord_n(:,1)+coord_n(:,2))./2;
                dir_vec = dir_vec - elem_centroid;
                dotpr = [n_tensor(1,1) n_tensor(2,2)]*dir_vec(1:2);
                if dotpr > 0
                    fii_new(1)=1;
                    fii_new(2)=1;
                else
                    fii_new(3)=1;
                    fii_new(4)=1;
                end
            case 34
                dir_vec = coord_n(:,4)-elem_centroid;
                dotpr = [n_tensor(1,1) n_tensor(2,2)]*dir_vec(1:2);
                if dotpr > 0
                    fii_new(4)=1;
                else
                    fii_new(1)=1;
                    fii_new(2)=1;
                    fii_new(3)=1;
                end
%             otherwise
%                 fii_new=fii_old;
        end
%      else
%         fii_new(:) = fii_old(:);
%     end
    %     Gradient of Phi:
    % *************************
    Phi_Grad(1,1) = dN_MACRO(1,:)*fii_new(:);  
    Phi_Grad(2,2) = dN_MACRO(2,:)*fii_new(:);         
    Phi_Grad(4,1) = Phi_Grad(2,2);
    Phi_Grad(4,2) = Phi_Grad(1,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%   Phi_Grad(:,:,iElem)=Phi_Grad_assessment(fii_new(:,iElem),dN_MACRO(:,:,iElem),e_VG);
    
%     if (n_tensor(1,1) *Phi_Grad(1,1) +n_tensor(2,2) *Phi_Grad(2,2) < 0)
%         n_tensor(:,:)=- n_tensor(:,:);
%     end
    
%end
