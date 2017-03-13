function reference_vector = SelectNormalREF_quad(EF_sides_SELEC,REF_VECT,coord_n_MACRO,smooth_dalpha_ELEM)


% FUNCION PARA DETERMINAR LA NORMAL AL TRACKING PATH
% **************************************************
%El REF_VECT sirve para obtener una dirección perpendicular a un plano formado por una cierta dirección sobre
%el plano x-y y vector saliente al plano x-y.
%Es una forma obtener el vector normal a una dirección en el plano x-y, pero cumple con la regla de la mano
%derecha.

switch EF_sides_SELEC
    
    case 12
        
        V12 = coord_n_MACRO(:,2)-coord_n_MACRO(:,1);
        L1 = norm(V12);
        FACTOR = abs(smooth_dalpha_ELEM(1))/(abs(smooth_dalpha_ELEM(1))+abs(smooth_dalpha_ELEM(2)));
        l_FACTOR = FACTOR*L1;
        % THETA = cos(V12(1)/L1); % angulo en radianes
        THETA = atan2(V12(2),V12(1)); 
        %if V12(1) >= 0
        %    if V12(2)<0
        %        THETA = THETA +1.5*pi;
        %    end
        %else
        %    if V12(2)>=0
        %        THETA = THETA + 0.5*pi;
        %    else
        %        THETA = THETA + pi;
        %    end
        %end
        
        COORD1 = coord_n_MACRO(:,1) + l_FACTOR*[cos(THETA);sin(THETA);0];
        
        V23 = coord_n_MACRO(:,3)-coord_n_MACRO(:,2);
        L2 = norm(V23);
        FACTOR = abs(smooth_dalpha_ELEM(2))/(abs(smooth_dalpha_ELEM(2))+abs(smooth_dalpha_ELEM(3)));
        l_FACTOR = FACTOR*L2;
        % THETA = acos(V23(1)/L2); % angulo en radianes
        THETA = atan2(V23(2),V23(1));
        
        COORD2 = coord_n_MACRO(:,2) + l_FACTOR*[cos(THETA);sin(THETA);0];
        
        
        DELTA_COORD = COORD2 - COORD1;
        
        %reference_vector = cross(DELTA_COORD/norm(DELTA_COORD),REF_VECT);        
        vector = cross(DELTA_COORD/norm(DELTA_COORD),REF_VECT);
        reference_vector = vector(1:2);
        
    case 13
        
        V12 = coord_n_MACRO(:,2)-coord_n_MACRO(:,1);
        L1 = norm(V12);
        FACTOR = abs(smooth_dalpha_ELEM(1))/(abs(smooth_dalpha_ELEM(1))+abs(smooth_dalpha_ELEM(2)));
        l_FACTOR = FACTOR*L1;
        % THETA = acos(V12(1)/L1); % angulo en radianes
        THETA = atan2(V12(2),V12(1));
        
        COORD1 = coord_n_MACRO(:,1) + l_FACTOR*[cos(THETA);sin(THETA);0];
        
        V34 = coord_n_MACRO(:,4)-coord_n_MACRO(:,3);
        L3 = norm(V34);
        FACTOR = abs(smooth_dalpha_ELEM(3))/(abs(smooth_dalpha_ELEM(3))+abs(smooth_dalpha_ELEM(4)));
        l_FACTOR = FACTOR*L3;
        % THETA = acos(V34(1)/L3); % angulo en radianes
        THETA = atan2(V34(2),V34(1));        
        
        COORD3 = coord_n_MACRO(:,3) + l_FACTOR*[cos(THETA);sin(THETA);0];
        
        
        DELTA_COORD = COORD3 - COORD1;
        
        %reference_vector = cross(DELTA_COORD/norm(DELTA_COORD),REF_VECT);
        vector = cross(DELTA_COORD/norm(DELTA_COORD),REF_VECT);
        reference_vector = vector(1:2);        
        
    case 14
        
        V12 = coord_n_MACRO(:,2)-coord_n_MACRO(:,1);
        L1 = norm(V12);
        FACTOR = abs(smooth_dalpha_ELEM(1))/(abs(smooth_dalpha_ELEM(1))+abs(smooth_dalpha_ELEM(2)));
        l_FACTOR = FACTOR*L1;
        % THETA = acos(V12(1)/L1); % angulo en radianes
        THETA = atan2(V12(2),V12(1));        
        
        COORD1 = coord_n_MACRO(:,1) + l_FACTOR*[cos(THETA);sin(THETA);0];
        
        V41 = coord_n_MACRO(:,1)-coord_n_MACRO(:,4);
        L4 = norm(V41);
        FACTOR = abs(smooth_dalpha_ELEM(4))/(abs(smooth_dalpha_ELEM(4))+abs(smooth_dalpha_ELEM(1)));
        l_FACTOR = FACTOR*L4;
        % THETA = acos(V41(1)/L4); % angulo en radianes
        THETA = atan2(V41(2),V41(1));
        
        COORD4 = coord_n_MACRO(:,4) + l_FACTOR*[cos(THETA);sin(THETA);0];
        
        
        DELTA_COORD = COORD4 - COORD1;
        
        %reference_vector = cross(DELTA_COORD/norm(DELTA_COORD),REF_VECT);
        vector = cross(DELTA_COORD/norm(DELTA_COORD),REF_VECT);
        reference_vector = vector(1:2);         
        
    case 23
        
        V23 = coord_n_MACRO(:,3)-coord_n_MACRO(:,2);
        L2 = norm(V23);
        FACTOR = abs(smooth_dalpha_ELEM(2))/(abs(smooth_dalpha_ELEM(2))+abs(smooth_dalpha_ELEM(3)));
        l_FACTOR = FACTOR*L2;
        % THETA = acos(V23(1)/L2); % angulo en radianes
        THETA = atan2(V23(2),V23(1));
        
        COORD2 = coord_n_MACRO(:,2) + l_FACTOR*[cos(THETA);sin(THETA);0];
        
        V34 = coord_n_MACRO(:,4)-coord_n_MACRO(:,3);
        L3 = norm(V34);
        FACTOR = abs(smooth_dalpha_ELEM(3))/(abs(smooth_dalpha_ELEM(3))+abs(smooth_dalpha_ELEM(4)));
        l_FACTOR = FACTOR*L3;
        % THETA = acos(V34(1)/L3); % angulo en radianes
        THETA = atan2(V34(2),V34(1));
        
        COORD3 = coord_n_MACRO(:,3) + l_FACTOR*[cos(THETA);sin(THETA);0];
        
        DELTA_COORD = COORD3 - COORD2;
        
        %reference_vector = cross(DELTA_COORD/norm(DELTA_COORD),REF_VECT);
        vector = cross(DELTA_COORD/norm(DELTA_COORD),REF_VECT);
        reference_vector = vector(1:2);         
        
    case 24
        
        V23 = coord_n_MACRO(:,3)-coord_n_MACRO(:,2);
        L2 = norm(V23);
        FACTOR = abs(smooth_dalpha_ELEM(2))/(abs(smooth_dalpha_ELEM(2))+abs(smooth_dalpha_ELEM(3)));
        l_FACTOR = FACTOR*L2;
        % THETA = acos(V23(1)/L2); % angulo en radianes
        THETA = atan2(V23(2),V23(1));
        
        COORD2 = coord_n_MACRO(:,2) + l_FACTOR*[cos(THETA);sin(THETA);0];
        
        V41 = coord_n_MACRO(:,1)-coord_n_MACRO(:,4);
        L4 = norm(V41);
        FACTOR = abs(smooth_dalpha_ELEM(4))/(abs(smooth_dalpha_ELEM(4))+abs(smooth_dalpha_ELEM(1)));
        l_FACTOR = FACTOR*L4;
        % THETA = acos(V41(1)/L4); % angulo en radianes
        THETA = atan2(V41(2),V41(1));        
        
        COORD4 = coord_n_MACRO(:,4) + l_FACTOR*[cos(THETA);sin(THETA);0];
        
        
        DELTA_COORD = COORD4 - COORD2;
        
        %reference_vector = cross(DELTA_COORD/norm(DELTA_COORD),REF_VECT);
        vector = cross(DELTA_COORD/norm(DELTA_COORD),REF_VECT);
        reference_vector = vector(1:2);         
        
    case 34
        
        V34 = coord_n_MACRO(:,4)-coord_n_MACRO(:,3);
        L3 = norm(V34);
        FACTOR = abs(smooth_dalpha_ELEM(3))/(abs(smooth_dalpha_ELEM(3))+abs(smooth_dalpha_ELEM(4)));
        l_FACTOR = FACTOR*L3;
        % THETA = acos(V34(1)/L3); % angulo en radianes
        THETA = atan2(V34(2),V34(1));        
        
        COORD3 = coord_n_MACRO(:,3) + l_FACTOR*[cos(THETA);sin(THETA);0];
        
        V41 = coord_n_MACRO(:,1)-coord_n_MACRO(:,4);
        L4 = norm(V41);
        FACTOR = abs(smooth_dalpha_ELEM(4))/(abs(smooth_dalpha_ELEM(4))+abs(smooth_dalpha_ELEM(1)));
        l_FACTOR = FACTOR*L4;
        % THETA = acos(V41(1)/L4); % angulo en radianes
        THETA = atan2(V41(2),V41(1));        
        
        COORD4 = coord_n_MACRO(:,4) + l_FACTOR*[cos(THETA);sin(THETA);0];
        
        
        DELTA_COORD = COORD4 - COORD3;
        
        %reference_vector = cross(DELTA_COORD/norm(DELTA_COORD),REF_VECT);
        vector = cross(DELTA_COORD/norm(DELTA_COORD),REF_VECT);
        reference_vector = vector(1:2);         
        
    otherwise
        
        %reference_vector = zeros(3,1);
        reference_vector = zeros(2,1);
        
end

