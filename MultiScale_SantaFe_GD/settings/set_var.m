function [xg,wg] = set_var(eltype,npg,fid)

%******************************************************************************************
%*  SETEO DE VARIABLES ASOCIADAS A REGLAS DE INTEGRACION SEGUN TIPO DE ELEMENTO FINITO    *
%*  PESOS DE GAUSS Y UBICACION DE LOS PUNTOS DE INTEGRACION                               *
%*                                                                                        *                  
%*  A.E. Huespe, P.J.Sanchez                                                              *                  
%*  CIMEC-INTEC-UNL-CONICET                                                               *                  
%******************************************************************************************

% **********************
% * VARIABLES GLOBALES *
% **********************
%global_var

% *******************************************
% * REGLA DE INTEGRACION S/TIPO DE ELEMENTO *
% *******************************************
switch eltype
   case {2,10,32}
      % TRIANGULO STANDARD DE 3 NODOS FORMULADO EN DESPLAZAMIENTOS (T1)
      if npg==1
          xg = [1/3  1/3];
          wg = 1/2; 
      elseif npg==3
          xg = [1/2  0 ; 1/2  1/2 ; 0  1/2];
          wg = [1/6;1/6;1/6];
      elseif npg==4;
          xg = [0.2 0.2 ; 0.6 0.2 ; 0.2 0.6 ; 1/3 1/3];
          wg = [25/96;25/96;25/96;-27/96];
      elseif npg==7
          a1 = 0.0597158717;
          b1 = 0.4701420641;
          a2 = 0.7974269853;
          b2 = 0.1012865073;
          xg = [b2 b2 ; b1 a1 ; a2 b2 ; b1 b1 ; b2 a2 ; a1 b1 ; 1/3 1/3];
          p135 = 0.1259391805/2;
          p246 = 0.1323941527/2;
          p7   = 0.225/2;
          wg = [p135;p246;p135;p246;p135;p246;p7];
      else
          fclose(fid);
          error('Lectura de datos: Coordenadas y pesos de PG: Incorrecta cantidad de PG para elementos tipo triángulo.');   
      end
   case {4,8,20,31,108}
      % CUADRANGULO DE 4 NODOS
      [xg,wg] = gauss_c(npg);
      wg = reshape(wg*wg',[],1);
      %Orden de la numeración de los puntos de gauss en el elemento.
      %N4 ______ N3
      %  | 3  4 |
      %  | 1  2 |
      %N1 ¯¯¯¯¯¯ N2
      %Ver la correspondencia de la numeración de los nodos (N#).
      m_indX = repmat((1:npg)',npg,1);
      m_indY = reshape(repmat(1:npg,npg,1),[],1);
      %Orden de la numeración de los puntos de gauss en el GiD es distinta (ver ejemplo para el caso de 4PGs
      %abajo), pero como no se usa los puntos de gauss internos del GiD (se ingresa las coordenadas de los PG
      %explícitamente) no debería haber problema.
      %N4 ______ N3
      %  | 4  3 |
      %  | 1  2 |
      %N1 ¯¯¯¯¯¯ N2  
      xg = [xg(m_indX),xg(m_indY)];
   case {21,22,23}
      % CUADRANGULO DE 4 NODOS con 6 PG
      [xg,wg] = gauss_c(npg);
      wg = reshape(wg*wg',[],1);
      m_indX = repmat((1:npg)',npg,1);
      m_indY = reshape(repmat(1:npg,npg,1),[],1);
      xg = [xg(m_indX),xg(m_indY)];
      if npg==1 xg=[xg;xg]; wg=[wg;wg]; end
   case 5
      % ELEMENTO DE BARRA DE 2 NODOS FORMULADA EN DESPLAZAMIENTOS (2D)
      xg = [0 0];
      wg = 1;
   case 7
      % HEXAEDRO STANDARD DE 8 NODOS FORMULADO EN DESPLAZAMIENTOS (C1)
      [xg,wg] = gauss_c(npg);
   otherwise
      fclose(fid);
      error('Lectura de datos: Coordenadas y pesos de PG: Tipo de elemento incorrecto.');
      %flags = 1 ; fclose (fid) ; return;           
end