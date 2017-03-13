function [m_LinCond,m_CteCond,m_InvCRR,doff,dofl,doffCondCte] = f_CondBord(m_CondBordLeid,xx,...
   ndoft,ndn,ndime,nSet,e_DatSet,m_ConecFront,RotacionGlobal,m_CteMinRestr)

   sinCond = 0;
   condCte = 1;
   condIgDesp = 2;
   condHomogDef = 3;
      
   m_CondBanda=[];
   m_CondBord=[];
   i_CondBanda=0;
   i_CondBorde=0;
   CondBand   =0;
   for i_CB= 1:size(m_CondBordLeid,1)
       if m_CondBordLeid(i_CB,2) == 55 || m_CondBordLeid(i_CB,2) == 51|| m_CondBordLeid(i_CB,2) == 15
             i_CondBanda=i_CondBanda+1;
             m_CondBanda(i_CondBanda,:)= m_CondBordLeid(i_CB,:);
             CondBand =1;
       else
             i_CondBorde=i_CondBorde+1;
             m_CondBord(i_CondBorde,:)= m_CondBordLeid(i_CB,:);
       end
   end
   m_CondBordLeid= m_CondBord;
   
   
   %Para condiciones de borde variable (segunda frontera)
   condHomogDef2 = 4;
   if iscell(m_ConecFront)
      m_ConecFront2 = m_ConecFront{2};
      m_ConecFront = m_ConecFront{1};
   end
   
   %% Tratamiento de los grados de libertad de las restricciones impuestas
   %Se ordena en la filas ndn y en las columnas la cantidad de nodos con restricci�n porque facilita
   %en algunas partes el indexado.
   %Se transforma la indicaci�n �nica de grados de libertad restringidos en columnas separadas.
   m_DirRestr = zeros(ndn,size(m_CondBordLeid,1));
   m_DirRestr(ndn,:) = m_CondBordLeid(:,2);
   for iGdl = 1:ndn-1
      decGdl = 10^(ndn-iGdl);
      m_DirRestrGdl = fix(m_DirRestr(ndn,:)/decGdl);
      m_DirRestr(iGdl,:) = m_DirRestrGdl;
      m_DirRestr(ndn,:) = m_DirRestr(ndn,:)-m_DirRestrGdl*decGdl;
   end   
   %Grados de libertad de los nodos en que se le impuso la restricci�n
   m_GdlNodRestr = reshape(f_DofElem(m_CondBordLeid(:,1)',ndn),ndn,[]);

   %Valores que se le pone a los nodos restringidos
   m_ValNodRestr = m_CondBordLeid(:,3:end)';

   %% Grados de libertad globales restringidos (f: fijos) y libres (l)
   if 0
      doff = m_GdlNodRestr(m_DirRestr~=sinCond);
      %Se hace la verificaci�n de condiciones de borde a nivel de grados de libertad en lugar de
      %nodos para que se pueda ingresar dos l�neas de condiciones de borde aplicadas a los mismos
      %nodos, pero grados de libertad distintos. Esto permite que valores aplicados a los grados de
      %libertad sean distintos en caso de periodicidad y a la vez que permite condiciones de
      %periodicidad con los de rigidez.
      if length(doff)~=length(unique(doff))
         error(['Lectura de datos: Condiciones de Borde: Determinaci�n de matrices: Los grados',...
            ' de libertad con restricci�n no deben estar repetidos.'])
      end
      m_esGdlLibre = true(ndoft,1);
      m_esGdlLibre(doff) = false;
      dofl = (1:ndoft)';
      dofl = dofl(m_esGdlLibre);
      nGdlf = length(doff);
   else
      doff = false(ndoft,1);
      %Si tiene alg�n grado de libertad repetido, no influye en la cantidad de true ingresado en la
      %matriz doff.
      doff(m_GdlNodRestr(m_DirRestr~=sinCond)) = true;
      dofl = ~doff;
      nGdlf = sum(doff);
      %Se hace la verificaci�n de condiciones de borde a nivel de grados de libertad en lugar de
      %nodos para que se pueda ingresar dos l�neas de condiciones de borde aplicadas a los mismos
      %nodos, pero grados de libertad distintos. Esto permite que valores aplicados a los grados de
      %libertad sean distintos en caso de periodicidad y a la vez que permite condiciones de
      %periodicidad con los de rigidez.
      if sum(sum(m_DirRestr~=sinCond))~=nGdlf
         error(['Lectura de datos: Condiciones de Borde: Determinaci�n de matrices: Los grados',...
            ' de libertad con restricci�n no deben estar repetidos.'])
      end
   end
   
   %% Desplazamientos prescriptos impuestos
   m_IndGdlRestrCondCte = m_DirRestr==condCte;
   [m_RestrLinCondCte,m_RestrCteCondCte,m_GdlRestrCondCte,m_GdlVincCondCte] = f_RestrCondCte(...
      m_GdlNodRestr(:,any(m_IndGdlRestrCondCte)),m_ValNodRestr(:,any(m_IndGdlRestrCondCte)),...
       m_IndGdlRestrCondCte(:,any(m_IndGdlRestrCondCte)),RotacionGlobal);
       
   %Determinaci�n de que grados de libertad son pertenecientes a la condici�n de movimientos r�gidos.
   if nargout>4
      %Solo sirve si se utiliza matrices l�gicas para llevar los grados de libertad libre y fijos.
      %Falta completar para el otro caso.
      doffCondCte = false(ndoft,1);
      doffCondCte(m_GdlNodRestr(m_IndGdlRestrCondCte)) = true;
   end
    
   %% Restricciones de igualdad de desplazamientos
   m_IndGdlRestrIgDesp = m_DirRestr==condIgDesp;
   [m_RestrLinIgDesp,m_RestrCteIgDesp,m_GdlRestrIgDesp,m_GdlVincIgDesp,m_GdlRestrCteIgDesp] = ...
      f_RestrIgDesp(m_GdlNodRestr(:,any(m_IndGdlRestrIgDesp)),...
      m_ValNodRestr(:,any(m_IndGdlRestrIgDesp)),m_IndGdlRestrIgDesp(:,any(m_IndGdlRestrIgDesp)),...
      xx(:,1:ndime),ndn,nSet,e_DatSet);
   
   %% Restricciones por homogenizaci�n de las deformaciones
   m_IndGdlRestrHom = m_DirRestr==condHomogDef;
   if any(any(m_IndGdlRestrHom))
      [m_RestrLinHom,m_RestrCteHom,m_GdlRestrHom,m_GdlVincHom,m_GdlRestrCteHom] = ...
         f_RestrCondHomog(...
         m_GdlNodRestr(:,any(m_IndGdlRestrHom)),m_ValNodRestr(:,any(m_IndGdlRestrHom)),...
         m_IndGdlRestrHom(:,any(m_IndGdlRestrHom)),xx(:,1:ndime),ndime,nSet,e_DatSet,m_ConecFront);
   else
      m_RestrLinHom = [];
      m_RestrCteHom = [];
      m_GdlRestrHom = [];
      m_GdlVincHom = [];
      m_GdlRestrCteHom = [];
   end
   
   %% Restricciones por homogenizaci�n de las deformaciones
   %m_IndGdlRestrHom = m_DirRestr==condHomogDef2;
   %Corregir lo siguiente para indicar que se impone la segunda condici�n de borde de
   %homogenizaci�n. Tambi�n corregir el parche siguiente. Hacer que las funciones se le pase un
   %vector con los gdl y los valores impuestos.
   if exist('m_ConecFront2','var')
      %Se impone cualquier valor
      m_IndGdlRestrHom = [true,false;true,true];
   else
      m_IndGdlRestrHom = false(ndn,2);
   end
   if any(any(m_IndGdlRestrHom))
      m_GdlNodRestrHom2 = zeros(ndn,2);
      m_GdlNodRestrHom2(m_IndGdlRestrHom) = 1:3;
      m_ValNodRestr = zeros(ndn,2);
      [m_RestrLinHom2,m_RestrCteHom2,m_GdlRestrHom2,m_GdlVincHom2,m_GdlRestrCteHom2] = ...
         f_RestrCondHomog(...
         m_GdlNodRestrHom2,m_ValNodRestr,...
         m_IndGdlRestrHom,xx(:,1:ndime),ndime,nSet,e_DatSet,m_ConecFront2,...
         m_CteMinRestr);
      %Se hace una iteraci�n por lo menos para asegurar que la parte correspondiente a la segunda
      %condici�n de borde de homogenizaci�n est� bien condicionada (determinante distinto de cero
      %alejado de cero). Igual puede fallar al combinarse con otras condiciones de borde.
      %Para no tener que hacer este ensamblaje, habr�a que elegir gdl arbitrarios para cada fila de
      %la matriz (la posici�n de las filas no cambian el resultado), y as� poder ensamblar sin
      %conocer los gdl fijos de esta condici�n, y luego tomar esas filas para determinar la 
      %matriz m_LinCondRestHom2.
      m_LinCondRestHom2 = sparse(m_GdlRestrHom2,m_GdlVincHom2,m_RestrLinHom2,3,ndoft);
      %Se asume que la conectividad de la frontera indicada en la primera es la frontera de la
      %malla.
      m_GdlRestMR = f_DetGdlFijoCBMR(m_LinCondRestHom2,m_ConecFront,doff,ndoft,ndn);
      doff(m_GdlRestMR) = true;
      dofl(m_GdlRestMR) = false;
      nGdlf = nGdlf+3;
      m_GdlRestrHom2 = repmat(m_GdlRestMR,length(m_GdlVincHom2)/3,1);
      %condMat = cond(full(m_LinCondRestHom2(:,m_GdlNodRestr(m_IndGdlRestrHom))));
      condMat = condest(m_LinCondRestHom2(:,m_GdlRestMR));
      fprintf('N�mero de Condici�n %f\n',condMat);
      %while 1/condMat<1e8*eps(1)    %condMat>1e-8/eps(1)
      %end
   else
      m_RestrLinHom2 = [];
      m_RestrCteHom2 = [];
      m_GdlRestrHom2 = [];
      m_GdlVincHom2 = [];
      m_GdlRestrCteHom2 = [];
   end
   
  
%    [m_RestrLinCondCte,m_RestrCteCondCte,m_GdlRestrCondCte,m_GdlVincCondCte] = f_RestrBanda(...
%       m_GdlNodRestr(:,any(m_IndGdlRestrCondCte)),m_ValNodRestr(:,any(m_IndGdlRestrCondCte)),...
%        m_IndGdlRestrCondCte(:,any(m_IndGdlRestrCondCte)));


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       

   
   %% Matrices de restricci�n
   m_LinCond = sparse([m_GdlRestrCondCte;m_GdlRestrIgDesp;m_GdlRestrHom;m_GdlRestrHom2],...
      [m_GdlVincCondCte;m_GdlVincIgDesp;m_GdlVincHom;m_GdlVincHom2],...
      [m_RestrLinCondCte;m_RestrLinIgDesp;m_RestrLinHom;m_RestrLinHom2],...
      ndoft,ndoft);
  m_CteCond = sparse([unique(m_GdlRestrCondCte);m_GdlRestrCteIgDesp;m_GdlRestrCteHom;m_GdlRestrCteHom2],...
      1,[m_RestrCteCondCte;m_RestrCteIgDesp;m_RestrCteHom;m_RestrCteHom2],ndoft,1);
  % Condensaci�n para que en funci�n de los grados de libertad sin restricci�n
  %m_InvCRR = inv(m_LinCond(doff,doff));
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  if CondBand
      
      %% Desplazamientos prescriptos impuestos por bandas de localizacion
      %
      %              Nodo 1 -------  Nodo 2       Matriz Boundary Conditions
      %                 |              |
      %                 |     Banda j  |              j    55   Nodo1  Nodo2
      %                 |              |              j    55   Nodo3  Nodo4
      %              Nodo 3 -------  Nodo 4           j    55   Nodo5  Nodo6
      %                 |              |              j    55   Nodo7  Nodo8
      %                 |    Banda j   |
      %                 |              |
      %              Nodo 5 -------  Nodo 6
      %                 |              |
      %                 |    Banda j   |
      %                 |              |
      %              Nodo 7 -------  Nodo 8
      %
      
      nband = max(m_CondBanda(:,1));
      CondBand=struct('matriz',[]);
      for i_band=1:nband
          CondBand(i_band).matriz = m_CondBanda(find(m_CondBanda(:,1)==i_band),:);
          % CondBand{i_band}.matrix= m_CondBanda(find(m_CondBanda(:,1)==i_band),:);
      end
      
      for i_band=1:nband
          node1= CondBand(i_band).matriz(1,3);
          node2= CondBand(i_band).matriz(1,4);
          doff1= [((node1-1)*ndn+1) ; node1*ndn];
          doff2= [((node2-1)*ndn+1) ; node2*ndn];
          nCond_Band= size(CondBand(i_band).matriz,1);
          
          for i_eq= 2:nCond_Band;
              
              node3= CondBand(i_band).matriz(i_eq,3);
              node4= CondBand(i_band).matriz(i_eq,4);
              doff3= [((node3-1)*ndn+1) ; node3*ndn];
              doff4= [((node4-1)*ndn+1) ; node4*ndn];
              
              doff(doff3)=1;
              nGdlf=nGdlf+2;
              m_LinCond(doff3(1),doff3(1))=1;m_LinCond(doff3(1),doff4(1))=-1;...
                  m_LinCond(doff3(1),doff2(1))=1;m_LinCond(doff3(1),doff1(1))=-1;
              
              m_LinCond(doff3(2),doff3(2))=1;m_LinCond(doff3(2),doff4(2))=-1;...
                  m_LinCond(doff3(2),doff2(2))=1;m_LinCond(doff3(2),doff1(2))=-1;
          end
      end
      dofl = ~doff;
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
   m_InvCRR = m_LinCond(doff,doff)\speye(nGdlf,nGdlf);
   %Ver si es m�s lento as� para matrices sparse que hacer dos veces la inversa.
   %m_LinCond = m_LinCond(doff,doff)\[-m_LinCond(doff,dofl),m_CteCond(doff)];
   %m_LinCond = m_LinCond(doff,doff)\-m_LinCond(doff,dofl);
   m_LinCond = m_InvCRR*-m_LinCond(doff,dofl);
   %Vector restricciones constantes de los grados de libertad fijos (doff).
   %m_CteCond = m_LinCond(:,end);
   %m_LinCond = m_LinCond(:,1:end-1);
   
end
 
%%
function esInt = f_esPtoIntElem2D(m_CoordNod,PtoInt)
   
   %Se asume que los nodos de las conectividades est�n ordenadas en sentido antihorario, es decir la
   %lista de coordenadas en m_CoordNod est�n ordenadas de esa forma. Esto es importante para la
   %verificaci�n que el punto sea interior.
   %Solo sirve cuando el elemento tiene forma convexa en la malla (en coordenadas globales).
   m_DifCoord = [diff(m_CoordNod);m_CoordNod(1,:)-m_CoordNod(end,:)];
   m_DifPtoInt = bsxfun(@minus,PtoInt,m_CoordNod);
   esInt = all((m_DifPtoInt(:,1).*m_DifCoord(:,2)-m_DifPtoInt(:,2).*m_DifCoord(:,1))<=0);
   
end

%%
function [m_RestrLin,m_RestrCte,m_GdlRestr,m_GdlVinc,m_GdlRestrCte] = f_RestrIgDesp(...
   m_GdlNodRestr,m_ValNodRestr,m_IndGdlRestr,xx,ndn,nSet,e_DatSet)

   % Restricciones de igualdad de desplazamientos
   %En forma general no se conoce la cantidad de grados de libertad que van estar vincunlados con
   %la igualdad de desplazamiento, ya que no se conoce el tipo de elemento donde coincide la
   %posici�n que se impone la misma.
   m_GdlVinc = [];
   m_GdlRestr = [];
   m_RestrLin = [];

   nNodRestr = size(m_GdlNodRestr,2);
   for iNodRestr = 1:nNodRestr
      m_Pos = m_ValNodRestr(:,iNodRestr);
      for iSet =  1:nSet
         eltype = e_DatSet(iSet).e_DatElem.eltype;
         conec = e_DatSet(iSet).conec;
         m_DofElem = e_DatSet(iSet).m_DofElem;
         %Cantidad de elementos en el Set.
         nElem = e_DatSet(iSet).nElem;
         for iElem = 1:nElem
            m_CoordNodElem = f_CoordElem(xx,conec(iElem,:))';
            if f_esPtoIntElem2D(m_CoordNodElem,m_Pos')
               switch eltype
                  case {2,32}
                     %Tri�ngulo de tres nodos con funciones de forma (para aprox. de coord.) lineales
                     m_FunFormX = [m_CoordNodElem,ones(3,1)]'\[m_Pos;1];
                  case {4,8,20,31,108}
                     %Cuadr�ngulo de cuatro nodos con funciones forma (para aprox. de coord.) 
                     %bilineales
                     m_FunFormX = [m_CoordNodElem,m_CoordNodElem(:,1).*m_CoordNodElem(:,2),...
                        ones(4,1)]'\[m_Pos;m_Pos(1)*m_Pos(2);1];
                  otherwise
                     error(['Lectura de datos: Condiciones de borde: Determinaci�n de matrices: ',...
                        'Restricciones de igualdad de desplazamiento: Tipo de elemento no ',...
                        'implementado.']);
               end
               %Se guarda como vector para usar el ensamblaje de la funci�n sparse.
               %Se considera que los grados de libertad del elemento est�n aproximados con las mismas
               %funciones de forma.
               m_GdliNodRestr = m_GdlNodRestr(m_IndGdlRestr(:,iNodRestr),iNodRestr);
               nGdlNodRestr = length(m_GdliNodRestr);
               %m_GdlElem = reshape(f_DofElem(conec(iElem,:),ndn),ndn,[]);
               m_GdlElem = reshape(m_DofElem(:,iElem),ndn,[]);
               %
               m_GdlRestr = [m_GdlRestr;repmat(m_GdliNodRestr,size(conec,2)+1,1)];
               m_GdlVinc = [m_GdlVinc;m_GdliNodRestr;reshape(m_GdlElem(m_IndGdlRestr(:,iNodRestr),...
                  :),[],1)];
               m_RestrLin = [m_RestrLin;ones(nGdlNodRestr,1);reshape(repmat(-m_FunFormX,1,...
                  nGdlNodRestr)',[],1)];
               %Se considera que el desplazamiento se eval�a s�lo con la funci�n de forma de un
               %elemento. Cuando coincide la posici�n con un lado contiguo de dos elementos se eval�a
               %solo sobre uno de ellos. (Esto no es correcto con mallas compuestas?).
               %Varias restricciones pueden estar sobre los mismos elementos
               break
            end
         end
      end
   end
   %Este no se utiliza pero deber�a devolver los valores que se ensamblan en la matriz de constante.
   %Es para si se quiere que los desplazamientos vinculados tenga una diferencia en su valor
   m_GdlRestrCte = m_GdlNodRestr(m_IndGdlRestr);
   m_RestrCte = zeros(length(m_GdlRestrCte),1);
   
end

%%
function [m_RestrLin,m_RestrCte,m_GdlRestr,m_GdlVinc] = f_RestrCondCte(m_GdlNodRestr,...
   m_ValNodRestr,m_IndGdlRestr,RotacionGlobal)
   
   angulo=(pi/180)*RotacionGlobal;
   cosAngulo=cos(angulo);
   sinAngulo=sin(angulo);

   % Restricciones debidas a desplazamientos prescriptos
   m_GdlRestr = m_GdlNodRestr(m_IndGdlRestr);
   m_GdlRestrTMP = repmat(m_GdlRestr',2,1);
   m_GdlRestr = m_GdlRestrTMP(:);
   m_RestrLin=zeros(size(m_GdlRestrTMP));
   for i=1:size(m_GdlRestrTMP,2)
       
      if mod(m_GdlRestrTMP(1,i),2)==0
          m_GdlRestrTMP(1,i)=m_GdlRestrTMP(1,i)-1;
          m_RestrLin(1,i)=-sinAngulo;
          m_RestrLin(2,i)=cosAngulo;
      else
          m_GdlRestrTMP(2,i)=m_GdlRestrTMP(2,i)+1;
          m_RestrLin(1,i)=cosAngulo;
          m_RestrLin(2,i)=sinAngulo;          
      end
   end
   
   %m_GdlVinc = m_GdlRestr;
   m_GdlVinc = m_GdlRestrTMP(:);
   m_RestrCte = m_ValNodRestr(m_IndGdlRestr);
   %m_RestrLin = ones(length(m_RestrCte),1);
   m_RestrLin = m_RestrLin(:);

end

%%
function [m_C,m_RestrCte,m_GdlRestr,m_GdlVinc,m_GdlRestrCte] = f_RestrCondHomog(...
   m_GdlNodRestr,m_ValNodRestr,m_IndGdlRestr,xx,ndime,nSet,e_DatSet,m_ConecFront,m_CteMinRestr)

    %           4   5
    %  11 o-----o---o
    %     |         |
    %     |         |
    %   6 o         |
    %     |         o 15
    %     |         |
    %   1 o--o------o 9
    %        8
    %
    % conec : lista de conectividades para pseudo-elementos de barra en la frontera, 
    % orden en sentido anti-horario. Ejemplo:
    % conec = [ 1 8 ; 8 9 ; 9 15 ; 15 5 ; 5 4 ; 4 11 ; 11 6 ; 6 1 ]
    %
    % xx    : lista de coordenadas de todos los nodos de la malla
    % xx    = [ x1 y1 z1 ; x2 y2 z2 ; x3 y3 z3 ; x4 y4 z4 ; ... ]

   %m_ValNodRestr: no se utiliza en esta restricci�n lineal.
   %Se tiene 3 ecuaciones de independientes de la restricci�n por homogenizaci�n de la deformaci�n
   %en un problema de dos dimensiones
   ngdLDep = 3;
  
   %m_ElType = [e_DatSet.eltype]';
   %if any(~(m_ElType==2|m_ElType==4|m_ElType==8))
   m_ElType = arrayfun(@(x)x.e_DatElem.eltype,e_DatSet);
   if any(~(m_ElType==2|m_ElType==4|m_ElType==8|m_ElType==20|m_ElType==31|m_ElType==32))
      error(['Lectura de datos: Condiciones de borde: Determinaci�n de matrices: Restricciones ',...
         'por homogenizaci�n de la deformaci�n: Tipo de elemento no implementado.'])
   end
   if ndime~=2
      error(['Lectura de datos: Condiciones de borde: Determinaci�n de matrices: Restricciones ',...
         'por homogenizaci�n de la deformaci�n: Implementado s�lo para problemas 2D.'])
   end
   if sum(sum(m_IndGdlRestr))~=ngdLDep
      error(['Lectura de datos: Condiciones de borde: Determinaci�n de matrices: Restricciones ',...
         'por homogenizaci�n de la deformaci�n: Est� implementado para que tenga tres ',...
         'restricciones de este tipo.'])
   end
   
   %% Determinaci�n de los grados de libertad dependientes (impuestos como restricci�n)
   m_GdlDep = m_GdlNodRestr(m_IndGdlRestr);
   
   %% Calculo de la normal al elemento
   %(esto es v�lido para elementos de frontera lineales, con n=cte a lo largo de la longitud de 
   %cada elemento)
   %Vector Tangente (tiene la longitud del elemento)
   tvector = xx(m_ConecFront(:,2),:)-xx(m_ConecFront(:,1),:);
   %Longitud del elemento lineal
   le = sqrt(sum(tvector.^2,2));
   %Vector Tangente normalizado
   tvector = bsxfun(@rdivide,tvector,le);
   %Vector Normal normalizado
   %Ac� se define que los elementos de frontera se ingresa en sentido horario para que la normal
   %quede hacia afuera pero como las restricciones de homogenizaci�n son ecuaciones homog�neas no
   %importa el signo que tenga las componentes de la normal ya que da el mismo resultado (es decir
   %no interesa en que direcci�n se ingrese los elementos mientras que sea la misma en todos ellos).
   nvector = [tvector(:,2),-tvector(:,1)];

   %% Determinacion de la matriz C global de homogenizacion
   [nelFront,nNodElFront] = size(m_ConecFront);   

   %Se considera que el n�mero de desplazamientos est� dictado por la dimensi�n del problema
   %(ndime).
   m_C = zeros(ngdLDep,nNodElFront*ndime,nelFront);
   m_GdlVinc = zeros(ngdLDep,nNodElFront*ndime,nelFront);
   
   %Se adopta en forma arbitraria que el primer grado de libertad donde se le impone restricci�n de
   %este en el programa tiene la ecuaci�n de homogenizaci�n de las componentes xx del tensor, la
   %segunda las componentes yy, y el tercero tiene las componentes xy.
   %Igual el orden no influye en que el sistema est� mal condicionada (cercano a matriz singular).
   m_GdlRestr = repmat(m_GdlDep,[1,nNodElFront*ndime,nelFront]);
   for iElemFront = 1:nelFront
    
      m_DofElemFront = f_DofElem(m_ConecFront(iElemFront,:),ndime);
    
      m_GdlVinc(:,:,iElemFront) = repmat(m_DofElemFront',ngdLDep,1);
    
      %Integral de la funci�n de forma sobre el borde del elemento (sobre el elemento de barra).
      %En el caso de elementos lineales la integral es la mitad de la longitud del elemento
      %(tri�ngulo de base le y altura 1).
      intN = le(iElemFront)/2;
      
      m_C(1,1:2:end,iElemFront) = nvector(iElemFront,1)*intN;
      m_C(2,2:2:end,iElemFront) = nvector(iElemFront,2)*intN;
      m_C(3,1:2:end,iElemFront) = nvector(iElemFront,2)*intN;
      m_C(3,2:2:end,iElemFront) = nvector(iElemFront,1)*intN;
         
   end
   
   %Se vectoriza las matrices para usarla con la funci�n sparse.
   m_GdlRestr = m_GdlRestr(:);
   m_GdlVinc = m_GdlVinc(:);
   m_C = m_C(:);
   
   %Este no se utiliza pero deber�a devolver los valores que se ensamblan en la matriz de constante.
   m_GdlRestrCte = m_GdlDep;
   if exist('m_CteMinRestr','var')
      %m_CteMinRestr: Vector conteniendo las ctes. de m�nima restricci�n ordenada seg�n las componentes xx, yy, y
      %xy.
      m_RestrCte = m_CteMinRestr(:);
   else
      m_RestrCte = zeros(length(m_GdlRestrCte),1);
   end

end