function varargout = f_GrafMalla(varargin)
   
   %Se asume que se está corriendo esta función con todos los directorios del programa agregado   
   %al path.
   %Devuelve los argumentos de salida: in,m_Coord,m_SetElem,e_DatSet,e_VG.   

   %Está pensado para exportar la malla como eps y después editarla en corel draw, ya que las exportaciones o
   %impresiones del GiD como eps es muy mala (hace archivos muy grandes ya que crea objetos innecesarios).
   %También se puede utilizar para dibujar la malla en MatLab cuando se hace debug.
   %Si se quiere imprimir arriba la zona localizada: hold on,conec = zeros(e_VGMicro.nElem,4);for iSet = 1:e_VGMicro.nSet,conec(e_DatSetMicro(iSet).m_IndElemSet,:) = e_DatSetMicro(iSet).conec;end,conec = conec(e_VarAuxPG(1,iElem).m_ElemLocCalc,:);[ConecUnique,~,m_pos] = unique(conec(:));m_Coord = e_DatMatSet.xx(ConecUnique,1:2);m_NroNod = 1:size(m_Coord,1);conec = reshape(m_NroNod(m_pos),[],4);patch('Faces',conec,'Vertices',m_Coord,'FaceColor','r'),hold off
   %Para imprimir los nodos frontera [~,m_Coord] = read_data('NomArch'); hold on,plot(m_Coord(m_ConecFront2(:,1),1),m_Coord(m_ConecFront2(:,1),2),'ro'),hold off y hold on,plot(m_Coord(m_ConecFront2(:,2),1),m_Coord(m_ConecFront2(:,2),2),'bo','MarkerSize',6),hold off
   
   %Se puede llamar la función con máximo dos argumentos, el nombre de archivos (string) y un argumento que
   %puede ser una matriz o una celda, ambas con filas igual al número de set. Se puede ingresar uno o dos 
   %argumentos.
   %f_GrafMalla(archDat)
   %f_GrafMalla(v_ColorSet)
   %f_GrafMalla(archDat,v_ColorSet).
   %Donde v_ColorSet puede ser una matriz o una celda. Si es una matriz, deber un color por set, cuando es una
   %celda, se debe ingresar todos los argumentos a aplicar a la función como strings separados. En la
   %dirección de las filas se coloca las opciones para cada set, mientras que dirección de las columnas cada
   %una las opciones que se imponer al set (en el caso de un tener set con menos opciones que el otro se
   %completa con matrices vacías). Si se envía una celda con una sola sola fila, se repite las opciones en
   %todos los sets (si la cantidad de set es mayor que 1).
   if nargin==1
      if ischar(varargin{1})
         archDat = varargin{1};
         [pathArchDat,nomArchDat,ext] = fileparts(archDat);
         nomArchDat = [nomArchDat,ext];
         genColor = 1;
      else
         [nomArchDat,pathArchDat] = uigetfile('*.mfl','Definir archivo de cálculo');r
         v_ColorSet = varargin{1};
         if iscell(v_ColorSet)
            c_ColorSet = v_ColorSet;
            esCell = 1;
         elseif isnumeric(v_ColorSet)
            m_ColorSet = v_ColorSet;
            esCell = 0;
         else
            error('Graficado de malla: No está definido este tipo de dato para el argumento.')
         end
      end
   elseif nargin==2
      %El primer argumento tiene que ser siempre el nombre del archivo y el segundo.
      archDat = varargin{1};
      if ~ischar(archDat)
         error(['Graficado de malla: El primer argumento tiene que ser un string indicado el directorio del ',...
            'archivo de la malla.'])         
      end
      [pathArchDat,nomArchDat,ext] = fileparts(archDat);
      nomArchDat = [nomArchDat,ext];
      %
      v_ColorSet = varargin{2};
      if iscell(v_ColorSet)
         c_ColorSet = v_ColorSet;
         esCell = 1;
      elseif isnumeric(v_ColorSet)
         m_ColorSet = v_ColorSet;
         esCell = 0;
      else
         error('Graficado de malla: No está definido este tipo de dato para el segundo argumento.')
      end      
   else
      error('Graficado de malla: Número de argumentos no definidos.')      
   end
   
   % Lectura de los datos de la malla   
   [in,m_Coord,m_SetElem,~,~,e_DatSet,e_VG] = read_data(nomArchDat,pathArchDat);
   nSet = e_VG.nSet;
   
   if exist('genColor','var')&&genColor
      %Grises según el set
      %(Notar que el set nSet es el blanco, ya que tiene valores unitarios)
      m_ColorSet = repmat((1:nSet)'/nSet,1,3);
      esCell = 0;
   end
   
   if esCell&&size(c_ColorSet,1)==1
      c_ColorSet = repmat(c_ColorSet,nSet,1);
   elseif ~esCell&&length(m_ColorSet)==1
      m_ColorSet = repmat(m_ColorSet,nSet,1);      
   end
   
   %Patch no reinicializa la figura
   clf('reset')
   
   for iSet = 1:nSet      

      m_Conec = e_DatSet(iSet).conec;
      
      %patch('Faces',m_Conec,'Vertices',m_Coord(:,1:2),'FaceVertexCData',m_Color,'FaceColor','flat');      
      %Para un solo color por set
      %patch('Faces',m_Conec,'Vertices',m_Coord(:,1:2),'FaceColor',m_Color)      
      %Sin superficie (sólo lineas)     
      %patch('Faces',m_Conec,'Vertices',m_Coord(:,1:2),'FaceColor','none')
      if esCell
         %Configuración según la indicada en la celda c_ColorSet por cada set
         patch('Faces',m_Conec,'Vertices',m_Coord(:,1:2),c_ColorSet{iSet,:})
      else
         %Un color por set.
         patch('Faces',m_Conec,'Vertices',m_Coord(:,1:2),'FaceColor',m_ColorSet(iSet,:))
      end      
      
   end
   
   %Se amplía los ejes para que se vea bien el patch.
   ampEjes = 1.02;
   coordMax = max(m_Coord(:,1:2));
   coordMin = min(m_Coord(:,1:2));
   xLim = [coordMin(1),coordMax(1)];
   yLim = [coordMin(2),coordMax(2)];
   %Se pone el axis equal antes de cambiar los límites de los ejes porque si no cambia los límites
   %de los ejes (pareciera que cambiar los límites no cambia la igualdad de los ejes). 
   axis equal
   axis([xLim*ampEjes-mean(xLim)*(ampEjes-1),yLim*ampEjes-mean(yLim)*(ampEjes-1)])
   axis off
   
   if nargout~=0
      varargout = {in,m_Coord,m_SetElem,e_DatSet,e_VG};
   end
   
end

%Cuando el Matlab (v2011a) exporta patchs a eps, y luego se importa al CorelDraw, se obtiene para cada Patch
%del MatLab un contorno cerrado que es el borde del patch (sin relleno) y otro contorno cerrado que es el Fill
%del patch, pero que no tiene borde, o uno muy fino. Una forma de unir la patch de todos los elementos en el
%CorelDraw, para simplificar los eps, es usar Buscar y Reemplazar, usar Rellenos-Color Uniforme, eligiendo el
%color del relleno que se quiere seleccionar. Ir a Organizar-Dar Forma-Límites, y crea una curva que sigue
%todo el contorno del relleno. Volver a seleccionar los rellenos, se los borra. Queda ese contorno, que
%sacándole el borde y poniendole relleno se recupera el color de los patchs. Cuidado que usar límites, parece
%que considera si la superficie tiene agujeros interiores, los descarta.
%Se puede reducir la cantidad de nodos que tiene el objeto resultante seleccionando Herramienta de Forma y
%después haciendo click en Reducir Nodos (cuidado, verificar si no elimina nodos importantes). Para utilizar
%esta herramienta hay que seleccionar todos los nodos del objeto combinado, para para ello se aprieta
%Seleccionar todos los nodos (que está en la barra de herramientas de la Herramienta de Forma).