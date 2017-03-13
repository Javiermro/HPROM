function f_GrafNodFront(m_Conec,m_ElLoc,m_Coord,m_ConecFront,color)
   
   %Patch no reinicializa la figura
   clf
   %Se amplía los ejes para que se vea bien el patch.
   ampEjes = 1.05;
   coordMax = max(m_Coord(:,1:2));
   coordMin = min(m_Coord(:,1:2));
   xLim = [coordMin(1),coordMax(1)];
   yLim = [coordMin(2),coordMax(2)];
   axis([xLim*ampEjes-mean(xLim)*(ampEjes-1),yLim*ampEjes-mean(yLim)*(ampEjes-1)])
   axis equal
   nEl = size(m_Conec,1);
   m_Color = [0.7,0.7,0.7];
   m_Color = repmat(m_Color,nEl,1);
   m_Color(m_ElLoc,:) = repmat([200,100,0]/255,length(m_ElLoc),1);
   patch('Faces',m_Conec,'Vertices',m_Coord(:,1:2),'FaceVertexCData',m_Color,'FaceColor','flat');
   %Para un solo color
   %patch('Faces',m_Conec,'Vertices',m_Coord(:,1:2),'FaceColor',[0.7,0.7,0.7])
   hold on
   %Se plotea los nodos extremos de los elementos de frontera
   %Nodo de partida
   %plot(m_Coord(m_ConecFront(:,1),1),m_Coord(m_ConecFront(:,1),2),'sr')
   %Nodo final
   %plot(m_Coord(m_ConecFront(:,2),1),m_Coord(m_ConecFront(:,2),2),'db')
   %pause
   %Se dibuja línea uniendo los elementos
%    for iElem = 1:length(m_ConecFront)
%       m_nodElem = m_ConecFront(iElem,2:3);
%       deltaCoordElem = m_Coord(m_nodElem(2),1:2)-m_Coord(m_nodElem(1),1:2);
%       quiver(m_Coord(m_nodElem(1),1),m_Coord(m_nodElem(1),2),deltaCoordElem(1),deltaCoordElem(2),...
%          '.','LineWidth',1.5,'k');
%          %'MaxHeadSize',0.9,'LineWidth',2.5)
%       %plot(m_Coord(m_nodElem,1),m_Coord(m_nodElem,2),'-','color',color,'linewidth',1.5)
%       %pause(0.5)
%    end
   hold off   
   
end