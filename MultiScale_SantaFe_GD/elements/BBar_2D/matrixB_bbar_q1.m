function [B,detJ] = matrixB_bbar_q1(coord_n,xg,e_DatElemSet,BH,e_VG)

ntens = e_VG.ntens;
dofpe = e_DatElemSet.dofpe;

% Matriz Bs (standard)
% ииииииииииииииииииии
[Bs,detJ] = matrixBs_bbar_q1(coord_n,xg,e_DatElemSet,e_VG);

% Matriz B corregida:
% иииииииииииииииииии
B      = zeros(ntens,dofpe);

B1     = [Bs(1,1) 0 ; 0 Bs(2,2) ; 0 0];
B2     = [Bs(1,3) 0 ; 0 Bs(2,4) ; 0 0];
B3     = [Bs(1,5) 0 ; 0 Bs(2,6) ; 0 0];
B4     = [Bs(1,7) 0 ; 0 Bs(2,8) ; 0 0];

B(1,:) = [2/3*B1(1,1) -1/3*B1(2,2) ...
          2/3*B2(1,1) -1/3*B2(2,2)...
          2/3*B3(1,1) -1/3*B3(2,2) ...
          2/3*B4(1,1) -1/3*B4(2,2)];
  

B(2,:) = [-1/3*B1(1,1) 2/3*B1(2,2) ...
          -1/3*B2(1,1) 2/3*B2(2,2) ...
          -1/3*B3(1,1) 2/3*B3(2,2) ...
          -1/3*B4(1,1) 2/3*B4(2,2)];
  
B(3,:) = [-1/3*B1(1,1) -1/3*B1(2,2) ...
          -1/3*B2(1,1) -1/3*B2(2,2) ...        
          -1/3*B3(1,1) -1/3*B3(2,2) ...        
          -1/3*B4(1,1) -1/3*B4(2,2)];
 

B(1:3,:)     = B(1:3,:) + 1/3*BH;
B(4:ntens,:) = Bs(4:ntens,:);