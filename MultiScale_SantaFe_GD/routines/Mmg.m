function M = Mmg(inn,xx,iel,conec,Eprop,e_VG)

    %global nElem nnod npe eltype M iel
    nElem = e_VG.nElem;
    nnod = e_VG.nnod;
    npe = e_VG.npe;
    eltype = e_VG.eltype;
    %iel = e_VG.iel;

    % Inicializaciones:
    % иииииииииииииииии
    M = zeros(nnod,1);

    % Ensamble:
    % ииииииии
    for i = 1:nElem;

       coord_n = f_CoordElem(xx,conec(i,:));

       ielem_material = iel(i);

       switch eltype
           case 2
               vol = volume_t1(coord_n,Eprop(ielem_material,2),e_VG);
           case 4
               vol = volume_q1(coord_n,Eprop(ielem_material,2),e_VG);
           case 5
               vol = volume_barra2D(coord_n,Eprop(ielem_material,2));
           case 7
               vol = volume_bbar_hexahedral(coord_n,e_VG);
           case 8
               vol = volume_bbar_q1(coord_n,Eprop(ielem_material,2),e_VG);
       end
       
       for j = 1:length(conec(i,:))
           M(conec(i,j)) = M(conec(i,j))+vol/npe;  
       end

    end
    
end