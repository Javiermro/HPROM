function plotq4(xyz,conec,col)
%plotq4(xy,conec,u,s)

if nargin == 2
    col = 'r.-';
end
if size(xyz,2)==2
    xyz(1,3) = 0;
end

nelem = size(conec,1);
hold on
for e=1:nelem
    c = conec(e,[1 2 3 4 1]);
    plot3(xyz(c,1),xyz(c,2),xyz(c,3),col)
end
