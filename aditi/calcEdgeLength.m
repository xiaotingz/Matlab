%calculating edge lengths of each triangle
function y = calcEdgeLength(Edge,coord)
%distance between vertices
for i =1:length(Edge)
%          Edgelength(i,1)=MyDist(coord(Edge(i,1)+1,1:3),coord(Edge(i,2)+1,1:3));%1st and 2nd vertices
    Edgelength(i, 1)=norm(coord(Edge(i,1)+1,1:3) - coord(Edge(i,2)+1,1:3));
end
y=[Edge,Edgelength];
end