file = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub1_recons.dream3d');
ntype = double(h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType'));
tri = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')).';
fl = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
numNeigh = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors'));
NeighborList = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));

centerGrain = 407;

% ###### correct quad points #####
mask = (any(fl==centerGrain, 2) & all(fl>0, 2));
data = [fl, tri];
faceTris = data(mask,:);

mask_Quads = (ntype(faceTris(:,3:5))==4);
tmp = faceTris(:,3:5);
Quads = unique(tmp(mask_Quads));

labelGs = zeros(length(Quads), 4);

for i = 1:length(Quads)
    mask = any(tri == Quads(i), 2);
    labels = fl(mask, :);
    mask_centerGF = any(labels == centerGrain, 2);
    tmp = labels(mask_centerGF,:);
    if length(unique(tmp)) < 4
        labelGs(i,:) = [unique(tmp)',0];
    else
        labelGs(i,:) = unique(tmp);
    end
end
wrongQuads = Quads(any(labelGs==0,2),:);

ntype(wrongQuads) = 3;


