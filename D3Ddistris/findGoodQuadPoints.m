file = ('/Users/xiaotingzhong/Desktop/Datas/synthetic/180502_CubicSingleEquiaxedOut.dream3d');
nodeTypes = double(h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType'));
TriNodes = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')).';
FL = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
numNeigh = double(h5read(file,'/DataContainers/SyntheticVolumeDataContainer/CellFeatureData/NumNeighbors'));
NeighborList = double(h5read(file,'/DataContainers/SyntheticVolumeDataContainer/CellFeatureData/NeighborList'));

GA = 22;
GB = 19;

mask = (FL(:,1) == GA & FL(:,2) == GB | FL(:,1) == GB & FL(:,2) == GA);
data = [FL, TriNodes];
faceTris = data(mask,:);

mask_Quads = (nodeTypes(faceTris(:,3:5))==4);
tmp = faceTris(:,3:5);
Quads = unique(tmp(mask_Quads));

labelGs = zeros(length(Quads), 4);
centerGrain = 22;
for i = 1:length(Quads)
    mask = any(TriNodes == Quads(i), 2);
    labels = FL(mask, :);
    mask_centerGF = any(labels == centerGrain, 2);
    tmp = labels(mask_centerGF,:);
    if length(unique(tmp)) < 4
        labelGs(i,:) = [unique(tmp)',0];
    else
        labelGs(i,:) = unique(tmp);
    end
end


