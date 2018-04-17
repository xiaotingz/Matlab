% ############################################################################
% NOTES
% - assuming the non-function 0 at the first element of numNeighbors has been deleted.
% ############################################################################
file = '/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_sub1Crop.dream3d';

numNeigh = h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors');
NeighborList = h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList');
numNeigh(1) = [];

keyGrain = 21;
keyNeighbors = getNeighList(keyGrain, numNeigh, NeighborList);

commonNeighPairs = [];
edges = [];
for i = 1:numNeigh(keyGrain)
    toLook = keyNeighbors(i);
    list2 = getNeighList(toLook, numNeigh, NeighborList);
    intersection = Intersect(keyNeighbors, list2);
    tmp = [int32(ones(size(intersection))).*toLook, intersection];
    commonNeighPairs = vertcat(commonNeighPairs, tmp);
    edges = vertcat(edges, intersection);
end

try 
    numEdges = length(edges)/2; 
catch rem(length(edges),2) ~= 0
    warning('#edges ~= 2n ');
end

