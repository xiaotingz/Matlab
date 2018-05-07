% ############################################################################
% NOTES
% - assuming the non-function 0 at the first element of numNeighbors has been deleted.
% Method
% - for each grain, get the neighbors' neighbor list(NLNs)
% - a face is defined by the grain and one of its neighbor,
%     the number of edges on one face = the #sharedElement between NL and NLN of that neighbor
% - compute sum(#intersection of NL and all NLNs) and #grain edges is half this summation
% ############################################################################
clear
file = '/Users/xiaotingzhong/Desktop/Datas/synthetic/180502_CubicSingleEquiaxedOut.dream3d';
% file = ('/Users/xiaotingzhong/Desktop/Datas/Jan.31 Austenite/180402_austenite_recons.dream3d');
% --------------------------- V6 structure ---------------------------
% numNeigh = h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors2');
% NeighborList = h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList2');
% --------------------------- synthetic structure ---------------------------
numNeigh = double(h5read(file,'/DataContainers/SyntheticVolumeDataContainer/CellFeatureData/NumNeighbors'));
NeighborList = double(h5read(file,'/DataContainers/SyntheticVolumeDataContainer/CellFeatureData/NeighborList'));
% --------------------------- V4 structure ---------------------------
% numNeigh = h5read(file, '/VoxelDataContainer/FIELD_DATA/NumNeighbors');
% NeighborList = h5read(file, '/VoxelDataContainer/FIELD_DATA/NeighborList');
% D = h5read(file, '/VoxelDataContainer/FIELD_DATA/EquivalentDiameters');
% D(1) = [];
numNeigh(1) = [];

numEdges = zeros(length(numNeigh), 1);
for i = 1:length(numNeigh)
% for i = 1:1
    keyGrain = i;
    keyNeighbors = getNeighList(keyGrain, numNeigh, NeighborList);

%     commonNeighPairs = [];
    edges = [];
    for j = 1:numNeigh(keyGrain)
        toLook = keyNeighbors(j);
        list2 = getNeighList(toLook, numNeigh, NeighborList);
        intersection = intersect(keyNeighbors, list2);
%         tmp = [int32(ones(size(intersection))).*toLook, intersection];
%         commonNeighPairs = vertcat(commonNeighPairs, tmp);
        edges = vertcat(edges, intersection);
    end

% -- Each edge is counted twice, in theory, the real number of edges should be
% exactly half of the length of our list.
    try 
        numEdges(i) = length(edges)/2; 
    catch rem(length(edges),2) ~= 0
        warning([num2str(3),' th grain, #edges ~= 2n']);
    end
end

