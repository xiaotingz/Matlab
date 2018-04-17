% file = '/Users/xiaotingzhong/Desktop/Datas/STO_1470/V6_MNK/subset1_fullrecon_MNK.dream3d';
file = '/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub2Bad_GBCD.dream3d';
% -------------------------- V6 commond --------------------------
num_of_neigh = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors'));
size = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumElements'));
centroids = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'));
% misA = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/MisorientationList'));
% CI = double(h5read(file,'/DataContainers/ImageDataContainer/CellData/Confidence Index'));
% -------------------------- V6 commond --------------------------
% num_of_neigh = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
% misA = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/MisorientationList'));
num_of_neigh(1) = [];
size(1) = [];
centroids(:,1) = [];


[maxNN, maxNN_ind] = max(num_of_neigh);
info = [1:length(num_of_neigh);num_of_neigh;size;centroids(1,:);centroids(2,:);centroids(3,:)]';
sorted = sortrows(info,-2);
% info = [centroid;0:length(num_of_neigh)].';
% sorted = sortrows(info,[1,2,3]);

% min(misA)

% misA = sort(misA);
% misA(1:5)

% histogram(CI)
% lowCI = CI(CI<0.1);
% badRate = length(lowCI)/length(reshape(CI, [],1));

% surfaceGrain = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures'));
% 1 - sum(surfaceGrain)/length(surfaceGrain);