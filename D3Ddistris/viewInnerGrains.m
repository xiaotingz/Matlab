file = '/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/A_forPeriphory.dream3d';
h5write(file, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN', int32([0,0,0]));
%% 

% % ---- criterion = 'centroidPos' | 'touchingFS' | 'numFaces' | 'NN_centoridPos' | 'NN_touchingFS' ----
X=434*0.15; Y=267*0.15; Z=100*0.2;
file = '/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/A_forPeriphory.dream3d';
centroids = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids')).';
centroids(1,:) = [];
active = zeros(length(centroids), 1);

criterion = 'centroidPos';
ID_ForCal = filterGrains(criterion, file, X,Y,Z);
active(ID_ForCal) = 1;
mask_2 = (centroids(:,2) < 267*0.15/3);
active(mask_2) = 1;
sum(active) 
h5write(file, '/DataContainers/ImageDataContainer/CellFeatureData/Active', int32([0;active]'));



