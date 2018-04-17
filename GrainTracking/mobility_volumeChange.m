file_An4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_mesh.dream3d');
file_An5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh.dream3d');

centroids_An4 = roundn(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
D_An4 = roundn(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters'),-5).';
numNeigh_An4 = double(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
neighList_An4 = double(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
% FL_An4 = double(h5read(file_An4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'));
% curvatures_An4 = h5read(file_An4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures');
% areas_An4 = roundn(h5read(file_An4,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-8);
% triNodes_An4 = double(h5read(file_An4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'));
% NCoords__An4 = double(h5read(file_An4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'));
centroids_An5 = roundn(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
D_An5 = roundn(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters'),-5).';
numNeigh_An5 = double(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
neighList_An5 = double(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
% FL_An5 = double(h5read(file_An5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'));
% curvatures_An5 = h5read(file_An5,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures');
% areas_An5 = roundn(h5read(file_An5,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-8);
% triNodes_An5 = double(h5read(file_An5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'));
% NCoords_An5 = double(h5read(file_An5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'));

centroids_An4(1,:) = [];
D_An4(1) = [];
numNeigh_An4(1) = [];
centroids_An5(1,:) = [];
D_An5(1) = [];
numNeigh_An5(1) = [];


