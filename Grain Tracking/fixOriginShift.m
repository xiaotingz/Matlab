load('lookUpTable_An4_An5.mat')
file_An4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_mesh.dream3d');
file_An5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh.dream3d');

centroids_An4 = roundn(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
centroids_An5 = roundn(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
surfG_An4 = h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
surfG_An5 = h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
origin_An4 = double(h5read(file_An4,'/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN'))';
origin_An5 = double(h5read(file_An5,'/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN'))';
centroids_An4(1,:) = [];
centroids_An5(1,:) = [];
surfG_An4(1,:) = [];
surfG_An5(1,:) = [];

diff_origin = origin_An4 - origin_An5;

% -- get corresponding centorid positions
centroids_An4InAn5Order = centroids_An4(lookUp(:,1),:);
centroids_An5tracked = centroids_An5(lookUp(:,2),:);
% -- filter out surface grains
surfG_tracked_An4 = logical(surfG_An4(lookUp(:,1)));
surfG_tracked_An5 = logical(surfG_An5(lookUp(:,2)));
surfG_eitherState = (surfG_tracked_An4 | surfG_tracked_An5);
centroids_An4InAn5Order = centroids_An4InAn5Order(surfG_eitherState,:);
centroids_An5tracked = centroids_An5tracked(surfG_eitherState,:);
% -- average origine difference
diff_centroids = centroids_An4InAn5Order - centroids_An5tracked;
aveDiffCentr = sum(diff_centroids)/length(diff_centroids);

file_newOrigin = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin.dream3d';
newOrigin = origin_An4 - aveDiffCentr;
h5write(file_newOrigin, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN', newOrigin)






