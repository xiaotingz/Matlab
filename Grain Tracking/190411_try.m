load('../data/181107_mig_piececorresp_comb.mat', 'tracked_uniqueface_an4', 'tracked_uniqueface_an5');
load('/Volumes/XIAOTING/Ni/working/181107_mig_piececorresp_comb.mat', 'face_onepiece');
tracked_uniqueface_an4 = tracked_uniqueface_an4(face_onepiece, :);
tracked_uniqueface_an5 = tracked_uniqueface_an5(face_onepiece, :);

%%
idx = 5;

x_to_y = X_to_Y_ot{idx};
face_node_info = getSingleFaceNodes(tracked_uniqueface_an4(idx, :), tracked_uniqueface_an5(idx, :));

visualizeFace(face_node_info, x_to_y);



%%
load('look_up_table_an4_an5.mat')
file_An4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin2_smooth.dream3d');
file_An5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');

centroids_An4 = roundn(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
centroids_An5 = roundn(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
surfG_An4 = h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
surfG_An5 = h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
centroids_An4(1,:) = [];
centroids_An5(1,:) = [];
surfG_An4(1,:) = [];
surfG_An5(1,:) = [];


% -- get corresponding centorid positions
centroids_An4InAn5Order = centroids_An4(look_up_table(:,1),:);
centroids_An5tracked = centroids_An5(look_up_table(:,2),:);
% -- filter out surface grains
surfG_tracked_An4 = logical(surfG_An4(look_up_table(:,1)));
surfG_tracked_An5 = logical(surfG_An5(look_up_table(:,2)));
tracked_bulkInBothState = ~(surfG_tracked_An4 | surfG_tracked_An5);
centroids_An4InAn5Order = centroids_An4InAn5Order(tracked_bulkInBothState,:);
centroids_An5tracked = centroids_An5tracked(tracked_bulkInBothState,:);
% -- average origine difference
diff_centroids = centroids_An4InAn5Order - centroids_An5tracked;
aveDiffCentr = sum(diff_centroids)/length(diff_centroids)

