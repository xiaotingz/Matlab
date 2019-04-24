%% ########################## Fix Origin ########################## 
% """
% Procedures
%     1. From Aditi's data, recalculate all grain statistics. 
%     2. Calculate avg_grain_centroid_shift from the tracked non-surface grains.
%     3. Apply avg_grain_centroid_shift to sample origin of an4.
% """
load('look_up_table_an4_an5.mat')
file_An4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/Gstats_recalc_fromAditi/An4new6_recalc.dream3d');
% file_An5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/Gstats_recalc_fromAditi/An5new6_recalc_cropforalign.dream3d');
file_An5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/Gstats_recalc_fromAditi/An5new6_recalc.dream3d');

centroids_An4 = roundn(h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
centroids_An5 = roundn(h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
surfG_An4 = h5read(file_An4,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
surfG_An5 = h5read(file_An5,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
origin_An4 = double(h5read(file_An4,'/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN'))';
origin_An5 = double(h5read(file_An5,'/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN'))';
step_size = double(h5read(file_An5,'/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/SPACING'))';
centroids_An4(1,:) = [];
centroids_An5(1,:) = [];
surfG_An4(1,:) = [];
surfG_An5(1,:) = [];

%% ##### Fix Origin By Average Grain Centroid Difference #####
diff_origin = origin_An4 - origin_An5;

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

% file_newOrigin = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/Gstats_recalc_fromAditi/An4new6_fixOrigin2.dream3d';
% new_origin = origin_An4 - aveDiffCentr;
% h5write(file_newOrigin, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN', new_origin)


%% ##### Fix Origin By Comparing Voxels #####
% shift = calcVolumeDisplacement(file_An4, file_An5, look_up_table);
diff_origin = origin_An4 - origin_An5;
origin_shift = shift .* step_size;
new_origin = origin_An5 + origin_shift;

file_newOrigin = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/Gstats_recalc_fromAditi/An4new6_fixOrigin_voxel.dream3d';
h5write(file_newOrigin, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN', new_origin)


%% ######################################### Visualize grains that are bulk in both states #########################################
% file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/Gstats_recalc_fromAditi/An4new6_recalc.dream3d');
% file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/Gstats_recalc_fromAditi/An5new6_recalc.dream3d');
% 
bulk_cell_an4 = h5read(file_an4,'/DataContainers/ImageDataContainer/CellData/FeatureIds');
bulk_cell_an5 = h5read(file_an5,'/DataContainers/ImageDataContainer/CellData/FeatureIds');
% 
% origin_an4 - origin_an5
tracked_id_an4 = look_up_table(:,1);
tracked_id_an5 = look_up_table(:,2);

bulk_grain_an4 = tracked_id_an4(tracked_bulkInBothState);
bulk_grain_an5 = tracked_id_an5(tracked_bulkInBothState);

mask = ismember(bulk_cell_an4, bulk_grain_an4);
bulk_cell_an4(mask) = 1;
bulk_cell_an4(~mask) = 0;
mask = ismember(bulk_cell_an5, bulk_grain_an5);
bulk_cell_an5(mask) = 1;
bulk_cell_an5(~mask) = 0;

h5write(file_an4, '/DataContainers/ImageDataContainer/CellData/CellsBulkInBothStates', bulk_cell_an4);
h5write(file_an5, '/DataContainers/ImageDataContainer/CellData/CellsBulkInBothStates', bulk_cell_an5);

