% """
% Procedures
%     1. From Aditi's data, recalculate all grain statistics. 
%     2. Calculate avg_grain_centroid_shift from the tracked non-surface grains.
%     3. Apply avg_grain_centroid_shift to sample origin of an4.
% """
%% ##### Load data ##### 
load('look_up_table_an4_an5.mat')
file_an4 = ('/Volumes/XIAOTING/Ni/an0-an4/An4new6.dream3d');
% file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/Gstats_recalc_fromAditi/An5new6_recalc_cropforalign.dream3d');
file_an5 = ('/Volumes/XIAOTING/Ni/an0-an4/An5new6.dream3d');

origin_an4 = double(h5read(file_an4,'/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN'))';
origin_an5 = double(h5read(file_an5,'/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN'))';
step_size = double(h5read(file_an5,'/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/SPACING'))';

%% ##### Fix Origin By Average Grain Centroid Difference #####
diff_origin = origin_an4 - origin_an5
avg_centroid_diff = calcOriginShiftFromTrackedGrains(file_an4, file_an5, look_up_table)

% file_new_origin = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/Gstats_recalc_fromAditi/An4new6_fixOrigin2.dream3d';
% new_origin = origin_an4 - avg_centroid_diff;
% h5write(file_new_origin, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN', new_origin)

%% ##### Fix Origin By Comparing Voxels #####
% shift = calcVolumeDisplacement(file_an4, file_an5, look_up_table);
diff_origin = origin_an4 - origin_an5;
origin_shift = shift .* step_size;
new_origin = origin_an5 + origin_shift;

% file_newOrigin = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/Gstats_recalc_fromAditi/An4new6_fixOrigin_voxel.dream3d';
% h5write(file_newOrigin, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN', new_origin)


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

