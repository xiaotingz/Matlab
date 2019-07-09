function avg_centroid_diff = calcOriginShiftFromTrackedGrains(file_a, file_b, look_up_table)
% ##################################################################
% * Notes
%     - This function calcualtes average grain centroid shift from the inner tracked grains.
%     - The results can be cross compared with results with calcVolumeDisplacement.m
% * Input
%      - look_up_table = [id_file_a, id_file_b]
% ##################################################################
% ------------------------------------------------------
% load('look_up_table_a_b.mat')
% file_a = ('/Volumes/XIAOTING/Ni/an0-an4/An4new6.dream3d');
% file_b = ('/Volumes/XIAOTING/Ni/an0-an4/An5new6.dream3d');
% ------------------------------------------------------

centroids_a = roundn(h5read(file_a,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
centroids_b = roundn(h5read(file_b,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
surf_grain_a = h5read(file_a,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
surf_grain_b = h5read(file_b,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
centroids_a(1,:) = [];
centroids_b(1,:) = [];
surf_grain_a(1,:) = [];
surf_grain_b(1,:) = [];

% -- get corresponding centorid positions
centroids_a_in_b_order = centroids_a(look_up_table(:,1),:);
centroids_b_tracked = centroids_b(look_up_table(:,2),:);

% -- filter out surface grains
surf_tracked_a = logical(surf_grain_a(look_up_table(:,1)));
surf_tracked_b = logical(surf_grain_b(look_up_table(:,2)));
mask_bulk_both_states = ~(surf_tracked_a | surf_tracked_b);
centroids_a_in_b_order = centroids_a_in_b_order(mask_bulk_both_states,:);
centroids_b_tracked = centroids_b_tracked(mask_bulk_both_states,:);

% -- average origine difference
diff_centroids = centroids_a_in_b_order - centroids_b_tracked;
avg_centroid_diff = sum(diff_centroids)/length(diff_centroids);

end












