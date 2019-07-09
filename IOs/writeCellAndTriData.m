% """
% This script rewrites a cell data and a triangle data to helper Paraview
% plot 'collection of grains' and 'cellular network'
% """

% ################### Filter Grains From Grain Centroid ###################
file_read = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
file_write = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth_curvSwapped_forParaview.dream3d';
dim = double(h5read(file_read,'/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS'));
grain_id = squeeze(h5read(file_read,'/DataContainers/ImageDataContainer/CellData/FeatureIds'));
origin = double(h5read(file_read,'/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/ORIGIN'));
step = double(h5read(file_read,'/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/SPACING'));
centroids = double(h5read(file_read,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'))';
centroids(1, :) = [];

% """
% equation of plane normal
% https://www.mathworks.com/matlabcentral/answers/75952-position-of-points-relative-to-a-plane
% """
coord_max = origin + dim.*step;
p = [origin(1), coord_max(2), coord_max(3)];
normal = [1, 1, 2];
normal = normal / norm(normal);

mask_grain_pos = ones(size(centroids, 1), 1);
for i = 1:size(centroids, 1)
    p_grain = centroids(i, :);
    if dot(p_grain - p, normal) < 0
        mask_grain_pos(i) = 0;
    end
end
mask_grain_pos = [-1; mask_grain_pos];
sum(mask_grain_pos)

mask_cell = mask_grain_pos(grain_id + 1);
mask_cell = reshape(mask_cell, 1, size(mask_cell, 1), size(mask_cell, 2), size(mask_cell, 3));

% h5write(file_write, '/DataContainers/ImageDataContainer/CellData/mask_pos', mask_cell);

%% ################### Rewrite triangle color based on FaceFeatureID ###################
file_read = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
file_write = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth_curvSwapped_forParaview.dream3d';

face_id = double(h5read(file_read,'/DataContainers/TriangleDataContainer/FaceData/FeatureFaceId'))';
% fl = h5read(file_read,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
feature_face_label = double(h5read(file_read,'/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels'))';
feature_face_label(1, :) = [];
num_faces = length(feature_face_label);
feature_face_id = (1:num_faces)';
feature_face_id_rand = feature_face_id(randperm(num_faces));

mask_surf_1 = all(feature_face_label<=0, 2);
feature_face_id_rand(mask_surf_1) = -2;
mask_surf_2 = sum(feature_face_label<=0, 2) == 1;
feature_face_id_rand(mask_surf_2) = -1;

face_id_rand = feature_face_id_rand(face_id);
h5write(file_write, '/DataContainers/TriangleDataContainer/FaceData/random_face_id', face_id_rand');








