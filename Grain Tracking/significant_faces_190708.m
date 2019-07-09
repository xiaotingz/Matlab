% """
% an4 ---> initial state
% an5 ---> final state
% """
% ########################### Prepare Data ###########################
file_an4 = '/Volumes/XIAOTINGbig/Ni_an0-an4/an0.dream3d';
file_an5 = '/Volumes/XIAOTINGbig/Ni_an0-an4/an1to0.dream3d';
load('/Volumes/XIAOTINGbig/Ni_an0-an4/look_up_table_an0_an1crop.mat');
fname = 'plot_significant_an0an1.mat';

fl_unique_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
fl_unique_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';

% ----- Filter Field Of View quasi faces -----
mask_true_an4 = all(fl_unique_an4 >= 0, 2) & any(fl_unique_an4 > 0, 2);
mask_true_an5 = all(fl_unique_an5 >= 0, 2) & any(fl_unique_an5 > 0, 2);
fl_unique_an4 = fl_unique_an4(mask_true_an4, :);
fl_unique_an5 = fl_unique_an5(mask_true_an5, :);

% ----- Identify true surface grains -----
mask_true_surf_an4 = any(fl_unique_an4 == 0, 2);
surf_grain_an4_true = setdiff(unique(fl_unique_an4(mask_true_surf_an4, :)), 0);
mask_true_surf_an5 = any(fl_unique_an5 == 0, 2);
surf_grain_an5_true = setdiff(unique(fl_unique_an5(mask_true_surf_an5, :)), 0);

% ----- Track faces -----
[tracked_uniqueface_an4_all, tracked_uniqueface_an5_all] = ...
    trackUniqueFace(file_an4, file_an5, look_up_table, 'use_all_true_faces');
% """ adjust the order of faces_tracked_all in faces so that the tracked faces have one-to-one corresp """
mask_tracked_an4 = ismember(fl_unique_an4, tracked_uniqueface_an4_all, 'rows');
mask_tracked_an5 = ismember(fl_unique_an5, tracked_uniqueface_an5_all, 'rows');
fl_unique_an4(mask_tracked_an4, :) = tracked_uniqueface_an4_all;
fl_unique_an5(mask_tracked_an5, :) = tracked_uniqueface_an5_all;

% ----- Filter small faces -----
eps_curv = 1;
eps_area = 7;
eps_min_ang = 10;
face_itg_curv_an4 = calcFaceItgCurv(file_an4, fl_unique_an4, 'as_given', eps_curv, eps_area, eps_min_ang);
face_itg_curv_an5 = calcFaceItgCurv(file_an5, fl_unique_an5, 'as_given', eps_curv, eps_area, eps_min_ang);

% """ faces should be big in both states to be tracked """
% mask_tracked_an4 = int32(mask_tracked_an4);
face_area_an4_all = face_itg_curv_an4(mask_tracked_an4,1);
face_area_an5_all = face_itg_curv_an5(mask_tracked_an5,1);
mask_good_track = face_area_an4_all > 20 & face_area_an5_all > 20;
mask_good_an4 = face_itg_curv_an4(:,1) > 20;
mask_good_an5 = face_itg_curv_an5(:,1) > 20;

fl_unique_an4 = fl_unique_an4(mask_good_an4, :);
fl_unique_an5 = fl_unique_an5(mask_good_an5, :);
face_itg_curv_an4 = face_itg_curv_an4(mask_good_an4, :);
face_itg_curv_an5 = face_itg_curv_an5(mask_good_an5, :);

mask_tracked_an4(mask_tracked_an4) = mask_good_track;
mask_tracked_an5(mask_tracked_an5) = mask_good_track;
mask_tracked_an4 = mask_tracked_an4(mask_good_an4, :);
mask_tracked_an5 = mask_tracked_an5(mask_good_an5, :);

% ########################### Determine Face Type ###########################
% ----- first determine surface / inner -----
% """
% type = 0: inside volume, at most one of the two grains are true (not D3D) free surface. 
% type = 1: on the true (not D3D) free surface. 
% type = 2: connected to the true (not D3D) free surface. 
% type = 3: the rest faces, namely connected to D3D, or FOV, free surface
% 0-3: tracked, 10-13: not tracked
% """
face_type_an4 = ones(size(fl_unique_an4, 1), 1) * 13;
face_type_an5 = ones(size(fl_unique_an5, 1), 1) * 13;
num_surf_grain_an4_true = sum(ismember(fl_unique_an4, surf_grain_an4_true), 2);
num_surf_grain_an5_true = sum(ismember(fl_unique_an5, surf_grain_an5_true), 2);
face_type_an4(num_surf_grain_an4_true == 0) = 10;
face_type_an5(num_surf_grain_an5_true == 0) = 10;
face_type_an4(any(fl_unique_an4==0, 2)) = 11;
face_type_an5(any(fl_unique_an5==0, 2)) = 11;
face_type_an4(num_surf_grain_an4_true == 2) = 12;
face_type_an5(num_surf_grain_an5_true == 2) = 12;

% ----- then determine tracked / not tracked (disappeared | appeared) -----
face_type_an4(mask_tracked_an4) = face_type_an4(mask_tracked_an4) - 10;
face_type_an5(mask_tracked_an5) = face_type_an5(mask_tracked_an5) - 10;


% % ----- display stats -----
% tmp = setdiff(fl_unique_an4(face_type_an4 == 0, :), tracked_uniqueface_an4_complete, 'rows');
% sum(ismember(tmp, tracked_uniqueface_an4_all, 'rows'))
% ismember(tmp, surf_grain_an4_true);


% ########################### Calculate Face Properties ###########################
% ----- centroid -----
face_centroids_an4 = calcFaceCentroid(file_an4, fl_unique_an4);
face_centroids_an5 = calcFaceCentroid(file_an5, fl_unique_an5);

% ----- store data -----
centroids_disappeard = face_centroids_an4(~ mask_tracked_an4, :);
face_type_disappeared = face_type_an4(~ mask_tracked_an4, :);
area_disappeared = face_itg_curv_an4(~ mask_tracked_an4, 1);
face_itg_curv_disappeared = face_itg_curv_an4(~ mask_tracked_an4, 1);

centroids_appeard = face_centroids_an5(~ mask_tracked_an5, :);
face_type_appeared = face_type_an5(~ mask_tracked_an5, :);
area_appeared = face_itg_curv_an4(~ mask_tracked_an5, 1);
face_itg_curv_appeared = face_itg_curv_an4(~ mask_tracked_an5, 2);

centroids_tracked = face_centroids_an4(mask_tracked_an4, :);
face_type_tracked = face_type_an4(mask_tracked_an4, :);
area_tracked_an4 = face_itg_curv_an4(mask_tracked_an4, 1);
face_itg_curv_tracked_an4 = face_itg_curv_an4(mask_tracked_an4, 2);
area_diff = face_itg_curv_an5(mask_tracked_an5, 1) - face_itg_curv_an4(mask_tracked_an4, 1);
face_itg_curv_diff = face_itg_curv_an5(mask_tracked_an5, 2) - face_itg_curv_an4(mask_tracked_an4, 2);


% ##################################### Grains Data #####################################
surf_grain_an4_d3d = h5read(file_an4, '/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures')';
surf_grain_an5_d3d = h5read(file_an5, '/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures')';
grain_centroids_an4 = h5read(file_an4, '/DataContainers/ImageDataContainer/CellFeatureData/Centroids')';
grain_centroids_an5 = h5read(file_an5, '/DataContainers/ImageDataContainer/CellFeatureData/Centroids')';
sizes_an4 = h5read(file_an4, '/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters')';
sizes_an5 = h5read(file_an5, '/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters')';
surf_grain_an4_d3d(1) = [];
surf_grain_an5_d3d(1) = [];
grain_centroids_an4(1, :) = [];
grain_centroids_an5(1, :) = [];
sizes_an4(1) = [];
sizes_an5(1) = [];

% ----- Get surface grains by D3D -----
tmp = (1:size(surf_grain_an4_d3d, 1))';
surf_grain_an4_d3d = tmp(boolean(surf_grain_an4_d3d));
tmp = (1:size(surf_grain_an5_d3d, 1))';
surf_grain_an5_d3d = tmp(boolean(surf_grain_an5_d3d));

% ----- disappeared & appeared grains -----
ids_an4 = (1:size(sizes_an4, 1))';
ids_an5 = (1:size(sizes_an5, 1))';
mask_disappear = ~ ismember(ids_an4, look_up_table(:,1));
mask_appear = ~ ismember(ids_an5, look_up_table(:,2));

% ----- store data -----
size_tracked_an4 = sizes_an4(look_up_table(:,1));
size_diff = sizes_an5(look_up_table(:,2)) - sizes_an4(look_up_table(:,1));
centroids_trakced = grain_centroids_an4(look_up_table(:,1), :);

size_disappear = sizes_an4(mask_disappear);
centroids_disappear = grain_centroids_an4(mask_disappear, :);
size_appear = sizes_an5(mask_appear);
centroids_appear = grain_centroids_an5(mask_appear, :);

save(fname, 'centroids_appeard', 'face_type_appeared', 'area_appeared', 'face_itg_curv_appeared', ...
    'centroids_disappeard', 'face_type_disappeared', 'area_disappeared', 'face_itg_curv_disappeared',...
    'centroids_tracked', 'face_type_tracked', 'area_tracked_an4', 'face_itg_curv_tracked_an4', ...
    'area_diff', 'face_itg_curv_diff', ...
    'size_tracked_an4', 'size_diff', 'centroids_trakced', ...
    'size_disappear', 'centroids_disappear', 'size_appear', 'centroids_appear');






