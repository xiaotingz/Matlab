% ##########################################################################
% * Notes
%     - This script is to prepare the geometric and topological features. 
%         1. Geometry Feature
%         2. Topology Feature
%         3. Dihedral Angles
%         4. Write txt File
%         5.1. Check geometry data
%         5.2. Check topology data
%     - Preparation of the topology features: see main.m in /Topologies
% ##########################################################################
% file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
% file_an5 = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4_Hsmooth.dream3d';
file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4.dream3d';
load('look_up_table_an4_an5crop.mat')
% load('/Volumes/XIAOTING/Ni/working/190621_tracked_faces_full.mat')
% tracked_uniqueface_an4 = tracked_uniqueface_an4_full;
% tracked_uniqueface_an5 = tracked_uniqueface_an5_full;

% """
% Threshold values for determining if a triangle is good. Use in calcFaceToGrainCentroidDist.m
% """
eps_curv = 1;
eps_area = 7;
eps_min_ang = 10;

% ##### get the faceLabels and their correpondence ##### 
% [faces_an4, faces_an5, face_corresp] = trackFace(file_an4, file_an5, look_up_table, 'use_complete_faces');
[tracked_uniqueface_an4_inner, tracked_uniqueface_an5_inner] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_inner_faces');
obj_faces_an4 = tracked_uniqueface_an4_inner;
obj_faces_an5 = tracked_uniqueface_an5_inner;

% ##### Make label_an4 & label_an5 one-to-one for the tracked faces #####
% """
% Recalc tracked_uniqueface_an4 To Keep Grains Consistent With track_uniquefaces_an5 
% """
% ----- make a real index table out of look_up_table -----
look_up_table = sortrows(look_up_table, 2);
look_up_table_2to1 = zeros(max(look_up_table(:,2)),2);
idx = 1;
for i = 1 : max(look_up_table(:,2))
    if look_up_table(idx,2) == i
        look_up_table_2to1(i,1) = look_up_table(idx,1);
        look_up_table_2to1(i,2) = look_up_table(idx,2);
        idx = idx + 1;
    else
        look_up_table_2to1(i,1) = NaN;
        look_up_table_2to1(i,2) = i;
    end
end
% ----- remake tracked_uniqueface_an4 -----
obj_faces_an4_sorted = zeros(size(obj_faces_an5));
for i = 1:size(obj_faces_an5, 1)
    for j = 1:size(obj_faces_an5, 2)
        obj_faces_an4_sorted(i, j) = look_up_table_2to1(obj_faces_an5(i,j), 1);
    end
end

%% #################################### 1. Geometry Feature  #################################### 
% ----- face_itg_abscurv -----
face_tmp_an4 = calcFaceItgCurv(file_an4, obj_faces_an4, 'abs', eps_curv, eps_area, eps_min_ang);
face_tmp_an5 = calcFaceItgCurv(file_an5, obj_faces_an5, 'abs', eps_curv, eps_area, eps_min_ang);

face_area_an4 = face_tmp_an4(:,1);
face_itg_abscurv_an4 = face_tmp_an4(:,2);
face_avg_abscurv_an4 = face_itg_abscurv_an4 ./ face_area_an4;

diff_tmp = face_tmp_an5 - face_tmp_an4;
face_area_diff = diff_tmp(:,1);
face_itg_abscurv_diff = diff_tmp(:,2);
face_avg_abscurv_diff = face_itg_abscurv_diff ./ face_area_diff;

% ----- Signed face_itg_curv, using left grain as resident -----
face_tmp_an4 = calcFaceItgCurv(file_an4, obj_faces_an4_sorted, 'signed_resident_left', eps_curv, eps_area, eps_min_ang);
face_tmp_an5 = calcFaceItgCurv(file_an5, obj_faces_an5, 'signed_resident_left', eps_curv, eps_area, eps_min_ang);

face_itgcurv_an4_left = face_tmp_an4(:,2);
face_itgcurv_an5_left = face_tmp_an4(:,2);

face_itgcurv_left_diff = face_itgcurv_an5_left - face_itgcurv_an4_left;

%%  #################################### 2. Topology Feature #################################### 
% """ go to ../Topologies/main to do calculation """
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_Lsmoth_Topologies.mat', ...
    'faces_an4', 'faces_an5', 'num_corners_an4', 'num_corners_an5', 'num_edges_an4', 'num_edges_an5');

% ----- Prepare data dictionary for faces -----
dict_corner_an4 = containers.Map('KeyType','char','ValueType','double');
dict_edge_an4 = containers.Map('KeyType','char','ValueType','double');
for i = 1:size(faces_an4, 1)
    key = mat2str(faces_an4(i, :));
    dict_corner_an4(key) = num_corners_an4(i);
    dict_edge_an4(key) = num_edges_an4(i);
end

dict_corner_an5 = containers.Map('KeyType','char','ValueType','double');
dict_edge_an5 = containers.Map('KeyType','char','ValueType','double');
for i = 1:size(faces_an5, 1)
    key = mat2str(faces_an5(i, :));
    dict_corner_an5(key) = num_corners_an5(i);
    dict_edge_an5(key) = num_edges_an5(i);
end

% ----- Get #corners and #edges for the obj_faces -----
face_corners_an4 = zeros(size(obj_faces_an4, 1), 1);
face_corners_an5 = zeros(size(obj_faces_an5, 1), 1);
face_edges_an4 = zeros(size(obj_faces_an4, 1), 1);
face_edges_an5 = zeros(size(obj_faces_an5, 1), 1);
for i = 1:size(obj_faces_an4, 1)
    key = mat2str(obj_faces_an4(i, :));
    face_corners_an4(i) = dict_corner_an4(key);
    face_edges_an4(i) = dict_edge_an4(key);

    key = mat2str(obj_faces_an5(i, :));
    face_corners_an5(i) = dict_corner_an5(key);
    face_edges_an5(i) = dict_edge_an5(key);
end
corner_diff = face_corners_an5 - face_corners_an4;
edge_diff = face_edges_an5 - face_edges_an4;


%%  #################################### 3. Dihedral Angles #################################### 
% % """
% % - Go to /Topologies. 
% % - Note use tracked_uniqueface_sorted 
% % """
% eps_area = 7;
% eps_curv = 1;
% eps_min_ang = 10;
% % [triple_line_an4, tl_info_an4] = findTripleLines(file_an4, eps_area, eps_curv, eps_min_ang);
% [triple_line_an5, tl_info_an5] = findTripleLines(file_an5, eps_area, eps_curv, eps_min_ang);
% [da_len_w_an4, da_num_w_an4] = calcGrainFaceDAs(obj_faces_an4_sorted, triple_line_an4, tl_info_an4);
% [da_len_w_an5, da_num_w_an5] = calcGrainFaceDAs(obj_faces_an5, triple_line_an5, tl_info_an5);
% 
% da_len_w_diff = da_len_w_an5 - da_len_w_an4;
% da_num_w_diff = da_num_w_an5 - da_num_w_an4;


%%  #################################### 4. Write txt File #################################### 
% area_thres = 20;
% ratio_thres_pos = 10;
% ratio_thres_neg = 1/ratio_thres_pos - 1;
% fa_diff_ratio = face_area_diff ./ face_area_an4;
% mask_good_face = (face_area_an4 > area_thres & face_area_an4 + face_area_diff > area_thres ...
%     & fa_diff_ratio < ratio_thres_pos & fa_diff_ratio > ratio_thres_neg);
% fileID = fopen('190718_Hsmooth_mask_good_face.txt','w');
% fprintf(fileID,'%s\n', 'mask_good_face');
% for i = 1:size(mask_good_face, 1)
%     fprintf(fileID, '%3d\n', mask_good_face(i));
% end
% fclose(fileID);

fileID = fopen('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_Lsmooth_geo_topo.txt','w');
fprintf(fileID,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n', ...
    'A_an4', 'fMs_abs_an4', 'fMs_signed_an4', 'avg_FabsavgH_an4', 'C_an4', ...
    'A_diff', 'fMs_abs_diff', 'fMs_signed_diff', 'avg_FabsavgH_diff', 'C_diff');
for i = 1:size(face_area_an4, 1)
    fprintf(fileID, '%6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f\n', ...
        face_area_an4(i), face_itg_abscurv_an4(i), face_itgcurv_an4_left(i), face_avg_abscurv_an4(i), face_corners_an4(i),...
        face_area_diff(i), face_itg_abscurv_diff(i), face_itgcurv_left_diff(i), face_avg_abscurv_diff(i), corner_diff(i));
end
fclose(fileID);

% ----- Write DA data -----
% fileID = fopen('190425_features_dihedral_ang.txt','w');
% fprintf(fileID,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n', ...
%     'da_len_w_an4_l', 'da_len_w_an4_r', 'da_len_w_an4_opp', 'da_num_w_an4_l', 'da_num_w_an4_r', 'da_num_w_an4_opp', ...
%     'da_len_w_diff_l', 'da_len_w_diff_r', 'da_len_w_diff_opp', 'da_num_w_diff_l', 'da_num_w_diff_r', 'da_num_w_diff_opp');
% for i = 1:length(face_area_an4)
%     fprintf(fileID, '%6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f\n', ...
%         da_len_w_an4(i, 1), da_len_w_an4(i, 2), da_len_w_an4(i, 3), da_num_w_an4(i, 1), da_num_w_an4(i, 2), da_num_w_an4(i, 3), ...
%         da_len_w_diff(i, 1), da_len_w_diff(i, 2), da_len_w_diff(i, 3), da_num_w_diff(i, 1), da_num_w_diff(i, 2), da_num_w_diff(i, 3));
% end
% fclose(fileID);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 5. CHECKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 5.1. Check geometry data
% 5.2. Check topology data

% ----- Check face_corresp -----
% """ 
% 1. Use look_up_tabke to convert between correpondences, see if match.
% 2. Plot the correspondence in paraview. This way area diff and itgcurv_diff can also be checked. 
% """
% ----- Prepare data -----
featurefacelabel_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
featurefacelabel_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
featurefacelabel_an4(1, :) = [];
featurefacelabel_an5(1, :) = [];
face_idx_an4 = (1:length(featurefacelabel_an4))';
face_idx_an5 = (1:length(featurefacelabel_an5))';

tri_fl = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
tri_curv =  roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)';
tri_area = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
tri_min_ang = roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
mask = all(tri_fl > 0, 2) & abs(tri_curv) < eps_curv & tri_area < eps_area & tri_min_ang > eps_min_ang;
tri_fl = tri_fl(mask, :);
tri_curv = tri_curv(mask);
tri_area = tri_area(mask);
tri_min_ang = tri_min_ang(mask);

num_tris_an4 = calcNumTrianglesOnFace(file_an4, obj_faces_an4, eps_curv, eps_area, eps_min_ang);


%% ######################################## 5.1. Check geometry data ########################################

% ----- Show data from above calculation -----
rng('shuffle');
% idx = randi(size(obj_faces_an4, 1));
% 7397 1695 9193 4088
idx = 4088;

label_an4 = obj_faces_an4_sorted(idx, :);
label_an5 = obj_faces_an5(idx, :);

mask_an4 = ismember(featurefacelabel_an4, sort(label_an4), 'rows');
mask_an5 = ismember(featurefacelabel_an5, sort(label_an5), 'rows');
featurefaceid_an4 = face_idx_an4(mask_an4);
featurefaceid_an5 = face_idx_an5(mask_an5);

disp('----------------------------------')
disp(['idx_corresp = ', num2str(idx)])
disp('-----')
disp('an4')
disp(['label = ', mat2str(label_an4), ',   face_feature_id = ', num2str(featurefaceid_an4)])
disp(['area_an4 = ', num2str(face_tmp_an4(idx, 1)), ...
        ',   itgcurv_abs_an4 = ', num2str(face_tmp_an4(idx, 2)), ...
        ',   itgcurv_left_an4 = ', num2str(face_itgcurv_an4_left(idx)), ...
        ',   num_tris_an4 = ', num2str(num_tris_an4(idx))]);
disp('an5')
disp(['label = ', mat2str(label_an5), ',   face_feature_id = ', num2str(featurefaceid_an5)])
disp(['area_an5 = ', num2str(face_tmp_an5(idx, 1)), ...
        ',   itgcurv_abs_an5 = ', num2str(face_tmp_an5(idx, 2)), ...
        ',   itgcurv_left_an5 = ', num2str(face_itgcurv_an5_left(idx))]);
disp('-----')
disp(['area_diff = ', num2str(face_area_diff(idx)), ',   itgcurv_diff = ', num2str(face_itg_abscurv_diff(idx))])



% ----- Calculate again from raw data -----
mask_1 = (tri_fl(:,1)==label_an4(1) & tri_fl(:,2)==label_an4(2));
mask_2 = (tri_fl(:,1)==label_an4(2) & tri_fl(:,2)==label_an4(1));
mask = mask_1 | mask_2;
itg_abscurv = sum(tri_area(mask) .* abs(tri_curv(mask)));
itg_curv_left = - sum(tri_area(mask_1) .* tri_curv(mask_1)) + sum(tri_area(mask_2) .* tri_curv(mask_2));
avg_curv = sum(tri_curv(mask)) / sum(mask);
avg_abscurv = sum(abs(tri_curv(mask))) / sum(mask);
area = sum(tri_area(mask));


disp('-----')
disp('check from raw')
disp(['area_an4 = ', num2str(area),  ',   itgcurv_abs_an4 = ', num2str(itg_abscurv), ...
        ',   itgcurv_left_an4 = ', num2str(itg_curv_left)]);
disp(['avgcurv_abs_an4 = ', num2str(avg_abscurv), ',   num_tris_on_face = ', num2str(sum(mask))]);
disp(' ')

    
%% ######################################## 5.2. Check topology data ########################################
% see /Topologies/main, section checks












