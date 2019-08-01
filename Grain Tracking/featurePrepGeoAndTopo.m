% ##########################################################################
% * Notes
%     - This script is to prepare the geometric and topological features. 
%         1. Geometry Feature
%         2. Topology Feature
%         3. Dihedral Angles
%         4. Write txt File
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
[faces_an4, faces_an5, face_corresp] = trackFace(file_an4, file_an5, look_up_table, 'use_complete_faces');
[tracked_uniqueface_an4, tracked_uniqueface_an5] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_complete_faces');

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
tracked_uniqueface_an4_sorted = zeros(size(tracked_uniqueface_an5));
for i = 1:size(tracked_uniqueface_an5, 1)
    for j = 1:size(tracked_uniqueface_an5, 2)
        tracked_uniqueface_an4_sorted(i, j) = look_up_table_2to1(tracked_uniqueface_an5(i,j), 1);
    end
end

%% #################################### 1. Geometry Feature  #################################### 
face_tmp_an4 = calcFaceItgCurv(file_an4, tracked_uniqueface_an4, 'as_given', eps_curv, eps_area, eps_min_ang);
face_tmp_an5 = calcFaceItgCurv(file_an5, tracked_uniqueface_an5, 'as_given', eps_curv, eps_area, eps_min_ang);

face_area_an4 = face_tmp_an4(:,1);
face_itg_abscurv_an4 = face_tmp_an4(:,2);
face_avg_abscurv_an4 = face_itg_abscurv_an4 ./ face_area_an4;

diff_tmp = face_tmp_an5 - face_tmp_an4;
face_area_diff = diff_tmp(:,1);
face_itg_abscurv_diff = diff_tmp(:,2);
face_avg_abscurv_diff = face_itg_abscurv_diff ./ face_area_diff;

%%  #################################### 2. Topology Feature #################################### 
load('../190621_Ni_crop_TopologyResult_uniqiue.mat', 'faces_an4', 'faces_an5', ...
    'num_corners_an4', 'num_corners_an5', 'num_nnface_avgcorner_an4', ...
    'num_nnface_avgcorner_an5', 'avg_nng_diff', 'max_nng_diff');

% ----- Get the faces to calculate -----
mask_uniqueface_an4 = ismember(faces_an4, tracked_uniqueface_an4, 'rows');
mask_uniqueface_an5 = ismember(faces_an5, tracked_uniqueface_an5, 'rows');

% ----- Calculate num_corners and corner_diff -----
num_corners_an4 = num_corners_an4(mask_uniqueface_an4, :);
num_corners_an5 = num_corners_an5(mask_uniqueface_an5, :);
corner_diff = num_corners_an5 - num_corners_an4;



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
% [da_len_w_an4, da_num_w_an4] = calcGrainFaceDAs(tracked_uniqueface_an4_sorted, triple_line_an4, tl_info_an4);
% [da_len_w_an5, da_num_w_an5] = calcGrainFaceDAs(tracked_uniqueface_an5, triple_line_an5, tl_info_an5);
% 
% da_len_w_diff = da_len_w_an5 - da_len_w_an4;
% da_num_w_diff = da_num_w_an5 - da_num_w_an4;


%%  #################################### 4. Write txt File #################################### 
area_thres = 20;
ratio_thres_pos = 10;
ratio_thres_neg = 1/ratio_thres_pos - 1;
fa_diff_ratio = face_area_diff ./ face_area_an4;
mask_good_face = (face_area_an4 > area_thres & face_area_an4 + face_area_diff > area_thres ...
    & fa_diff_ratio < ratio_thres_pos & fa_diff_ratio > ratio_thres_neg);
fileID = fopen('190718_Hsmooth_mask_good_face.txt','w');
fprintf(fileID,'%s\n', 'mask_good_face');
for i = 1:size(mask_good_face, 1)
    fprintf(fileID, '%3d\n', mask_good_face(i));
end
fclose(fileID);

fileID = fopen('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_Lsmooth_geo_topo.txt','w');
fprintf(fileID,'%s, %s, %s, %s, %s, %s, %s, %s\n', ...
    'A_an4', 'fMs_an4', 'avg_FabsavgH_an4', 'C_an4', ...
    'A_diff', 'fMs_diff', 'avg_FabsavgH_diff', 'C_diff');
for i = 1:size(face_area_an4, 1)
    fprintf(fileID, '%6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f\n', ...
        face_area_an4(i), face_itg_abscurv_an4(i), face_avg_abscurv_an4(i), num_corners_an4(i),...
        face_area_diff(i), face_itg_abscurv_diff(i), face_avg_abscurv_diff(i), corner_diff(i));
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


