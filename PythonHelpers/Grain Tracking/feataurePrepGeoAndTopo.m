% ##########################################################################
% * Notes
%     - This script is to prepare the geometric and topological features. 
% ##########################################################################
file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin2_smooth.dream3d');
file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
load('look_up_table_an4_an5.mat')

% ##### get the faceLabels and their correpondence ##### 
[faces_an4, faces_an5, face_corresp] = trackFace(file_an4, file_an5, look_up_table, 'use_complete_faces');
[tracked_uniqueface_an4, tracked_uniqueface_an5] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_complete_faces');

%% #################################### Geometry Feature  #################################### 
face_tmp_an4 = calcFaceItgCurv(file_an4, tracked_uniqueface_an4, 'unique_faces');
face_tmp_an5 = calcFaceItgCurv(file_an5, tracked_uniqueface_an5, 'unique_faces');

face_area_an4 = face_tmp_an4(:,1);
face_itg_curv_an4 = face_tmp_an4(:,2);
face_avg_curv_an4 = face_itg_curv_an4 ./ face_area_an4;

diff_tmp = face_tmp_an5 - face_tmp_an4;
face_area_diff = diff_tmp(:,1);
face_itg_curv_diff = diff_tmp(:,2);
face_avg_curv_diff = face_itg_curv_diff ./ face_area_diff;


%%  #################################### Topology Feature #################################### 
load('../Topologies/181101_Ni_TopologyResult_uniqiue.mat', 'faces_an4', 'faces_an5', ...
    'num_corners_an4', 'num_corners_an5', 'num_nnface_avgcorner_an4', 'num_nnface_avgcorner_an5');

% ----- Get the faces to calculate -----
mask_uniqueface_an4 = ismember(faces_an4, tracked_uniqueface_an4, 'rows');
mask_uniqueface_an5 = ismember(faces_an5, tracked_uniqueface_an5, 'rows');

% ----- Calculate num_corners and corner_diff -----
num_corners_an4 = num_corners_an4(mask_uniqueface_an4, :);
num_corners_an5 = num_corners_an5(mask_uniqueface_an5, :);
corner_diff = num_corners_an5 - num_corners_an4;

% ----- Calculate num_corners - avg(num_corners_nearest_neigh_face) -----
num_nnface_avgcorner_an4 = num_nnface_avgcorner_an4(mask_uniqueface_an4);
num_nnface_avgcorner_an5 = num_nnface_avgcorner_an5(mask_uniqueface_an5);
c_nnc_an4 = num_corners_an4 - num_nnface_avgcorner_an4;
c_nnc_an5 = num_corners_an5 - num_nnface_avgcorner_an5;

c_nnc_diff = c_nnc_an5 - c_nnc_an4;


%%  #################################### Write txt File #################################### 
fileID = fopen('190408_features_geo_topo.txt','w');
fprintf(fileID,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n', ...
    'A_an4', 'fMs_an4', 'avg_FabsavgH_an4', 'C_an4', 'CnnC_an4', ...
    'A_diff', 'fMs_diff', 'avg_FabsavgH_diff', 'C_diff', 'CnnC_diff');
for i = 1:length(face_area_an4)
    fprintf(fileID, '%6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f\n', ...
        face_area_an4(i), face_itg_curv_an4(i), face_avg_curv_an4(i), num_corners_an4(i), c_nnc_an4(i),...
        face_area_diff(i), face_itg_curv_diff(i), face_avg_curv_diff(i), corner_diff(i), c_nnc_diff(i));
end
fclose(fileID);
