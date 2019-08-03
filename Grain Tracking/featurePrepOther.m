% ##########################################################################
% * Notes
%     - Misorientaion distance is calculated from misorientation angle. 
%     - Distance to coherent twin is calculated as defined by Moraweic.
%         1. Identify (inner) complete faces 
%         2. Identify broken faces
%         3.1. Identify twins, from D3D
%         3.2. Myself distance from twins
%         4. Identify distance from target normals
%         5. Write txt file
% ##########################################################################

% load('/Users/xiaotingzhong/Documents/Matlab/Grain Tracking/data/181108_rfvec.mat', 'rfvecs_an4', 'rfvecs_an5');
% file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin3_Hsmooth.dream3d';
% file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4_Hsmooth.dream3d';
file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4.dream3d';
load('look_up_table_an4_an5crop.mat')
eps_curv = 1;
eps_area = 7;
eps_min_ang = 10;

% ------------- get faces to track -------------
% """ Calculate all inner faces and record which ones are complete """
[tracked_uniqueface_an4_inner, tracked_uniqueface_an5_inner] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_inner_faces');
obj_faces_an4 = tracked_uniqueface_an4_inner;
obj_faces_an5 = tracked_uniqueface_an5_inner;


%% #################################### 1. Identify (inner) complete faces ####################################
[tracked_uniqueface_an4_complete, tracked_uniqueface_an5_complete] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_complete_faces');
complete = ismember(sort(tracked_uniqueface_an4_inner, 2), sort(tracked_uniqueface_an4_complete, 2), 'rows');



%% #################################### 2. Identify broken faces ####################################
is_one_piece = checkIfOnepiece(file_an4, obj_faces_an4)
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_one_piece.mat', 'is_one_piece');



%% #################################### 3.1. Identify twins, from D3D ####################################
[is_twin_an4, incoherence_an4] = calcTwinsFromD3D(file_an4, obj_faces_an4, eps_curv, eps_area, eps_min_ang);
[is_twin_an5, incoherence_an5] = calcTwinsFromD3D(file_an5, obj_faces_an5, eps_curv, eps_area, eps_min_ang);



%% #################################### 3.2. Myself distance from twins ####################################
% % ----- BAD: by Rodrigues vector -----
% % rfvec_twin = [1,1,1]/norm([1,1,1]) * tand(60/2);
% % % """ Note 'AvgEAs' | 'AvgEulerAngles' """
% % [rfvecs_an4]  = getFaceRFvecs(file_an4, obj_faces_an4);
% % [rfvecs_an5]  = getFaceRFvecs(file_an5, obj_faces_an5);
% % dist_twin_an4 = vecnorm(rfvecs_an4 - rfvec_twin, 2, 2);
% % dist_twin_an5 = vecnorm(rfvecs_an5 - rfvec_twin, 2, 2);

% % ----- By misorientation matrix -----
% dg_obj = zeros(3, 3, 6);
% dg_obj(:,:,1) = AAToG(60, [1, 1, 1]);     % sigma3
% dg_obj(:,:,2) = AAToG(36.86, [1, 0, 0]);  % sigma5
% dg_obj(:,:,3) = AAToG(38.21, [1, 1, 1]);  % sigma7
% dg_obj(:,:,4) = AAToG(38.94, [1, 1, 0]);  % sigma9
% dg_obj(:,:,5) = AAToG(31.59, [1, 1, 0]);  % sigma27a
% dg_obj(:,:,6) = AAToG(35.43, [2, 1, 0]);  % sigma27b
% 
dists_miso = calcDistFromMisorientation(file_an4, obj_faces_an4, dg_obj);

load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_Lsmoth_miso_dists.mat')

%% #################################### 4. Identify distance from target normals ####################################
% target_normals = [1,1,1;1,1,0;1,0,0;1,1,2];
% target_normals = target_normals ./ vecnorm(target_normals, 2, 2);
[dists_norm_an4, std_norm_an4] = calcDistFromTargetNormals(file_an4, obj_faces_an4, target_normals, eps_curv, eps_area, eps_min_ang);
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_normal_dists.mat', 'face_dist_avgs', 'face_dist_stds');



%% #################################### 5. Write txt file ####################################
not_twin_an4 = ~is_twin_an4;
dist_twin_an4 = dists_miso(:,1);

fileID = fopen('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_Lsmooth_otherinfo.txt','w');
fprintf(fileID,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n', ...
    'not_twin_an4', 'dist_twin_an4', 'incoherence_an4', ...
    'dist_avg_111', 'dist_std_111', 'dist_avg_110', 'dist_std_110', ...
    'dist_avg_100', 'dist_std_100', 'dist_avg_211', 'dist_std_211', ...
    'mask_complete', 'mask_onepiece');
for i = 1:length(not_twin_an4)
    fprintf(fileID, '%6d, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6d, %6d\n', ...
        not_twin_an4(i), dist_twin_an4(i), incoherence_an4(i), ...
        face_dist_avgs(i, 1), face_dist_stds(i, 1), face_dist_avgs(i, 2), face_dist_stds(i, 2), ...
        face_dist_avgs(i, 3), face_dist_stds(i, 3), face_dist_avgs(i, 4), face_dist_stds(i, 4), ...
        complete(i), is_one_piece(i)) ;
end
fclose(fileID);

