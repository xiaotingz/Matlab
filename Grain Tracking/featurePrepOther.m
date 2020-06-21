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
%         6. rv_vecs & PCA rv_vecs
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
is_one_piece = checkIfOnepiece(file_an4, obj_faces_an4);
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
dg_obj = zeros(3, 3, 6);
dg_obj(:,:,1) = AAToG(60, [1, 1, 1]);     % sigma3
dg_obj(:,:,2) = AAToG(36.86, [1, 0, 0]);  % sigma5
dg_obj(:,:,3) = AAToG(38.21, [1, 1, 1]);  % sigma7
dg_obj(:,:,4) = AAToG(38.94, [1, 1, 0]);  % sigma9
dg_obj(:,:,5) = AAToG(31.59, [1, 1, 0]);  % sigma27a
dg_obj(:,:,6) = AAToG(35.43, [2, 1, 0]);  % sigma27b



% dists_miso = calcDistFromMisorientation(file_an4, obj_faces_an4, dg_obj);
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_Lsmoth_miso_dists.mat')

%% #################################### 4. Identify distance from target normals ####################################
% target_normals = [1,1,1;1,1,0;1,0,0;1,1,2];
% target_normals = target_normals ./ vecnorm(target_normals, 2, 2);
% [dists_norm_an4, std_norm_an4, tri_fls, tri_dists] = calcDistFromTargetNormals(file_an4, obj_faces_an4, target_normals, eps_curv, eps_area, eps_min_ang);
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_normal_dists.mat', 'face_dist_avgs', 'face_dist_stds');



%% #################################### 5. Write txt file ####################################
% not_twin_an4 = ~is_twin_an4;
dist_twin_an4 = dists_miso(:,1);
is_one_piece = is_one_piece_an4 == 1 & is_one_piece_an5 == 1;

fileID = fopen('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_Lsmooth_otherinfo_normdist_formal.txt','w');
fprintf(fileID,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n', ...
    'mask_complete', 'mask_onepiece', 'dist_twin_an4', ...
    'dist_avg_111', 'dist_std_111', 'dist_avg_110', 'dist_std_110', ...
    'dist_avg_100', 'dist_std_100', 'dist_avg_211', 'dist_std_211');
for i = 1:length(not_twin_an4)
    fprintf(fileID, '%6d, %6d, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f\n', ...
        complete(i), is_one_piece(i), dist_twin_an4(i),  ...
        face_dist_avgs(i, 1), face_dist_stds(i, 1), face_dist_avgs(i, 2), face_dist_stds(i, 2), ...
        face_dist_avgs(i, 3), face_dist_stds(i, 3), face_dist_avgs(i, 4), face_dist_stds(i, 4)) ;
end
fclose(fileID);


%% #################################### 5. rv_vecs & PCA rv_vecs ####################################
% """
% Ref: 
%     - 1) http://www.cs.cmu.edu/~guestrin/Class/10701-S05/slides/dimensionality.pdf
%     - 2) https://stats.stackexchange.com/questions/229092/how-to-reverse-pca-and-reconstruct-original-variables-from-several-principal-com
%     - 3) https://www.jeremyjordan.me/principal-components-analysis/
%     - 4) CS229 Lecture notes, Andrew Ng, PCA
% Notes:
%     - PCA is to compress the number of data points (n), but what we want is to compress the data dimensions (d).
%     - To compress the dimensions, first find the principle component then project along that direction.
% """

rf_miso_an4 = calcMisorientationAsRFvec(obj_faces_an4, file_an4);
% rf_miso_an5 = calcMisorientationAsRFvec(obj_faces_an5, file_an5);


% ------- Normalization, necessary -------
% """ 
% Note PCA normalization is within each column or feature. (Ref 4)
% std(data, 1) uses N instead of N-1
% ""
rf_miso_an4_mean = mean(rf_miso_an4);
rf_miso_an4_norm = rf_miso_an4 - rf_miso_an4_mean;
rf_miso_std = std(rf_miso_an4_norm);
rf_miso_an4_norm = rf_miso_an4_norm ./ rf_miso_std;

% ------- PCA -------
% """ https://www.mathworks.com/help/stats/pca.html """
[rf_pc, rf_score, rf_latent] = pca(rf_miso_an4_norm);

% ------- Compute explained variance ------- 
rf_explained_var = rf_latent ./ sum(rf_latent);
disp(['explained variance = ', mat2str(round(rf_explained_var, 4))])

% ------- Projection ------- 
rf_proj = rf_miso_an4_norm * rf_pc;

% ------- double check with SVD -------
% % """ 
% % rf_pc = v should hold: 
% %     https://stats.stackexchange.com/questions/134282/relationship-between-svd-and-pca-how-to-use-svd-to-perform-pca
% % """
% % [u, s, v] = svd(rf_miso_an4_norm);
% % rf_proj_2 = u * s;
% % rf_proj_3 = rf_miso_an4_norm * v;


% ------- Restore matrix from projection -------
proj_mat = rf_pc(:,1) * rf_pc(:,1)';
rf_restore = (rf_miso_an4_norm * proj_mat).* rf_miso_std + rf_miso_an4_mean;
pca_score = 1 - sum(norm((rf_miso_an4 - rf_restore), 2)) / sum(norm(rf_miso_an4, 2));
disp(['pca_score = ', num2str(pca_score)])
% mae = sum(abs(rf_restore - rf_miso_an4_norm));
% disp(['MAE = ', mat2str(mae)])


fileID = fopen('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_Lsmooth_otherinfo_miso_rf.txt','w');
fprintf(fileID,'%s, %s, %s, %s\n', 'rf_1', 'rf_2', 'rf_3', 'rf_pc_1');
for i = 1:size(rf_miso_an4, 1)
    fprintf(fileID, '%6.3f, %6.3f, %6.3f, %6.3f\n', ...
        rf_miso_an4(i, 1), rf_miso_an4(i, 2), rf_miso_an4(i, 3), rf_proj(i, 1));
end
fclose(fileID);









