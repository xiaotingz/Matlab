% function mig_signs = calcFaceMigSign(file_an4, file_an5, tracked_uniqueface_an4, tracked_uniqueface_an5, ...
% data_grain_an4, data_grain_an5, data_face_an4, data_face_an5)
% ##########################################################################
% * Input
%     - tracked_uniqueface_ = [n,2]  
%           returned by trackFace.m or trackUniqueFace.m
%     - data_grain_ = [m, k]
%           m = #grains, k = #data_grain passed in
%           [size, F, or F - <Fnn>, integral curvature], which  actually includes grain property and grain neighborhood property
%     - data_face = [n, 4]
%           [label_A, label_B, area, itg_curv/area]
%           returned by calcFaceCurvature.m in /Grain Curvature
% * Output
%     - mig_signs = [n, k]
%           n = #grain_faces, k = #data_grain passed in
% * NOTE
%     - This function is to calculate the sign of face migration. 
%     - sign = 1    if face moving towards the grain with lower data_grain; 
%       sign = -1   if face moving toward the grain with higher data_grain;
%       sign = 0    if the two grains swapped data_grain level.
%     - Dependence: calcFaceToGrainCentroidDist.m

% """
% TODO: there are two ways to define data_grain relavance: just an4 only, or
% use data_grain difference. now is using an4 only. 
% """
% ##########################################################################
% ----------------------- load debug data -----------------------
% file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
% file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/181107.mat', 'tracked_uniqueface_an4', 'tracked_uniqueface_an5');
% load('look_up_table_an4_an5.mat')

% ##### run  /Grain Curvature/G_F_mF.m to get the following data: from data_grain & F_mF_diff #####
% data_grain_an4 = [grain_size_an4, grain_F_an4, grain_itg_curv_an4, F_avgFnn_an4];
% data_grain_an5 = [grain_size_an5, grain_F_an5, grain_itg_curv_an5, F_avgFnn_an5];
eps_motion = 1e-2;

% ---------------------------------------------------------------


%% #################################### Caculate sign, grain properties ####################################
% ##### Recalc tracked_uniqueface_an4 To Keep Grains Consistent With track_uniquefaces_an5 #####
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
tracked_uniqueface_an4 = zeros(size(tracked_uniqueface_an5));
for i = 1:size(tracked_uniqueface_an5, 1)
    for j = 1:size(tracked_uniqueface_an5, 2)
        tracked_uniqueface_an4(i, j) = look_up_table_2to1(tracked_uniqueface_an5(i,j), 1);
    end
end

% #################################### Check If Grain Face Moved Left ####################################
% ----- Absolute distance -----
% dists_f_g_an4 = calcFaceToGrainCentroidDist(file_an4, tracked_uniqueface_an4);
% dists_f_g_an5 = calcFaceToGrainCentroidDist(file_an5, tracked_uniqueface_an5);
% ----- Convert distance to portion -----
total_dists_an4 = sum(dists_f_g_an4, 2);
total_dists_an5 = sum(dists_f_g_an5, 2);
dists_f_g_an4 = dists_f_g_an4 ./ total_dists_an4;
dists_f_g_an5 = dists_f_g_an5 ./ total_dists_an5;

move_left = zeros(size(dists_f_g_an4,1), 1);
for i = 1:length(dists_f_g_an4)
    if dists_f_g_an4(i, 1) < dists_f_g_an5(i, 1) - eps_motion
        move_left(i) = -1;
    elseif dists_f_g_an4(i, 1) > dists_f_g_an5(i, 1) + eps_motion
        move_left(i) = 1;
    else
        move_left(i) = 0;
    end
end

%% #################################### Calc energy_gradient Defined By Grain Properties & Grain Neighborhood Properties ####################################
% """
% - eps = 0.05 * min(data_grain_left, data_grain_right) is a reasonable
%   choice for converting to binary values. 
%       However, probably should keep the values of grain_property_diff as
%       input of regression is probably better. 
% - diff = data_grain_an4_right - data_grain_an4_left becase move left <-> left grain small
% """
data_grain_an4_left = zeros(size(dists_f_g_an4,1), size(data_grain_an4, 2));
data_grain_an4_right = zeros(size(dists_f_g_an4,1), size(data_grain_an4, 2));
left_small = zeros(size(dists_f_g_an4,1), size(data_grain_an4, 2));
for j = 1:size(data_grain_an4_left, 2)
    for i = 1:size(data_grain_an4_left, 1)
        g_id_an4_left = tracked_uniqueface_an4(i, 1);
        g_id_an4_right = tracked_uniqueface_an4(i, 2);
        data_grain_an4_left(i, j) = data_grain_an4(g_id_an4_left, j);
        data_grain_an4_right(i, j) = data_grain_an4(g_id_an4_right, j);
        
        eps = min(data_grain_an4_left(i, j), data_grain_an4_right(i, j)) * 0.05;
        if data_grain_an4_left(i, j) < data_grain_an4_right(i, j) - eps
            left_small = 1;
        elseif data_grain_an4_left(i, j) > data_grain_an4_right(i, j) + eps
            left_small = -1;
        else
            left_small = 1;
        end
    end
end
data_grain_an4_diff = data_grain_an4_right - data_grain_an4_left;
% --- data_grain = [size, F, integral curvature] ---
data_grain_an4_diff(:,3) = - data_grain_an4_diff(:,3);

%% #################################### Calc energy_gradient Defined By Grain Face Properties ####################################
% """
% Hypothesis: if grain face convex when taking the left grain as resident grain, then grain face should move left.
% Positive cases: move left; left grain small; face area decrease
% """  
face_itgcurv_an4_left = zeros(size(tracked_uniqueface_an4, 1), 1);
face_area_an4_m_an5 = zeros(size(tracked_uniqueface_an4, 1), 1);
for i = 1:size(tracked_uniqueface_an4, 1)
    % ------------------ Grain Face Integral Curvature ------------------
    % ----- get itg_curv for the specific grain face -----
    mask_an4_1 = (tracked_uniqueface_an4(i, 1) == data_face_an4(:,1) & tracked_uniqueface_an4(i, 2) == data_face_an4(:,2));
    mask_an4_2 = (tracked_uniqueface_an4(i, 2) == data_face_an4(:,1) & tracked_uniqueface_an4(i, 1) == data_face_an4(:,2));
    if sum(mask_an4_1) > 0
        face_itg_curv_1 = data_face_an4(mask_an4_1, 3)*data_face_an4(mask_an4_1, 4);
        face_area_an4_1 = data_face_an4(mask_an4_1, 3);
    else
        face_itg_curv_1 = 0;
        face_area_an4_1 = 0;
    end
    if sum(mask_an4_2) > 0
        face_itg_curv_2 = data_face_an4(mask_an4_2, 3)*data_face_an4(mask_an4_2, 4);
        face_area_an4_2 = data_face_an4(mask_an4_2, 3);
    else
        face_itg_curv_2 = 0;
        face_area_an4_2 = 0;
    end
    % ----- adjust sign of curvature -----
    if tracked_uniqueface_an4(i, 1) == data_face_an4(mask_an4_1, 1)
        face_itgcurv_an4_left(i)  = - face_itg_curv_1 + face_itg_curv_2;
    else
        face_itgcurv_an4_left(i)  = face_itg_curv_1 - face_itg_curv_2;
    end
    
    % ------------------ Grain Face Area Diff ------------------
    % ----- get face_area -----
    face_area_an4 = face_area_an4_1 + face_area_an4_2;
    mask_an5_1 = (tracked_uniqueface_an5(i, 1) == data_face_an5(:,1) & tracked_uniqueface_an5(i, 2) == data_face_an5(:,2));
    mask_an5_2 = (tracked_uniqueface_an5(i, 2) == data_face_an5(:,1) & tracked_uniqueface_an5(i, 1) == data_face_an5(:,2));
    if sum(mask_an5_1) > 0
        face_area_an5_1 = data_face_an5(mask_an5_1, 3);
    else
        face_area_an5_1 = 0;
    end
    if sum(mask_an5_2) > 0
        face_area_an5_2 = data_face_an5(mask_an5_2, 3);
    else
        face_area_an5_2 = 0;
    end
    face_area_an5 = face_area_an5_1 + face_area_an5_2;
    % ----- get face_area_diff -----
    face_area_an4_m_an5(i) = face_area_an4 - face_area_an5;
end


%% #################################### Calc energy_gradient Defined By Grain Face Neighborhood Properties ####################################
% """
% Hypothesis: geometry/topology changes in neighbors of grain faces will drag of the face to move.
% """  
% ------------------ calculate area for all grain faces ------------------
[faces_an4_all, faces_an5_all, face_corresp] = trackFace(file_an4, file_an5, look_up_table, 0);
face_tmp_an4 = calcFaceItgCurv(file_an4, faces_an4_all, 'all_faces');
face_tmp_an5 = calcFaceItgCurv(file_an5, faces_an5_all, 'all_faces');

% ------------------ calculate area change of grain faces ------------------
face_corresp = sortrows(face_corresp, 1);
tracked_faces_an4 = face_corresp(:,1);
face_corresp = containers.Map(face_corresp(:,1),face_corresp(:,2));
area_diff_an5_m_an4 = - face_tmp_an4(:,1);
for i = 1:size(faces_an4_all, 1)
    if ismember(i, tracked_faces_an4)
        area_diff_an5_m_an4(i) = area_diff_an5_m_an4(i) + face_tmp_an5(face_corresp(i), 1);
    end
end

% ------------------ keep only unique face data in an4 ------------------
faces_an4_all = sort(faces_an4_all, 2);
[~, idx_unique_an4] = unique(faces_an4_all, 'rows');
faces_an4_all = faces_an4_all(idx_unique_an4, :);
area_diff_an5_m_an4 = area_diff_an5_m_an4(idx_unique_an4);

% ------------------ max & avg area change for the neighboring faces ------------------
fnnf_maxdec_r_m_l = zeros(size(tracked_uniqueface_an4, 1), 1);
fnnf_avgdec_r_m_l = zeros(size(tracked_uniqueface_an4, 1), 1);
for i = 1:size(tracked_uniqueface_an4, 1)
    nnf_an4_left = nn_faces_an4{i, 1};
    mask_lnnf_an4 = ismember(faces_an4_all, nnf_an4_left, 'rows');
    nnf_area_maxdec_left = min(area_diff_an5_m_an4(mask_lnnf_an4, :));
    nnf_area_avgdec_left = sum(area_diff_an5_m_an4(mask_lnnf_an4, :))/sum(mask_lnnf_an4);
    nnf_an4_right = nn_faces_an4{i, 2};
    mask_rnnf_an4 = ismember(faces_an4_all, nnf_an4_right, 'rows');
    nnf_area_maxdec_right = min(area_diff_an5_m_an4(mask_rnnf_an4, :));
    nnf_area_avgdec_right = sum(area_diff_an5_m_an4(mask_rnnf_an4, :))/sum(mask_rnnf_an4);
    
    fnnf_maxdec_r_m_l(i) = nnf_area_maxdec_right - nnf_area_maxdec_left;
    fnnf_avgdec_r_m_l(i) = nnf_area_avgdec_right - nnf_area_avgdec_left;
end
 

%% #################################### Write CSV file ####################################



fileID = fopen('190318_mig_signs.txt','w');
fprintf(fileID,'%s ,%s, %s, %s, %s, %s, %s, %s\n','gsize_diff_an4', 'gf_diff_an4', 'g_itgcurv_diff_an4', ...
                'face_itgcurv_an4_left', 'face_area_an4_m_an5', 'fnnf_maxdec_r_m_l', 'fnnf_avgdec_r_m_l', 'move_left');
for i = 1:length(move_left)
    fprintf(fileID, '%6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6d\n', data_grain_an4_diff(i,1), ...
        data_grain_an4_diff(i,2), data_grain_an4_diff(i,3), face_itgcurv_an4_left(i), face_area_an4_m_an5(i), ...
        fnnf_maxdec_r_m_l(i), fnnf_avgdec_r_m_l(i), move_left(i));
end
fclose(fileID);
    






