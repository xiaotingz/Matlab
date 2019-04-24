% ##########################################################################
% * Input
%     - tracked_uniqueface_ = [n,2]  
%           returned by trackFace.m or trackUniqueFace.m
%     - data_grain_ = [m, k]
%           m = #grains, k = #data_grain passed in
%           [size, F, F - <Fnn>, integral curvature], which  actually includes grain property and grain neighborhood property
%     - data_face = [n, 4]
%           [label_A, label_B, area, itg_curv/area]
%           returned by calcFaceCurvature.m in /Grain Curvature
% * Output
%     - all differences / decreases are calculated as right_grain - left_grain
%     - data_grain_an4_diff = [gV_diff_an4, gF_diff_an4, gFnnF_diff_an4, gMs_diff_an4]
%     - fMs_an4_left
%     - FnnF_maxdec, FnnF_avgdec'
% * NOTE
%     - This function is to calculate the energy gradient features.
% """
% TODO: there are two ways to define data_grain relavance: just an4 only, or
% use data_grain difference. now is using an4 only. 
% """
% ##########################################################################
% ----------------------- load debug data -----------------------
% file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
% file_an5 = '/Volumes/XIAOTING/Ni/An5new6_Hsmooth.dream3d';
% load('/Volumes/XIAOTING/Ni/working/181107.mat', 'tracked_uniqueface_an4', 'tracked_uniqueface_an5');
% file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/D3Ds/An4new6_fixOrigin3_smooth.dream3d');
% file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/D3Ds/An5new6_smooth.dream3d');
file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin3_Hsmooth.dream3d');
file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_Hsmooth.dream3d');
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/181107.mat', 'tracked_uniqueface_an4', 'tracked_uniqueface_an5');
load('look_up_table_an4_an5.mat')

% """
% Threshold values for determining if a triangle is good. Use in calcFaceToGrainCentroidDist.m
% """
eps_curv = 1;
eps_area = 7;
eps_min_ang = 10;

% % """
% % run  /Grain Curvature/G_F_mF.m to get the following data: from data_grain & F_mF_diff 
% % """
% data_face_an4 = data_face;
% data_grain_an4 = [data_grain(:,2), data_grain(:,3), F_mF_diff, data_grain(:,5)];
data_grain_an4 = [data_grain_an4(:,2), data_grain_an4(:,3), F_mF_diff_an4, data_grain_an4(:,5)];
data_grain_an5 = [data_grain_an5(:,2), data_grain_an5(:,3), F_mF_diff_an5, data_grain_an5(:,5)];
% clearvars -except data_face_an4 data_grain_an4 file_an4 
% ---------------------------------------------------------------
%% #################################### Make label_an4 & label_an5 one-to-one ####################################
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

%% #################################### Check If Grain Face Moved Left ####################################
% ----- Absolute distance -----
dists_f_g_an4 = calcFaceToGrainCentroidDist(file_an4, tracked_uniqueface_an4);
dists_f_g_an5 = calcFaceToGrainCentroidDist(file_an5, tracked_uniqueface_an5);
% ----- Convert distance to portion -----
total_dists_an4 = sum(dists_f_g_an4, 2);
total_dists_an5 = sum(dists_f_g_an5, 2);
dists_f_g_an4 = dists_f_g_an4 ./ total_dists_an4;
dists_f_g_an5 = dists_f_g_an5 ./ total_dists_an5;

eps_motion = 1e-2;
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


%% #################################### energy_gradient by Grain Properties & Grain Neighborhood Properties ####################################
% """
% - eps = 0.05 * min(data_grain_left, data_grain_right) is a reasonable
%   choice for converting to categorical values. 
%       However, probably should keep the values of grain_property_diff as
%       input of regression is probably better. 
% - data_grain = [size, F, F-<Fnn>, integral curvature]
% """
data_grain_an4_left = data_grain_an4(tracked_uniqueface_an4(:,1), :);
data_grain_an4_right = data_grain_an4(tracked_uniqueface_an4(:,2), :);

% """
% - diff = right - left, becase move left <-> left grain small
% - size big <-> itg_curv small (more negative)
% """
data_grain_an4_diff = data_grain_an4_right - data_grain_an4_left;
data_grain_an4_diff(:,4) = - data_grain_an4_diff(:,4);


%% #################################### energy_gradient by Grain Face Curvature ####################################
% """
% Hypothesis: if grain face convex when taking the left grain as resident grain, then grain face should move left.
% Positive cases: move left; left grain small; face area decrease
% """  
face_itgcurv_an4_left = zeros(size(tracked_uniqueface_an4, 1), 1);
tmp_face_itgcurv_an4_left = zeros(size(tracked_uniqueface_an4, 1), 1);

for i = 1:size(tracked_uniqueface_an4, 1)
    % ----- get itg_curv for the specific grain face -----
    % """
    % The reason for the long code was some faces have only [A, B] but not [B, A]. []*[] rises error.
    % """
    mask_an4_1 = (tracked_uniqueface_an4(i, 1) == data_face_an4(:,1) & tracked_uniqueface_an4(i, 2) == data_face_an4(:,2));
    mask_an4_2 = (tracked_uniqueface_an4(i, 2) == data_face_an4(:,1) & tracked_uniqueface_an4(i, 1) == data_face_an4(:,2));
    if sum(mask_an4_1) > 0
        face_itg_curv_1 = data_face_an4(mask_an4_1, 3)*data_face_an4(mask_an4_1, 4);
        face_area_an4_1 = data_face_an4(mask_an4_1, 3);
    else
        face_itg_curv_1 = 0.0;
        face_area_an4_1 = 0.0;
    end
    if sum(mask_an4_2) > 0
        face_itg_curv_2 = data_face_an4(mask_an4_2, 3)*data_face_an4(mask_an4_2, 4);
        face_area_an4_2 = data_face_an4(mask_an4_2, 3);
    else
        face_itg_curv_2 = 0.0;
        face_area_an4_2 = 0.0;
    end
    face_itgcurv_an4_left(i)  = - face_itg_curv_1 + face_itg_curv_2;
end


%% #################################### energy_gradient by Grain Face Neighborhood Properties ####################################
% """
% Hypothesis: geometry/topology changes in neighbors of grain faces will drag of the face to move.
% Can't assign sign to curvature because the sign of untracked grain faces would be random.
% """  
% load('/Volumes/XIAOTING/Ni/working/Grain Tracking/data/190221_mig_sign.mat', 'nn_faces_an4', 'nn_faces_an5')
load('/Users/xiaotingzhong/Documents/Matlab/Grain Tracking/data/190221_mig_sign.mat', 'nn_faces_an4', 'nn_faces_an5')

% ------------------ calculate area for all grain faces ------------------
[faces_an4_all, faces_an5_all, face_corresp] = trackFace(file_an4, file_an5, look_up_table, 0);
face_tmp_an4 = calcFaceItgCurv(file_an4, faces_an4_all, 'all_faces', eps_curv, eps_area, eps_min_ang);
face_tmp_an5 = calcFaceItgCurv(file_an5, faces_an5_all, 'all_faces', eps_curv, eps_area, eps_min_ang);


% ------------------ calculate area change of grain faces ------------------
% """
% diff = an5 - an4;
% """
face_corresp = sortrows(face_corresp, 1);
tracked_faces_an4 = face_corresp(:,1);
face_corresp = containers.Map(face_corresp(:,1),face_corresp(:,2));
area_diff = - face_tmp_an4(:,1);
itgcurv_diff = - face_tmp_an4(:,2);
for i = 1:size(faces_an4_all, 1)
    if ismember(i, tracked_faces_an4)
        area_diff(i) = area_diff(i) + face_tmp_an5(face_corresp(i), 1);
        itgcurv_diff(i) = itgcurv_diff(i) + face_tmp_an5(face_corresp(i), 2);
    end
end

%%
% ------------------ keep only unique face data in an4 ------------------
% faces_an4_all = sort(faces_an4_all, 2);
[~, idx_unique_an4] = unique(faces_an4_all, 'rows');
faces_an4_all = faces_an4_all(idx_unique_an4, :);
area_diff = area_diff(idx_unique_an4);
itgcurv_diff = itgcurv_diff(idx_unique_an4);

%%
% ------------------ max & avg area change for the neighboring faces ------------------
fnnf_area_maxdec = zeros(size(tracked_uniqueface_an4, 1), 1);
fnnf_area_avgdec = zeros(size(tracked_uniqueface_an4, 1), 1);
fnnf_itgcurv_maxdec = zeros(size(tracked_uniqueface_an4, 1), 1);
fnnf_itgcurv_avgdec = zeros(size(tracked_uniqueface_an4, 1), 1);
for i = 1:size(tracked_uniqueface_an4, 1)
    nnf_an4_left = nn_faces_an4{i, 1};
    mask_lnnf_an4 = ismember(faces_an4_all, nnf_an4_left, 'rows');
    nnf_area_maxdec_left = min(area_diff(mask_lnnf_an4, :));
    nnf_area_avgdec_left = sum(area_diff(mask_lnnf_an4, :))/sum(mask_lnnf_an4);
    nnf_itgcurv_maxdec_left = min(itgcurv_diff(mask_lnnf_an4, :));
    nnf_itgcurv_avgdec_left = sum(itgcurv_diff(mask_lnnf_an4, :))/sum(mask_lnnf_an4);
    nnf_an4_right = nn_faces_an4{i, 2};
    mask_rnnf_an4 = ismember(faces_an4_all, nnf_an4_right, 'rows');
    nnf_area_maxdec_right = min(area_diff(mask_rnnf_an4, :));
    nnf_area_avgdec_right = sum(area_diff(mask_rnnf_an4, :))/sum(mask_rnnf_an4);
    nnf_itgcurv_maxdec_right = min(itgcurv_diff(mask_rnnf_an4, :));
    nnf_itgcurv_avgdec_right = sum(itgcurv_diff(mask_rnnf_an4, :))/sum(mask_rnnf_an4);
    
    fnnf_area_maxdec(i) = nnf_area_maxdec_right - nnf_area_maxdec_left;
    fnnf_area_avgdec(i) = nnf_area_avgdec_right - nnf_area_avgdec_left;
    fnnf_itgcurv_maxdec(i) = nnf_itgcurv_maxdec_right - nnf_itgcurv_maxdec_left;
    fnnf_itgcurv_avgdec(i) = nnf_itgcurv_avgdec_right - nnf_itgcurv_avgdec_left;
end
 



%% #################################### Write txt file ####################################
% """
% - Gradient defined as differences / decreases from (right - left). 
% - data_grain_diff = [size, F, F-Fnn, integral curvature]
% - fMs_an4_left: face itg_curv. the most important information being its sign.
% """
fileID = fopen('190421_features_energygrad.txt','w');
fprintf(fileID,'%s ,%s, %s ,%s, %s, %s, %s, %s, %s, %s\n', 'move_left', 'gV_diff_an4', 'gF_diff_an4', 'gFnnF_diff_an4', 'gMs_diff_an4', ...
                'fMs_an4_left', 'FnnF_A_maxdec', 'FnnF_A_avgdec', 'FnnF_fMs_maxdec', 'FnnF_fMs_avgdec');
for i = 1:length(data_grain_an4_diff)
    fprintf(fileID, '%6d, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f \n', ...
        move_left(i), data_grain_an4_diff(i,1), data_grain_an4_diff(i,2), data_grain_an4_diff(i,3), data_grain_an4_diff(i,4), ...
        face_itgcurv_an4_left(i), fnnf_area_maxdec(i), fnnf_area_avgdec(i), fnnf_itgcurv_maxdec(i), fnnf_itgcurv_avgdec(i));
end
fclose(fileID);




%% ###################### Display data & Make plot ##########################
% load('/Volumes/XIAOTING/Ni/working/Grain Tracking/projection/190404_proj_dists.mat', ...
%     'face_onepiece', 'mig_localnorm_nearest', 'mig_localnorm_ot');
load('/Volumes/XIAOTING/Ni/working/Grain Tracking/projection/190417_dist_left_from_an4normproj.mat');
load('/Volumes/XIAOTING/Ni/working/Grain Tracking/tmp_190416.mat')

% daspect([1 1 1]);
% ##### plot the centroids #####
rng('shuffle');
% i = randi(size(face_onepiece, 1), 1);  
i = 1638;
idx = face_onepiece(i);

eps = 0.5;
if dist_left(idx) > eps
    memo = '   (move left)';
elseif dist_left(idx) < -eps
    memo = '   (move right)';
else
    memo = '   (stable)';
end
if gV_diff_an4(idx) > 0
    memo2 = 'right big,   ';
else
    memo2 = 'left big,   ';
end

disp(' ');
disp('#############################')
disp(['Pair ', num2str(idx), '  ;  one_piece_id = ', num2str(i)]);
dispFacePairInfo(file_an4, file_an5, tracked_uniqueface_an4, tracked_uniqueface_an5, idx)
% disp(['mig_sign_norm_nearest = ',num2str(mig_localnorm_nearest(i, 3:4))])
% disp(['mig_sign_norm_ot = ',num2str(mig_localnorm_ot(i, 3:4))])
disp('---------------')
disp(['dist_left_norm_proj = ',num2str(dist_left(idx)), memo])
disp('---------------')
disp(['mig_left_centrod = ', num2str(move_left(idx)), ...
    ',  left_centrod_an4 = ',num2str(dists_f_g_an4(idx, 1)), ',  left_centrod_an5 = ',num2str(dists_f_g_an5(idx, 1))])
disp('---------------')
disp([memo2, 'gV_diff_an4 = ', num2str(gV_diff_an4(idx)), ',  fMs_an4_left = ',num2str(fMs_an4_left(idx, 1))])
disp('---------------')
disp(['gF_diff_an4 = ', num2str(gF_diff_an4(idx)), ',  gMs_diff_an4 = ',num2str(gMs_diff_an4(idx, 1))])


figure
plotCentroids(file_an4, file_an5, tracked_uniqueface_an4, tracked_uniqueface_an5, idx)
% if move_left(idx) > 0
%     disp('move left')
%     title('move left', 'fontSize', 18)
% elseif move_left(idx) < 0
%     disp('move right')
%     title('move right', 'fontSize', 18)
% else
%     disp('stable')
%     title('stable', 'fontSize', 18)
% end

%% ##### save plot #####
f_name = ['centroids_pair_', num2str(idx)];
print(f_name, '-dtiff', '-r300')

%% ###### Unify Curvature Direction For Paraview ######
file = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin2_smooth_curveswaped.dream3d';
curv = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures')';
fl = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
mask_inner = all(fl > 0, 2);
fl_inner = fl(mask_inner, :);
curv_inner = curv(mask_inner, :);

% """ 
% if obj_grain on left of face_label, tri_curve = - tri_curve
% """
for i = 1:size(tracked_uniqueface_an4, 1)
    disp(i)
%     mask = ismember(fl, tracked_uniqueface_an4(i,:), 'rows');
    mask = (fl_inner(:,1) == tracked_uniqueface_an4(i,1) & fl_inner(:,2) == tracked_uniqueface_an4(i,2));
    curv_inner(mask) = - curv_inner(mask);
end
    
curv(mask_inner) = curv_inner;
h5write(file, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures', curv');

