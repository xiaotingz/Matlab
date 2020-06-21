% ##########################################################################
% * NOTE
%     - This function is to calculate the energy gradient features, which
%     is defined property difference between the left and right grain of the target face. 
%         The choice of left and right is artibary but consistent across
%         different properties. 
% * Sections
%     1. Calcualte Grain Face Moved Left
%     2. energy_gradient by Grain Properties
%     3. energy_gradient by Signed Grain Face Curvature
%     4. Write txt file
%     5. Checks
%         5.1. Check Migartion & Grain data
%         5.2. Check Face Diff
%         5.3. Unify Curvature Direction For Paraview
% * Inputs
%     - tracked_uniqueface_ = [n,2]  
%           returned by trackFace.m or trackUniqueFace.m
%     - data_grain_ = [m, k]
%           m = #grains, k = #data_grain passed in
%           [size, F, F - <Fnn>, integral curvature], which  actually includes grain property and grain neighborhood property
%     - data_face = [n, 4]
%           [label_A, label_B, area, itg_curv/area]
%           returned by calcFaceCurvature.m in /Grain Curvature
% """
% TODO: there are two ways to define data_grain relavance: just an4 only, or
% use data_grain difference. now is using an4 only. 
% """
% ##########################################################################
file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4.dream3d';
load('look_up_table_an4_an5crop.mat')

% """
% Threshold values for determining if a triangle is good. Use in calcFaceToGrainCentroidDist.m
% """
eps_curv = 1;
eps_area = 7;
eps_min_ang = 10;

%%
% % """
% % run  /Grain Curvature/G_F_mF.m to get the following data: from data_grain & F_mF_diff 
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_Lsmoth_data_grain.mat')
% % """
data_grain_an4 = [data_grain_an4(:,2), data_grain_an4(:,3), F_mF_diff_an4, data_grain_an4(:,5)];
data_grain_an5 = [data_grain_an5(:,2), data_grain_an5(:,3), F_mF_diff_an5, data_grain_an5(:,5)];



%% ##### Make label_an4 & label_an5 one-to-one #####
[tracked_uniqueface_an4_inner, tracked_uniqueface_an5_inner] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_inner_faces');
obj_faces_an4 = tracked_uniqueface_an4_inner;
obj_faces_an5 = tracked_uniqueface_an5_inner;

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
% ----- sort obj_faces_an4 it keep exact corresp with obj_faces_an5 -----
obj_faces_an4_sorted = zeros(size(obj_faces_an5));
for i = 1:size(obj_faces_an5, 1)
    for j = 1:size(obj_faces_an5, 2)
        obj_faces_an4_sorted(i, j) = look_up_table_2to1(obj_faces_an5(i,j), 1);
    end
end

% ----- check faces_an4_sorted -----
good_sort = sum(ismember([obj_faces_an4_sorted(:,1), obj_faces_an5(:,1)], look_up_table, 'rows')) == size(obj_faces_an4, 1) && ...
    sum(ismember([obj_faces_an4_sorted(:,2), obj_faces_an5(:,2)], look_up_table, 'rows')) == size(obj_faces_an4, 1);
if ~ good_sort
    disp('faces_an4_sorted were wrong!')
end

%% #################################### 1. Calcualte Grain Face Moved Left ####################################
% ----- Absolute distance -----
dists_f_g_an4 = calcFaceToGrainCentroidDist(file_an4, obj_faces_an4_sorted, eps_curv, eps_area, eps_min_ang);
dists_f_g_an5 = calcFaceToGrainCentroidDist(file_an5, obj_faces_an5, eps_curv, eps_area, eps_min_ang);
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
    elseif isnan(dists_f_g_an4(i, 1)) || isnan(dists_f_g_an5(i, 1))
        move_left(i) = NaN;
    else
        move_left(i) = 0;
    end
end


%% #################################### 2. energy_gradient by Grain Properties ####################################
% """
% - eps = 0.05 * min(data_grain_left, data_grain_right) is a reasonable
%   choice for converting to categorical values. 
%       However, probably should keep the values of grain_property_diff as
%       input of regression is probably better. 
% - data_grain = [size, F, F-<Fnn>, integral curvature]
% """
data_grain_an4_left = data_grain_an4(obj_faces_an4_sorted(:,1), :);
data_grain_an4_right = data_grain_an4(obj_faces_an4_sorted(:,2), :);
data_grain_an5_left = data_grain_an5(obj_faces_an5(:,1), :);
data_grain_an5_right = data_grain_an5(obj_faces_an5(:,2), :);

% """
% - diff = right - left, becase move left <-> left grain small
% - size big <-> itg_curv small (more negative)
% """
data_grain_an4_diff = data_grain_an4_right - data_grain_an4_left;
data_grain_an4_diff(:,4) = - data_grain_an4_diff(:,4);


%% #################################### 3. energy_gradient by Signed Grain Face Curvature ####################################
% """
% Hypothesis: if grain face convex when taking the left grain as resident grain, then grain face should move left.
% Positive cases: move left; left grain small; face area decrease
% """  
face_itgcurv_an4 = calcFaceItgCurv(file_an4, obj_faces_an4_sorted, 'signed_resident_left', eps_curv, eps_area, eps_min_ang);
face_itgcurv_an5 = calcFaceItgCurv(file_an5, obj_faces_an5, 'signed_resident_left', eps_curv, eps_area, eps_min_ang);
face_itgcurv_an4_left = face_itgcurv_an4(:, 2);
face_itgcurv_an5_left = face_itgcurv_an5(:, 2);

%% #################################### 4. Write txt file ####################################
% """
% - Gradient defined as differences / decreases from (right - left). 
% - data_grain_diff = [size, F, F-Fnn, integral curvature]
% - fMs_an4_left: face itg_curv. the most important information being its sign.
% """

dist_f_g_diff = dists_f_g_an5 - dists_f_g_an4;
% ----------------------------------- Full faces -----------------------------------
fileID = fopen('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_Hsmooth_energygrad.txt','w');
fprintf(fileID,'%s, %s, %s, %s, %s, %s\n', 'dist_f_g_diff', 'gV_diff_an4', 'gF_diff_an4', ...
    'gFnnF_diff_an4', 'gMs_diff_an4', 'fMs_an4_left');
for i = 1:length(data_grain_an4_diff)
    fprintf(fileID, '%6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f \n', ...
        dist_f_g_diff(i), data_grain_an4_diff(i,1), data_grain_an4_diff(i,2), ...
        data_grain_an4_diff(i,3), data_grain_an4_diff(i,4), face_itgcurv_an4_left(i));
end
fclose(fileID);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 5. CHECKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 1. Check Migartion & Grain data
% 2. Check Face Diff 
% 3. Unify Curvature Direction For Paraview

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

%% ######################################## 5.1. Check Migartion & Grain data ########################################
% load('/Volumes/XIAOTING/Ni/working/181107_mig_piececorresp_comb.mat', 'face_onepiece');
% load('/Volumes/XIAOTING/Ni/working/Grain Tracking/projection/190417_dist_left_from_an4normproj.mat');
file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4.dream3d';

gV_diff_an4 = data_grain_an4_diff(:,1);
fMs_an4_left = face_itgcurv_an4_left;
fMs_an5_left = face_itgcurv_an5_left;

% daspect([1 1 1]);
% ##### plot the centroids #####
rng('shuffle');
idx = randi(size(obj_faces_an4_sorted, 1), 1);  
% i = 1638;

% eps = 0.5;
% if dist_left(idx) > eps
%     memo = '   (move left)';
% elseif dist_left(idx) < -eps
%     memo = '   (move right)';
% else
%     memo = '   (stable)';
% end

disp(' ');
disp('#############################')
disp(['Pair ', num2str(idx)]);
dispFacePairInfo(file_an4, file_an5, obj_faces_an4_sorted, obj_faces_an5, idx)
% disp(['mig_sign_norm_nearest = ',num2str(mig_localnorm_nearest(i, 3:4))])
% disp(['mig_sign_norm_ot = ',num2str(mig_localnorm_ot(i, 3:4))])
% disp('---------------')
% disp(['dist_left_norm_proj = ',num2str(dist_left(idx)), memo])
disp('---------------')
disp(['mig_left_centroid = ', num2str(move_left(idx)), ...
    ',  left_centroid_an4 = ',num2str(dists_f_g_an4(idx, 1)), ',  left_centroid_an5 = ',num2str(dists_f_g_an5(idx, 1))])


figure
plotCentroids(file_an4, file_an5, obj_faces_an4_sorted, obj_faces_an5, idx)
if move_left(idx) > 0
    disp('move left')
    title('move left', 'fontSize', 18)
elseif move_left(idx) < 0
    disp('move right')
    title('move right', 'fontSize', 18)
else
    disp('stable')
    title('stable', 'fontSize', 18)
end


% ######################################## 5.2. Check Face Itg Curv By Calculate From Raw ########################################
% rng('shuffle');
% idx = randi(size(faces_an4_sorted, 1), 1); 

% ----- Display data from above calculation -----
if gV_diff_an4(idx) > 0
    memo2 = 'right big,   ';
else
    memo2 = 'left big,   ';
end

disp('---------------')
disp([memo2, 'gV_diff_an4 = ', num2str(gV_diff_an4(idx)), ',  fMs_an4_left = ',num2str(fMs_an4_left(idx, 1)), ',  fMs_an5_left = ',num2str(fMs_an5_left(idx, 1))])
% disp('---------------')
% disp(['gF_diff_an4 = ', num2str(gF_diff_an4(idx)), ',  gMs_diff_an4 = ',num2str(gMs_diff_an4(idx, 1))])


% % ----- Display data from raw calculation ----
% label_an4 = obj_faces_an4_sorted(idx, :);
% mask_1 = (tri_fl(:,1)==label_an4(1) & tri_fl(:,2)==label_an4(2));
% mask_2 = (tri_fl(:,1)==label_an4(2) & tri_fl(:,2)==label_an4(1));
% 
% itg_curv_left = - sum(tri_area(mask_1) .* tri_curv(mask_1)) + sum(tri_area(mask_2) .* tri_curv(mask_2));
% 
% disp(['check from raw:  fMs_an4_left = ', num2str(itg_curv_left)]);


%% ##### save plot #####
f_name = ['centroids_pair_', num2str(idx)];
print(f_name, '-dtiff', '-r300')




%% ######################################## 5.3. Unify Curvature Direction For Paraview ########################################
file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth_forParaview.dream3d';
curv = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures')';
fl = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
mask_inner = all(fl > 0, 2);
fl_inner = fl(mask_inner, :);
curv_inner = curv(mask_inner, :);

% """ 
% if obj_grain on left of face_label, tri_curve = - tri_curve
% """
obj_faces = obj_faces_an4_sorted;
for i = 1:size(obj_faces, 1)
%     mask = ismember(fl, tracked_uniqueface_an4(i,:), 'rows');
    mask = (fl_inner(:,1) == obj_faces(i,1) & fl_inner(:,2) == obj_faces(i,2));
    curv_inner(mask) = - curv_inner(mask);
end
    
curv(mask_inner) = curv_inner;
h5write(file, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures', curv');



