% ##########################################################################
% * Required data
%     - Go to /Topology and prepare
%         triple_line_full
%         num_corners
% * Notes
%     - This script is to prepare features related to the mean field defined
%       the nearest neighbors (faces or grains) of the target face. 
%         1. Connection / Neighborhood (Face), Initial Status
%                 1.1. C - <Cnn>
%                 1.2. GF - <GFnn> 
%                 1.3. fraction of neighboring faces: twin & positive curvature
%         2. Connection / Neighborhood (Face), Change
%                 2.1. sum(DAnn)
%                 2.2. fraction of neighboring faces: growed, disappeared, appeared
%         3. Connection / Neighborhood (Grain) Changes
%                 3.1. fraction of neighboring grains: disappeared, appeared
%                 3.2. max(\delta Fnn), min(\delta Fnn) & avg(\delta Fnn)
%         4. Write txt File
%         5. Checks
%     - This script is closely related to featurePrepGeoAndTopo.m
%       Preparation of the topology features, see main.m in /Topologies
% ##########################################################################
% file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
% file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4.dream3d';
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190624_Ni_crop_TopologyResult_uniqiue.mat', ...
%     'triple_line_full_an4')
load('/Volumes/XIAOTING/Ni/working/190624_Ni_crop_TopologyResult_uniqiue.mat', ...
    'triple_line_full_an4', 'num_corners_an4', 'triple_line_full_an5', 'num_corners_an5')
file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixedOrigin_smooth.dream3d';
file_an5 = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4.dream3d';
load('look_up_table_an4_an5crop.mat')
eps_curv = 1;
eps_area = 7;
eps_min_ang = 10;

% ---------------------------- Prepare data ----------------------------
% """ Remember triple_line_full and num_corners from Topologies/main.m"""
num_neigh_an4 = h5read(file_an4, '/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')';
num_neigh_an5 = h5read(file_an5, '/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')';
num_neigh_an4(1) = [];
num_neigh_an5(1) = [];

faces_an4_inner = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
faces_an5_inner = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
faces_an4_inner = faces_an4_inner(all(faces_an4_inner > 0, 2), :);
faces_an5_inner = faces_an5_inner(all(faces_an5_inner > 0, 2), :);

nn_faces_an4 = getLeftRightFaceConnections(faces_an4_inner, triple_line_full_an4, file_an4, 'use_inner_faces');
nn_faces_an5 = getLeftRightFaceConnections(faces_an5_inner, triple_line_full_an5, file_an5, 'use_inner_faces');

face_itgcurv_an4 = calcFaceItgCurv(file_an4, faces_an4_inner, 'as_given', eps_curv, eps_area, eps_min_ang);
face_itgcurv_an5 = calcFaceItgCurv(file_an5, faces_an5_inner, 'as_given', eps_curv, eps_area, eps_min_ang);
face_area_an4 = face_itgcurv_an4(:, 1);
face_itgcurv_an4 = face_itgcurv_an4(:, 2);
face_area_an5 = face_itgcurv_an5(:, 1);
face_itgcurv_an5 = face_itgcurv_an5(:, 2);
dg_obj(:,:,1) = AAToG(60, [1, 1, 1]);
twin_an4 = calcDistFromMisorientation(file_an4, faces_an4_inner, dg_obj);
twin_an4 = int32(twin_an4 < 5);

[tracked_uniqueface_an4_complete, tracked_uniqueface_an5_complete] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_complete_faces');
[tracked_uniqueface_an4_inner, tracked_uniqueface_an5_inner] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_inner_faces');

%%
% ---------------------------- Prepare dictionaries, faces ----------------------------
dict_fl_init2final = containers.Map();
dict_fl_final2init = containers.Map();
for i = 1:size(tracked_uniqueface_an4_inner, 1)
    key_an4 = mat2str(tracked_uniqueface_an4_inner(i, :));
    key_an5 = mat2str(tracked_uniqueface_an5_inner(i, :));
    dict_fl_init2final(key_an4) = tracked_uniqueface_an5_inner(i, :);
    dict_fl_final2init(key_an5) = tracked_uniqueface_an4_inner(i, :);
end

dict_corner_an4 = containers.Map('KeyType','char','ValueType','double');
dict_area_an4 = containers.Map('KeyType','char','ValueType','double');
dict_itgfacecurv_an4 = containers.Map('KeyType','char','ValueType','double');
dict_twin_an4 = containers.Map('KeyType','char','ValueType','int32');
dict_neigh_faces_an4 = containers.Map();
for i = 1:size(faces_an4_inner, 1)
    key = mat2str(faces_an4_inner(i, :));
    dict_corner_an4(key) = num_corners_an4(i);
    dict_area_an4(key) = face_area_an4(i);
    dict_itgfacecurv_an4(key) = face_itgcurv_an4(i);
    dict_twin_an4(key) = twin_an4(i);
    dict_neigh_faces_an4(key) = [nn_faces_an4{i, 1}; nn_faces_an4{i, 2}];
end

dict_corner_an5 = containers.Map('KeyType','char','ValueType','double');
dict_area_an5 = containers.Map('KeyType','char','ValueType','double');
dict_itgfacecurv_an5 = containers.Map('KeyType','char','ValueType','double');
dict_neigh_faces_an5 = containers.Map();
for i = 1:size(faces_an5_inner, 1)
    key = mat2str(faces_an5_inner(i, :));
    dict_corner_an5(key) = num_corners_an5(i);
    dict_area_an5(key) = face_area_an5(i);
    dict_itgfacecurv_an5(key) = face_itgcurv_an5(i);
    dict_neigh_faces_an5(key) = [nn_faces_an5{i, 1}; nn_faces_an5{i, 2}];
end

% ---------------------------- Prepare dictionaries, grains ----------------------------
dict_g_init2final = containers.Map(look_up_table(:,1), look_up_table(:,2));
dict_g_final2init = containers.Map(look_up_table(:,2), look_up_table(:,1));
%%  #################################### 1. Connection / Neighborhood (Face), Initial Status #################################### 
% """
% Only focus on faces in the initial state. Namely, the disappeared faces
% will be tracked but the newly appeared faces won't.
% """
obj_faces_an4 = tracked_uniqueface_an4_complete;

% ---------------------------- Initialize variables ----------------------------
% ----- 1.1. C - <Cnn> -----
% """ Similar to F - <Fnn> for grains, but replace F with C (#corners of face) """ 
avg_cnn = zeros(size(obj_faces_an4, 1), 1);
c_cnn = zeros(size(obj_faces_an4, 1), 1);

% ----- 1.2. GF - <GFnn> ----- 
% """ Similar to F - <Fnn> for grains, but replace F with GF (absolute value of face integral curvature) """ 
avg_gfnn = zeros(size(obj_faces_an4, 1), 1);
gf_gfnn = zeros(size(obj_faces_an4, 1), 1);

% ----- 1.3. fraction of neighboring faces: twin & positive curvature -----
pos_frac = zeros(size(obj_faces_an4, 1), 1);
twin_frac = zeros(size(obj_faces_an4, 1), 1);


% % ---------------------------- Loop all neighbors in initial state to calculate properties ----------------------------
for i = 1:size(obj_faces_an4, 1)
    key_obj_face = mat2str(obj_faces_an4(i, :));
    neighbors = dict_neigh_faces_an4(key_obj_face);
    n = size(neighbors, 1);
    
    % ----- loop neighbors -----
    total_corners = 0;
    total_itgfacecurv = 0;
    num_pos_curv_neigh = 0;
    num_twins = 0;
    for j = 1:n
        key_neigh = mat2str(neighbors(j, 1:2));
        total_corners = total_corners + dict_corner_an4(key_neigh);
        total_itgfacecurv = total_itgfacecurv + dict_itgfacecurv_an4(key_neigh);
        if dict_itgfacecurv_an4(key_neigh) > 0
            num_pos_curv_neigh = num_pos_curv_neigh + 1;
        end
        num_twins = num_twins + dict_twin_an4(key_neigh);
    end
    
    % ----- record values -----
    avg_cnn(i) = total_corners / n;
    c_cnn(i) = dict_corner_an4(key_obj_face) - avg_cnn(i);
    avg_gfnn(i) = total_itgfacecurv / n;
    gf_gfnn(i) = dict_itgfacecurv_an4(key_obj_face) - avg_gfnn(i);
    pos_frac(i) = num_pos_curv_neigh / n;
    twin_frac(i) = num_twins / n;
end



%%  #################################### 2. Connection / Neighborhood (Face), Change #################################### 
% """
% Not only the disappeared faces, but also the appeared faces need to be tracked. 
% """
obj_faces_an4 = tracked_uniqueface_an4_complete;
obj_faces_an5 = tracked_uniqueface_an5_complete;

% ---------------------------- Initialize variables ----------------------------
% ----- 2.1. sum(DAnn) -----
% """ sum(area change of nearest neighbor, or connected, faces) """
total_dann = zeros(size(obj_faces_an4, 1), 1);

% ----- 2.2. fraction of neighboring faces: growed, disappeared, appeared -----
grow_frac = zeros(size(obj_faces_an4, 1), 1);
disappeared_frac = zeros(size(obj_faces_an4, 1), 1);
appeared_frac = zeros(size(obj_faces_an4, 1), 1);

% ---------------------------- Loop all neighbors in both states to calculate properties ----------------------------
for i = 1:size(obj_faces_an4, 1)
    key_an4 = mat2str(obj_faces_an4(i, :));
    key_an5 = mat2str(obj_faces_an5(i, :));
    neighbors_an4 = dict_neigh_faces_an4(key_an4);
    neighbors_an5 = dict_neigh_faces_an5(key_an5);
    
    total_neigh_area_change = 0;
    num_grow_faces = 0;
    num_appear_faces = 0;
    num_disappear_faces = 0;
    tracked_neigh_an5 = [];
    
    % ----- Loop neighbors_an4 to check initial neighbors-----
    for j = 1:size(neighbors_an4, 1)
        key_neigh_an4 = mat2str(neighbors_an4(j, 1:2));
        
        % ----- Tracked neighbor face -----
        if isKey(dict_fl_init2final, key_neigh_an4)
            fl_an5 = dict_fl_init2final(key_neigh_an4);
            tracked_neigh_an5 = [tracked_neigh_an5; fl_an5];
            key_neigh_an5 = mat2str(fl_an5);
            area_diff = dict_area_an5(key_neigh_an5) - dict_area_an4(key_neigh_an4);
            
            total_neigh_area_change = total_neigh_area_change + area_diff;
            if area_diff > 0
                num_grow_faces = num_grow_faces + 1;
            end
        % ----- Disappeared neighbor face -----
        else
            total_neigh_area_change = total_neigh_area_change - dict_area_an4(key_neigh_an4);
            num_disappear_faces = num_disappear_faces + 1;
        end
    end
    
    % ----- Loop neighbors_an5 to check new neighbors-----
    % """
    % Note all new neighbors are new faces.
    % """
    if isempty(tracked_neigh_an5)
        appeared_candidates = neighbors_an5(:, 1:2);
    else
        appeared_candidates = setdiff(neighbors_an5(:, 1:2), tracked_neigh_an5(:, 1:2), 'rows');
    end
    n = size(neighbors_an4, 1) + size(appeared_candidates, 1);
    for j = 1:size(appeared_candidates, 1)
        key_neigh_an5 = mat2str(appeared_candidates(j, :));
        
        % ----- New neighbor from a tracked face -----
        if isKey(dict_fl_final2init, key_neigh_an5)
            fl_an4 = dict_fl_final2init(key_neigh_an5);
            key_neigh_an4 = mat2str(fl_an4);
            area_diff = dict_area_an5(key_neigh_an5) - dict_area_an4(key_neigh_an4);
            
            total_neigh_area_change = total_neigh_area_change + area_diff;
            if area_diff > 0
                num_grow_faces = num_grow_faces + 1;
            end
            
        % ----- New neighbor from a new appeared face -----
        else
            total_neigh_area_change = total_neigh_area_change + dict_area_an5(key_neigh_an5);
            num_appear_faces = num_appear_faces + 1;
        end
    end
    
    % ----- record values -----
    total_dann(i) = total_neigh_area_change;
    grow_frac(i) = num_grow_faces / n;
    disappeared_frac(i) = num_disappear_faces / n;
    appeared_frac(i) = num_appear_faces / n;
end


%%  #################################### 3. Connection / Neighborhood (Grain), Change #################################### 
% """
% The hypothesis is that large change of grains will introduce interuption
% in the local mean field.
% """
obj_faces_an4 = tracked_uniqueface_an4_complete;
obj_faces_an5 = tracked_uniqueface_an5_complete;

% ------------------ 3.1. fraction of neighboring grains: disappeared, appeared ------------------
disappeared_frac_g = zeros(size(obj_faces_an4, 1), 1);
appeared_frac_g = zeros(size(obj_faces_an4, 1), 1);

% ------------------ 3.2. max(\delta Fnn), min(\delta Fnn) & avg(\delta Fnn) ------------------
% """ 
% similar to Topology/findFaceLocalTopologyChange.m 
%     - max_dfnn(max_nndiff): max(diff(F) of connected grains)
%     - avg_dfnn(avg_nndiff): average(diff(F) of connected grains)
% """
max_dfnn_g = zeros(size(obj_faces_an4, 1), 1);
min_dfnn_g = zeros(size(obj_faces_an4, 1), 1);
avg_dfnn_g = zeros(size(obj_faces_an4, 1), 1);

for i = 1:size(obj_faces_an4, 1)
    key_an4 = mat2str(obj_faces_an4(i, :));
    key_an5 = mat2str(obj_faces_an5(i, :));
    neighbors_an4 = dict_neigh_faces_an4(key_an4);
    neighbors_an4 = unique(neighbors_an4(:, 1:2));
    neighbors_an5 = dict_neigh_faces_an5(key_an5);
    neighbors_an5 = unique(neighbors_an5(:, 1:2));
    
    num_neigh_g_disappeared = 0;
    num_neigh_g_appeared = 0;
    max_neigh_g_df = -Inf;
    min_neigh_g_df = Inf;
    total_neigh_g_df = 0;
    tracked_neigh_an5 = [];
    
    % ----- Loop neighbors_an4 to check initial neighbors-----
    for j = 1:size(neighbors_an4, 1)
        id_an4 = neighbors_an4(j);
        f_an4 = num_neigh_an4(id_an4);
        
        % ----- Tracked neighbor grain -----
        if isKey(dict_g_init2final, id_an4)
            id_an5 = dict_g_init2final(id_an4);
            tracked_neigh_an5 = [tracked_neigh_an5; id_an5];
            
            f_diff = num_neigh_an5(id_an5) - f_an4;
            total_neigh_g_df = total_neigh_g_df + f_diff;
            if f_diff > max_neigh_g_df
                max_neigh_g_df = f_diff;
            end
            if f_diff < min_neigh_g_df
                min_neigh_g_df = f_diff;
            end
        % ----- Disappeared neighbor grain -----
        else
            num_neigh_g_disappeared = num_neigh_g_disappeared + 1;
            if - f_an4 < min_neigh_g_df
                min_neigh_g_df = - f_an4;
            end
            total_neigh_g_df = total_neigh_g_df - f_an4;
        end
    end

    % ----- Loop neighbors_an5 to check final neighbors-----
    if isempty(tracked_neigh_an5)
        appeared_candidates = neighbors_an5;
    else
        appeared_candidates = setdiff(neighbors_an5, tracked_neigh_an5);
    end
    n = size(neighbors_an4, 1) + size(appeared_candidates, 1);
    for j = 1:size(appeared_candidates, 1)
        id_an5 = appeared_candidates(j);
        f_an5 = num_neigh_an5(id_an5);
        
        % ----- New neighbor from a tracked face -----
        if isKey(dict_g_final2init, id_an5)
            id_an4 = dict_g_final2init(id_an5);
            
            f_diff = f_an5 - num_neigh_an4(id_an4);
            total_neigh_g_df = total_neigh_g_df + f_diff;
            if f_diff > max_neigh_g_df
                max_neigh_g_df = f_diff;
            end
            if f_diff < min_neigh_g_df
                min_neigh_g_df = f_diff;
            end
        % ----- New neighbor from a new appeared face -----
        else
            num_neigh_g_appeared = num_neigh_g_appeared + 1;
            if  f_an5 > max_neigh_g_df
                max_neigh_g_df = f_an5;
            end
            total_neigh_g_df = total_neigh_g_df + f_an5;
        end
    end
    
    disappeared_frac_g(i) = num_neigh_g_disappeared;
    appeared_frac_g(i) = num_neigh_g_appeared;
    max_dfnn_g(i) = max_neigh_g_df;
    min_dfnn_g(i) = min_neigh_g_df;
    avg_dfnn_g(i) = total_neigh_g_df / n;
end

%%  #################################### 4. Write txt file #################################### 
fileID = fopen('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_Lsmooth_mean_field.txt','w');
fprintf(fileID,'%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s\n', ...
    'c_cnn', 'gf_gfnn', 'pos_frac', 'twin_frac', ...
    'total_dann', 'grow_frac', 'disappeared_frac', 'appeared_frac', ...
    'disappeared_frac_g', 'appeared_frac_g', 'max_dfnn_g', 'min_dfnn_g', 'avg_dfnn_g');
for i = 1:size(obj_faces_an4, 1)
    fprintf(fileID, '%6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f, %6.3f,%6.3f, %6.3f, %6.3f, %6.3f, %6.3f\n', ...
        c_cnn(i), gf_gfnn(i), pos_frac(i), twin_frac(i), ...
        total_dann(i), grow_frac(i), disappeared_frac(i), appeared_frac(i), ...
        disappeared_frac_g(i), appeared_frac_g(i), max_dfnn_g(i), min_dfnn_g(i), avg_dfnn_g(i));
end
fclose(fileID);


%%  #################################### 5. Checks #################################### 
% ###### 4.1 Check if F-Fnn relationship applies to grain_itg_curv instead of faces ######
% data_grain is calculate from /Grain Curvature/
% load('tmp.mat')
data_grain = data_grain_an4;
file = '/Volumes/XIAOTING/Ni/An4new6_fixedOrigin_smooth.dream3d';
% """ data_grain = [grainId, grainDiameter, #Faces, #edges, IntegralGrainCurvature] """
mean_field_id = 5;
% start = -300;
% width = 10;
% step = 60;
start = -80;
width = 1;
step = 160;

figure()
plotMeanField(file, data_grain, start, step, width)

% xlim([-30, 30])
% ylim([-4, 4])
