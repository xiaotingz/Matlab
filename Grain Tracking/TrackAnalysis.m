% ##################################################
% Contents
% * Make a easy to check dictionary 
% * Surface grains
%     - The Complete Tracked Grains
%     - The Grains that Growed From Being Inner to Touch Surface
% * Plot Triangle Curvature 
% * Get twins
% * Check Piecewise and Twins
% * Get FeatureFaceId for the Faces Satisfying a Certain Condition
% * Distorted Triangles
% * Visualization
% ##################################################
file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
load('look_up_table_an4_an5.mat')

%% ##### Make a easy to check dictionary #####
data = [face_corresp, face_area_diff, FCurv_diff];
data = sortrows(data, [4, 3]);

face_map_An4 = containers.map(face_corresp(:,1), num2cell(faces_an4(face_corresp(:,1),:),2));
face_map_An5 = containers.map(face_corresp(:,2), num2cell(faces_an5(face_corresp(:,2),:),2));

%% ##### The Complete Tracked Grains #####
surf_grain_an4 = h5read(file_an4,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
surf_grain_an5 = h5read(file_an5,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
surf_grain_an4(1) = [];
surf_grain_an5(1) = [];

mask_lookup = logical([surf_grain_an4(look_up_table(:,1)), surf_grain_an5(look_up_table(:,2))]);
lookup_CG = look_up_table(all(~mask_lookup, 2), :);


%% ##### The Grains that Growed From Being Inner to Touch Surface #####
[faces_an4, faces_n5, face_corresp_all] = TrackFace(file_an4, file_an5, look_up_table, false);
faceInfo_all = [faces_an4(face_corresp_all(:,1),:), faces_an5(face_corresp_all(:,2),:)];

surf_grain_an4 = h5read(file_an4,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
surf_grain_an5 = h5read(file_an5,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
surf_grain_an4(1) = [];
surf_grain_an5(1) = [];
surf_grain_info_all = [surf_grain_an4(faceInfo_all(:,1:2)), surf_grain_an5(faceInfo_all(:,3:4))];

% tmp1 = (SGinfo_all(:,1)==1 & SGinfo_all(:,2)==1 & SGinfo_all(:,3)==0 & SGinfo_all(:,4)==0);
tmp2 = (surf_grain_info_all(:,1)==0 & surf_grain_info_all(:,2)==0 & surf_grain_info_all(:,3)==1 & surf_grain_info_all(:,4)==1);
mask_growtosurf = tmp2;

% FItgCurvs_An4 = faceCurvatureForTrack(file_An4, faces_An4);
% FItgCurvs_An5 = faceCurvatureForTrack(file_An5, faces_An5);
diff_tmp = face_integ_curv_an5(face_corresp_all(:,2),:) - face_integ_curv_an4(face_corresp_all(:,1),:);
face_area_diff = diff_tmp(mask_growtosurf,1);


%% ##### Plot Triangle Curvature Distribution #####
set(0,'defaultAxesLabelFontSize',1.1)
set(0,'defaultAxesFontSize',19)

facelabel_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_curve_an4 =  abs(roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)).';
mask_an4 = ~(any(facelabel_an4 <= 0, 2) | tri_curve_an4 > 100);
tri_curve_an4 = tri_curve_an4(mask_an4);
facelabel_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_curve_an5 =  abs(roundn(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)).';
mask_an5 = ~(any(facelabel_an5 <= 0, 2) | tri_curve_an5 > 100);
tri_curve_an5 = tri_curve_an5(mask_an5);

histogram(tri_curve_an4,100)
hold on
histogram(tri_curve_an5,100)
set(gca, 'YScale', 'log')
avg_tri_curv_an4 = sum(tri_curve_an4)/length(tri_curve_an4);
avg_tri_curv_an5 = sum(tri_curve_an5)/length(tri_curve_an4);
legend(['aveTriCurv, An4 = ', num2str(avg_tri_curv_an4)], ['aveTriCurv, An5 = ', num2str(avg_tri_curv_an5)])
xlabel('Triangle Curvature, \mum^{-1}')
ylabel('# Triangles')


%% ##### Get Twins #####
% """
% * Note that twin is decided by the ||rf_vec_diff||, calculated
% rf_vec_diff(60@111, 55@111)=0.0568, so chose thres = 0.5
% * Also, the twin judged from an4 is different from the twin judged from
% an5. It's arguable that which should be used but define twin as twin_an4 | twin_an5
% """
EA_an4 = rad2deg(h5read(file_an4,'/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles')).';
EA_an5 = rad2deg(h5read(file_an5,'/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles')).';
O = CrysSym;

[tracked_uniqueface_an4, tracked_uniqueface_an5] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_complete_faces');
g1 = eye(3);
g2 = AAToG(60, [1,1,1]);
[rf_vec_obj, ~] = dgInFZ(g1, g2, O);
thres = 0.05;

rf_vec_an4 = zeros(length(tracked_uniqueface_an4), 3);
rf_vec_an5 = rf_vec_an4;
for i = 1:length(tracked_uniqueface_an4)
    g1_an4 = EAtoG(EA_an4(tracked_uniqueface_an4(i,1), :));
    g2_an4 = EAtoG(EA_an4(tracked_uniqueface_an4(i,2), :));
    g1_an5 = EAtoG(EA_an5(tracked_uniqueface_an5(i,1), :));
    g2_an5 = EAtoG(EA_an5(tracked_uniqueface_an5(i,2), :));    
    [rf_vec_an4(i, :), ~] = dgInFZ(g1_an4, g2_an4, O);
    [rf_vec_an5(i, :), ~] = dgInFZ(g1_an5, g2_an5, O);
end

rf_vec_obj = repmat(rf_vec_obj, length(rf_vec_an4), 1);
dg_diff_an4 = rf_vec_an4 - rf_vec_obj;
dg_diff_an4 = sum(dg_diff_an4.*dg_diff_an4, 2);
dg_diff_an5 = rf_vec_an5 - rf_vec_obj;
dg_diff_an5 = sum(dg_diff_an5.*dg_diff_an5, 2);
mask_twin_an4 = (dg_diff_an4 < thres);
mask_twin_an5 = (dg_diff_an5 < thres);
mask_twin = mask_twin_an4 | mask_twin_an5;

clear g1 g2 g1_an4 g2_an4 g1_an5 g2_an5 i O

tracked_facenotwin_an4 = tracked_uniqueface_an4(~mask_twin, :);
tracked_facenotwin_an5 = tracked_uniqueface_an5(~mask_twin, :);
disp(['#twins_an4 = ', num2str(sum(mask_twin_an4)), '; #twins_an5 = ', num2str(sum(mask_twin_an5))])
disp(['#twins_either =', num2str(sum(mask_twin))]);


%% ##### Check Piecewise and Twins #####
% """
% * Run the previous section (Get Twins first)
% * There are many piecewise surfaces that are not twins, checked they ----
% """

% ----- Get piecewise faces (by subgrah) -----
load('180927_tmp');
load('180822.mat', 'facelabel_an4', 'facelabel_an5', 'tri_node_an4', 'tri_node_an5');
load('180828_piecewise_face', 'face_piecewise');

twin_id = (1:length(tracked_uniqueface_an4))';
twin_id = twin_id(mask_twin);

% ----- Check size of the piecewise non-twins faces -----
seclarge_piecesize = zeros(size(face_piecewise));
for i = 1: length(face_piecewise)
    % ----- Make subgraph from the piecewise surfaces -----
    idx = face_piecewise(i);
    obj_facelabel_1 = tracked_uniqueface_an4(idx, :);
    obj_facelabel_2 = tracked_uniqueface_an5(idx, :);
    mask_objface_1 = (facelabel_an4(:,1) == obj_facelabel_1(1) & facelabel_an4(:,2) == obj_facelabel_1(2));
    mask_objface_2 = (facelabel_an5(:,1) == obj_facelabel_2(1) & facelabel_an5(:,2) == obj_facelabel_2(2));
    face_tri_nodeid_1 = tri_node_an4(mask_objface_1, :);
    face_unique_nodeid_1 = unique(face_tri_nodeid_1);
    face_tri_nodeid_2 = tri_node_an5(mask_objface_2, :);
    face_unique_nodeid_2 = unique(face_tri_nodeid_2);

    [subgraph_1, subgraph_2] = findSubgraph(face_unique_nodeid_1, face_unique_nodeid_2, face_tri_nodeid_1, face_tri_nodeid_2);
    
    % ----- Record max size of the two 2nd largest piece -----
    tmp1 = 0;
    tmp2 = 0;
    if length(unique(subgraph_1)) > 1
        tmp1 = histcounts(subgraph_1);
        tmp1 = sort(tmp1);
        tmp1 = tmp1(end - 1);
    elseif length(unique(subgraph_2)) > 1
        tmp2 = histcounts(subgraph_2);
        tmp2 = sort(tmp2);
        tmp2 = tmp2(end - 1);
    end
    seclarge_piecesize(i) = max(tmp1, tmp2);        
end

[piecewise_normal, idx_normal] = setdiff(face_piecewise, twin_id);
[piecewise_twin, idx_twin] = intersect(face_piecewise, twin_id);
disp(['#face_piecewise = ', num2str(length(face_piecewise)), '; #face_piecewise that are twins = ', num2str(length(piecewise_twin))])

piecesize_vari_normal = [piecewise_normal, seclarge_piecesize(idx_normal)];
piecesize_vari_twin = [piecewise_twin, seclarge_piecesize(idx_twin)];

clearvars -except piecesize_vari_normal piecesize_vari_twin
%% ----- Display face_pair info ----- 
% idx_piece = 569;
% idx = piecewise_notwin(idx_piece);
idx = all_multipiece_sym(101);

dispFacePairInfo(file_an4, file_an5, tracked_uniqueface_an4, tracked_uniqueface_an5, idx)

x_to_y = X_to_Y{idx};
obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
% face_node_info = getSingleFaceNodes(file_an4, obj_facelabel_an4, file_an5, obj_facelabel_an5);
face_node_info = getSingleFaceNodes(obj_facelabel_an4, obj_facelabel_an5);
visualizeFace(face_node_info, x_to_y)
hold off

% x_min = min(face_node_info{4,2}(:,1));
% x_max = max(face_node_info{4,2}(:,1));
% y_min = min(face_node_info{4,2}(:,2));
% y_max = max(face_node_info{4,2}(:,2));
% z_min = min(face_node_info{4,2}(:,3));
% z_max = max(face_node_info{4,2}(:,3));



%% ##### Get FeatureFaceId for the Faces Satisfying a Certain Condition #####
clc

file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_mesh.dream3d');
file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh.dream3d');

face_id_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
face_id_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
num_tri_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceFeatureData/NumTriangles')';
num_tri_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceFeatureData/NumTriangles')';
surf_grain_an4 = h5read(file_an4,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
surf_grain_an5 = h5read(file_an5,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures').';
face_id_an4(1,:) = [];
face_id_an5(1,:) = [];
face_id_an4 = sort(face_id_an4, 2);
face_id_an5 = sort(face_id_an5, 2);
surf_grain_an4(1) = [];
surf_grain_an5(1) = [];

tmp1 = [(1:length(FCentrs_diff))', FCentrs_diff];
sort_tmp = sortrows(tmp1, -2);

cond = true;
while cond
    obj_face = round(rand()*14008);
    tracked_uniqueface_an4 = sort(tracked_uniqueface_an4, 2);

    obj_label_an4 = tracked_uniqueface_an4(sort_tmp(obj_face, 1), :);
    [idx_an4, ~] = find(face_id_an4(:,1) == obj_label_an4(1) & face_id_an4(:,2) == obj_label_an4(2));
    obj_label_an5 = tracked_uniqueface_an5(sort_tmp(obj_face, 1), :);
    [idx_an5, ~] = find((face_id_an5(:,1) == obj_label_an5(1) & face_id_an5(:,2) == obj_label_an5(2)) | (face_id_an5(:,1) == obj_label_an5(2) & face_id_an5(:,2) == obj_label_an5(1)));
    if ~isempty(face_id_an4(idx_an4, :)) && ~isempty(face_id_an5(idx_an5, :))
        cond1 = any([surf_grain_an4(face_id_an4(idx_an4, :)); surf_grain_an5(face_id_an5(idx_an5, :))]==1);
        cond2 = (num_tri_an4(idx_an4)>300 || num_tri_an5(idx_an5)>300);
%         cond = (cond1 || cond2);
        numTri_diff = abs(num_tri_an4(idx_an4) - num_tri_an4(idx_an4));
        cond3 = numTri_diff > 0.5*num_tri_an4(idx_an4) || numTri_diff > 0.5*num_tri_an5(idx_an5);
        cond = (cond1 || cond2 || cond3);
    end
end
disp(['FaceId in An4:   ', num2str(idx_an4), ',   size=', num2str(num_tri_an4(idx_an4))])
disp(['FaceLabel in An4:  [', num2str(face_id_an4(idx_an4, :)), ']', ]);
disp(['FaceId in An5:   ', num2str(idx_an5), ',   size=', num2str(num_tri_an5(idx_an5))])
disp(['FaceLabel in An5:  [', num2str(face_id_an5(idx_an5, :)), ']', ]);


% [face_coords_an4, unique_facelabel_an4] = makeFaceCoords(file_an4);
% [face_coords_an5, unique_facelabel_an5] = makeFaceCoords(file_an5);
% unique_facelabel_corresp = trackUniqueFace(unique_facelabel_an4, unique_facelabel_an5, look_up_table);


%% ##### Record Distorted Triangles #####
% ----- the objective is to check if the an5 triangles found by correspondence are distorted -----
load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822_FaceCorresp.mat');
load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822.mat');


distort_tri_nodeid_an4 = cell(length(tracked_uniqueface_an4), 1);
distort_tri_nodeid_an5 = cell(length(tracked_uniqueface_an4), 1);
min_angle_diff = cell(length(tracked_uniqueface_an4), 1);
parfor i = 1:length(tracked_uniqueface_an4)
    % ----- get the object face triangles and nodes -----
    obj_facelabel_an4 = tracked_uniqueface_an4(i, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(i, :);
    x_to_y = X_to_Y{i};
    [distort_tri_nodeid_an4{i}, distort_tri_nodeid_an5{i}, min_angle_diff{i}] = getDistortCorrespTris(obj_facelabel_an4, obj_facelabel_an5, x_to_y, ...
        'min_angle_diff', 20, facelabel_an4, facelabel_an5, node_coord_an4, node_coord_an5, tri_node_an4, tri_node_an5);
    disp(i)
end

%% ##### Finer Classify the Distorted Tirangles #####

% ----- edge lengths -----
% load('180826_longEdges.mat')
% load('180822_FaceCorresp.mat')

% longedge_length = [];
edge_10to20_faces = [];
edge_20to50_faces = [];
edge_longer50_faces = [];
for i = 1:length(min_angle_diff)
%     longedge_length = [longedge_length; longestedge{i}];
    if any(min_angle_diff{i} > 10) && all(min_angle_diff{i} <= 20) 
        edge_10to20_faces = [edge_10to20_faces; i];
    end
    if any(min_angle_diff{i} > 20) && all(min_angle_diff{i} <= 50) 
        edge_20to50_faces = [edge_20to50_faces; i];
    end
    if any(min_angle_diff{i} > 50)
        edge_longer50_faces = [edge_longer50_faces; i];
    end
end

% ----- min_angle_diff, portion of distorted triangles -----
total_num_tri_onepiece = 0;
num_tri_minadiff20_onepiece = 0;
num_tri_minadiff30_onepiece = 0;
face_onepice = (1:length(tracked_uniqueface_an4))';
face_onepice(face_piecewise) = [];
parfor i = 1:length(face_onepice)
    idx = face_onepice(i);
    obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
    mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
    total_num_tri_onepiece = total_num_tri_onepiece + sum(mask_objface_an4);
    num_tri_minadiff20_onepiece = num_tri_minadiff20_onepiece + size(min_angle_diff{idx},1);
    num_tri_minadiff30_onepiece = num_tri_minadiff30_onepiece + sum(min_angle_diff{idx}>30);
end

total_num_tri_all = 0;
num_tri_minadiff20_all = 0;
num_tri_minadiff30_all = 0;
parfor i = 1:length(tracked_uniqueface_an4)
    obj_facelabel_an4 = tracked_uniqueface_an4(i, :);
    mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
    total_num_tri_all = total_num_tri_all + sum(mask_objface_an4);
    num_tri_minadiff20_all = num_tri_minadiff20_all + size(min_angle_diff{i},1);
    num_tri_minadiff30_all = num_tri_minadiff30_all + sum(min_angle_diff{i}>30);
end

%% ##### Visualization #####
file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');

% idx = edge_10to20_faces(randi(length(edge_10to20_faces)));
idx = 1;
x_to_y = X_to_Y{idx};
obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
face_node_info = getSingleFaceNodes(file_an4, obj_facelabel_an4, file_an5, obj_facelabel_an5);
obj_node_an4 = distort_tri_nodeid_an4{idx};
obj_node_an5 = distort_tri_nodeid_an5{idx};
edge_length = min_angle_diff{idx};
% mask_longer10 = (edge_length > 10);
% obj_node_an4 = obj_node_an4(mask_longer10, :);
% obj_node_an5 = obj_node_an5(mask_longer10, :);


% ----- Corresp -----
% visualizeFace(face_node_info, x_to_y)

% ----- Distort triangles -----
visualizeFace(face_node_info, x_to_y, obj_node_an4, obj_node_an5, 'distort_tri')
hold off





