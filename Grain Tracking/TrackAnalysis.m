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

face_label_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_curve_an4 =  abs(roundn(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)).';
mask_an4 = ~(any(face_label_an4 <= 0, 2) | tri_curve_an4 > 100);
tri_curve_an4 = tri_curve_an4(mask_an4);
face_label_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_curve_an5 =  abs(roundn(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)).';
mask_an5 = ~(any(face_label_an5 <= 0, 2) | tri_curve_an5 > 100);
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


%% ##### Plot Triangle Curvature Distribution #####
% get the coords of one face pair

i = 1;
idx_an4 = face_corresp(i, 1);
idx_an5 = face_corresp(i, 2);
label_an4 = faces_an4(idx_an4, :);
label_an5 = faces_an5(idx_an5, :);

face_label_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_nodes_an4 = 1 + double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
num_coords_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
face_label_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_nodes_an5 = 1 + double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
num_coords_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';

mask_an4 = ((face_label_an4(:,1) == label_an4(1) & face_label_an4(:,2) == label_an4(2)) | (face_label_an4(:,1) == label_an4(2) & face_label_an4(:,2) == label_an4(1)));
mask_an5 = ((face_label_an5(:,1) == label_an5(1) & face_label_an5(:,2) == label_an5(2)) | (face_label_an5(:,1) == label_an5(2) & face_label_an5(:,2) == label_an5(1)));

this_face_nodes_an4 = tri_nodes_an4(mask_an4,:);
this_face_nodes_an5 = tri_nodes_an5(mask_an5,:);
this_face_nodes_an4 = reshape(this_face_nodes_an4, [], 1);
this_face_nodes_an5 = reshape(this_face_nodes_an5, [], 1);
tri_node_coords_an4 = num_coords_an4(this_face_nodes_an4,:);
tri_node_coords_an5 = num_coords_an5(this_face_nodes_an5,:);

tmp_nodes_an4 = tri_nodes_an4(mask_an4,:);
tmp_nodes_an5 = tri_nodes_an5(mask_an5,:);


%% ##### Visualize Grain Face in MATLAB #####
trisurf( tmp_nodes_an4, num_coords_an4(:,1), num_coords_an4(:,2), num_coords_an4(:,3) );
rotate3d on
hold on
trisurf( tmp_nodes_an5, num_coords_an5(:,1), num_coords_an5(:,2), num_coords_an5(:,3),'Facecolor',[218/255 168/255 32/255] );
rotate3d on




%% ##### Find the FaceID of A Face Between Two Grains #####
file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_mesh.dream3d');
file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh.dream3d');
face_id = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FeatureFaceId')';
face_label = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
face_label_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
face_label_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
face_id_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FeatureFaceId')';
face_id_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceData/FeatureFaceId')';

obj_face = 1;


[~, ~, faceCorresp_CF] = TrackFace(file_an4, file_an5, look_up_table, true);

face_info_CF = [C_An4, faces_an4(faceCorresp_CF(:,1),:)];
face_info_CF = sortrows(face_info_CF);

GA = face_info_CF(68,2);
GB = face_info_CF(68,3);
mask = ((face_label(:,1) == GA & face_label(:,2) == GB) | (face_label(:,1) == GB & face_label(:,2) == GA));
% unique(face_id(mask))

tmp = [(1:length(FCentrs_diff))', FCentrs_diff];
sort_tmp = sortrows(tmp, -2);

obj_label_an5 = tracked_facelabel_an5(sort_tmp(obj_face, 1), :);
mask = (face_label_an5(:,1) == obj_label_an5(1) & face_label_an5(:,2) == obj_label_an5(2));
unique(face_label_an5(mask, :), 'rows')
unique(face_id_an5(mask))
obj_label_an4 = tracked_facelabel_an4(sort_tmp(obj_face, 1), :);
mask = (face_label_an4(:,1) == obj_label_an4(1) & face_label_an4(:,2) == obj_label_an4(2));
unique(face_label_an4(mask, :), 'rows')
unique(face_id_an4(mask))



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

tmp = [(1:length(FCentrs_diff))', FCentrs_diff];
sort_tmp = sortrows(tmp, -2);

cond = true;
while cond
    obj_face = round(rand()*14008);
    trackedFL_An4 = sort(trackedFL_An4, 2);

    obj_label_an4 = trackedFL_An4(sort_tmp(obj_face, 1), :);
    [idx_an4, ~] = find(face_id_an4(:,1) == obj_label_an4(1) & face_id_an4(:,2) == obj_label_an4(2));
    obj_label_an5 = trackedFL_An5(sort_tmp(obj_face, 1), :);
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


%% ##### Check If the an5 Triangles Found By Correspondence are Distorted #####
load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822_FaceCorresp.mat');
load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822.mat');


longedge_tri_nodeid_an4 = cell(length(tracked_uniqueface_an4), 1);
longedge_tri_nodeid_an5 = cell(length(tracked_uniqueface_an4), 1);
longestedge = cell(length(tracked_uniqueface_an4), 1);
for i = 1:2
    % ----- get the object face triangles and nodes -----
    obj_facelabel_an4 = tracked_uniqueface_an4(i, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(i, :);
    x_to_y = X_to_Y{i};
    [longedge_tri_nodeid_an4{i}, longedge_tri_nodeid_an5{i}, longestedge{i}] = getLongEdgeCorrespTris(obj_facelabel_an4, obj_facelabel_an5, x_to_y, 6);
end


% idx = 1;
% % ----- get the object face triangles and nodes -----
% obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
% obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
% x_to_y = X_to_Y{idx};
% 
% face_node_info = getSingleFaceNodes(file_an4, obj_facelabel_an4, file_an5, obj_facelabel_an5);
% plot_trinode_an4 = node_coord_an4(longedge_tri_nodeid_an4, :);
% plot_trinode_an5 = node_coord_an5(longedge_tri_nodeid_an5, :);
% visualizeFace(face_node_info, x_to_y, plot_trinode_an4, plot_trinode_an5, 'distort_tri');







