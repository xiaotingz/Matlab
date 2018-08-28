% function migration = calcMigrationFromCorresp(file_an4, obj_facelabel_an4, obj_facelabel_an5, x_to_y)
% ############################################################################
% * Input 
%     - obj_facelabel_ = [2,1], facelabel of the objective face. 
%         The full facelabel_list is returned by TrackUniqueFace.m
%     - x_to_y = [n, 1], the corresponden ce
%         Returned by main_TrackNodes.m or solveNodeCorresp.m
%  * Output
%     - migration = migration distance of this grain face
% NOTES
%     - There are three models:
%           direct migration distance, projected along local triangle normal, and truncated prism
%     - The basic idea behind the truncated prism model is migration_distance = swept_volume/surface_area
%           The geometric model is tri-prism with oblique bottom and top surfaces
% ############################################################################
% ------------------ load data for debug --------------------
load('180822_FaceCorresp');
load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822.mat');
file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
idx = 654;
obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
x_to_y = X_to_Y{idx};
% -----------------------------------------------------------
% load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822.mat');
face_tri_normal_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FaceNormals')';
mask_inner = all(facelabel_an4 > 0, 2);
face_tri_normal_an4 = face_tri_normal_an4(mask_inner, :);


% ##### Get Id and Doord of the Objective Triangles #####
mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));
face_tri_nodeid_an4 = tri_node_an4(mask_objface_an4, :);
face_tri_nodeid_an5 = tri_node_an5(mask_objface_an5, :);
face_tri_normal_an4 = face_tri_normal_an4(mask_objface_an4, :);
face_unique_nodeid_an4 = unique(face_tri_nodeid_an4);
face_unique_nodeid_an5 = unique(face_tri_nodeid_an5);


% ##### Construct Triangle on face_an5 Following the Correspondence #####
map = containers.Map(face_unique_nodeid_an4, face_unique_nodeid_an5(x_to_y));
face_tri_5from4 = zeros(size(face_tri_nodeid_an4));
for i = 1:size(face_tri_5from4, 1)
    for j = 1:3
        face_tri_5from4(i,j) = map(face_tri_nodeid_an4(i,j));
    end
    
% ----- If two points correspond to one point, then can't form a prism -----
    if length(unique(face_tri_5from4(i,:))) < 3
        face_tri_5from4(i,:) = NaN;
    end
 
end
mask_bad_triprism = any(isnan(face_tri_5from4), 2);
face_tri_nodeid_an4(mask_bad_triprism, :) = [];
face_tri_normal_an4(mask_bad_triprism, :) = [];
face_tri_5from4(mask_bad_triprism, :) = [];

% ##### Go Through Triangle Pairs, Calculate Height for each Tri-prims from Triangle Pair #####
h = zeros(size(face_tri_5from4, 1), 1);
for i = 1:length(h)
    tri1_nodes = node_coord_an4(face_tri_nodeid_an4(i, :), :);
    tri2_nodes = node_coord_an5(face_tri_5from4(i, :), :);
    DT = delaunayTriangulation([tri1_nodes; tri2_nodes]);
    [C, v] = convexHull(DT);
    area1 = calcArea(tri1_nodes);
    area2 = calcArea(tri2_nodes);
    h(i) = heightUsingTruncatedConeModel(area1, area2, v);
end
tri_migration_prismheight = sum(h)/length(h);


% ##### Direct Migration Distance #####
face_unique_nodeid_5from4 = zeros(size(face_unique_nodeid_an4, 1), 1);
for i = 1:length(face_unique_nodeid_an4)
    face_unique_nodeid_5from4(i) = map(face_unique_nodeid_an4(i));
end
unique_nodecoord_an4 = node_coord_an4(face_unique_nodeid_an4, :);
unique_nodecoord_an5 = node_coord_an5(face_unique_nodeid_5from4, :);
d_direct = vecnorm(unique_nodecoord_an5 - unique_nodecoord_an4, 2, 2);
tri_migration_direct = sum(d_direct)/length(d_direct);


% ##### Migration by Projection #####
% ----- a node corresp along normals of all triangles that it belongs to -----
face_tri_normal_an4 = repmat(face_tri_normal_an4, 3, 1);
d_tmp = node_coord_an4(face_tri_nodeid_an4(:), :) - node_coord_an5(face_tri_5from4(:), :);
d_proj = abs(dot(d_tmp, face_tri_normal_an4, 2));
d_proj = sum(reshape(d_proj, [], 3), 2) / 3;
tri_migration_projection = sum(d_proj)/length(d_proj);


% 
% trisurf(C,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3));
% rotate3d on
% 
% scatter3(tri1_nodes(:,1), tri1_nodes(:,2), tri1_nodes(:,3))
% hold on
% scatter3(tri2_nodes(:,1), tri2_nodes(:,2), tri2_nodes(:,3))
% 





