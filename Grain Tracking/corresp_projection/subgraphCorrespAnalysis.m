%% ############################################################## 
% Contents
% - Number of Pieces of the Subgraphs
% - Trial: Project Nodes Onto Basis Plane
% - Solve Subgraph Corresp, One Complete Plane
% ############################################################## 
%% ############################### Pieces of the Subgraphs ############################### 
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/trials/180822.mat');
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/trials/180822_FaceCorresp.mat');
% load('180828_piecewise_face')

% ----- find the piecewise faces and their number of pieces -----
face_piecewise = [];
num_pieces = [];
for i = 1:length(tracked_uniqueface_an4)
    obj_facelabel_an4 = tracked_uniqueface_an4(i, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(i, :);
    
    mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
    mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));

    face_tri_nodeid_an4 = tri_node_an4(mask_objface_an4, :);
    face_unique_nodeid_an4 = unique(face_tri_nodeid_an4);
    face_tri_nodeid_an5 = tri_node_an5(mask_objface_an5, :);
    face_unique_nodeid_an5 = unique(face_tri_nodeid_an5);
    [subgraph_an4, subgraph_an5] = findSubgraph(face_unique_nodeid_an4, face_unique_nodeid_an5, face_tri_nodeid_an4, face_tri_nodeid_an5);
    
    m = length(unique(subgraph_an4));
    n = length(unique(subgraph_an5));
    if m > 1 || n > 1
        face_piecewise = [face_piecewise; i];
        num_pieces = [num_pieces; max(subgraph_an4), max(subgraph_an5)];
    end
end

all_multipiece = face_piecewise(all(num_pieces>1, 2));
tmp = num_pieces(all(num_pieces>1, 2), :);
all_multipiece_asym = all_multipiece(tmp(:,1) ~= tmp(:,2));
all_multipiece_sym = all_multipiece(tmp(:,1) == tmp(:,2));

%% ----- the many-piece faces, are all pieces small -----
amp_asym_size = ones(length(all_multipiece_asym), 2);
for i = 1:length(all_multipiece_asym)
    idx = all_multipiece_asym(i);
    obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
    
    mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
    mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));

    face_tri_nodeid_an4 = tri_node_an4(mask_objface_an4, :);
    face_unique_nodeid_an4 = unique(face_tri_nodeid_an4);
    face_tri_nodeid_an5 = tri_node_an5(mask_objface_an5, :);
    face_unique_nodeid_an5 = unique(face_tri_nodeid_an5);
    [subgraph_an4, subgraph_an5] = findSubgraph(face_unique_nodeid_an4, face_unique_nodeid_an5, face_tri_nodeid_an4, face_tri_nodeid_an5);
    
    cnt_subgraph_an4 = histcounts(subgraph_an4);
    cnt_subgraph_an5 = histcounts(subgraph_an5);
    
    amp_asym_size(i, 1) = max(cnt_subgraph_an4);
    amp_asym_size(i, 2) = max(cnt_subgraph_an5);
end

sum(any(amp_asym_size < 10, 2))





%% ############################### Check, Visualize Subgraph ###############################
x_to_y = X_to_Y{idx};
face_node_info = getSingleFaceNodes(file_an4, file_an5, tracked_uniqueface_an4(idx,:), tracked_uniqueface_an5(idx,:));
visualizeFace(face_node_info, x_to_y)
hold on


% ----- subgraph corresps -----
coord = facenode_coord_an4;
id = subgraph_id_an4;
obj = 1;
scatter3(coord(id==obj, 1), coord(id==obj, 2), coord(id==obj, 3), 'MarkerFaceColor', 'r')
coord = facenode_coord_an5;
id = subgraph_id_an5;
obj = 1;
scatter3(coord(id==obj, 1), coord(id==obj, 2), coord(id==obj, 3), 'MarkerFaceColor', 'r')

coord = facenode_coord_an4;
id = subgraph_id_an4;
obj = 2;
scatter3(coord(id==obj, 1), coord(id==obj, 2), coord(id==obj, 3), 'MarkerFaceColor', 'g')
coord = facenode_coord_an5;
id = subgraph_id_an5;
obj = 1;
scatter3(coord(id==obj, 1), coord(id==obj, 2), coord(id==obj, 3), 'MarkerFaceColor', 'g')

daspect([1 1 1])

%% ##################### Detail Subgraphs, Normal Direction #####################
colors = get(gca,'colororder');
mask_an4 = mask_tri_subgraph_an4==1;
mask_an5 = mask_tri_subgraph_an5==2;

figure
trisurf(facetri_nodeid_an4(mask_an4, :), node_coord_an4(:,1), node_coord_an4(:,2), node_coord_an4(:,3),'Facecolor',colors(1, :), 'Facealpha', 0.3, 'edgealpha', 0.3);
hold on
pos_an4 = node_coord_an4(facetri_nodeid_an4(mask_an4, 1), :);
norm_an4 = facetri_normal_an4(mask_an4, :);
quiver3(pos_an4(:,1), pos_an4(:,2), pos_an4(:,3), norm_an4(:,1), norm_an4(:,2), norm_an4(:,3), 'color',colors(1, :))
trisurf(facetri_nodeid_an5(mask_an5, :), node_coord_an5(:,1), node_coord_an5(:,2), node_coord_an5(:,3),'Facecolor',colors(3, :), 'Facealpha', 0.3, 'edgealpha', 0.3);
hold on
pos_an5 = node_coord_an5(facetri_nodeid_an5(mask_an5, 1), :);
norm_an5 = facetri_normal_an5(mask_an5, :);
quiver3(pos_an5(:,1), pos_an5(:,2), pos_an5(:,3), norm_an5(:,1), norm_an5(:,2), norm_an5(:,3), 'color',colors(3, :))

% a = sum(facetri_normal_an5(mask_tri_subgraph_an5==2, :))/2;
% b = sum(facetri_normal_an4(mask_tri_subgraph_an4==1, :))/sum(mask_tri_subgraph_an4==1);
% quiver3(pos_an4(1,1), pos_an4(1,2), pos_an4(1,3),a(1),a(2),a(3), 3)
% quiver3(pos_an5(1,1), pos_an5(1,2), pos_an5(1,3),b(1),b(2),b(3), 3)

daspect([1 1 1])
rotate3d on

%% ##################### 2D nodes & boundary #####################
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/trials/180822.mat');
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/trials/180822_FaceCorresp.mat');
file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/D3Ds/An5new6_smooth.dream3d';
idx = 6719;

x_to_y = X_to_Y{idx};
face_node_info = getSingleFaceNodes(file_an4, file_an5, tracked_uniqueface_an4(idx,:), tracked_uniqueface_an5(idx,:));
visualizeFace(face_node_info, x_to_y)
set(gca,'visible','off')


%%
% ------------ Read info ------------
facelabel_an4 = double(h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
facelabel_an5 = double(h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
normal_an4 = double(h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
normal_an5 = double(h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
facelabel_an4 = sort(facelabel_an4, 2);
facelabel_an5 = sort(facelabel_an5, 2);

mask_obj_face_an4 = ( facelabel_an4(:,1) == tracked_uniqueface_an4(idx,1) & ...
    facelabel_an4(:,2) == tracked_uniqueface_an4(idx,2) ); 
mask_obj_face_an5 = ( facelabel_an5(:,1) == tracked_uniqueface_an5(idx,1) & ...
    facelabel_an5(:,2) == tracked_uniqueface_an5(idx,2) ); 

face_normals_an4 = normal_an4(mask_obj_face_an4, :);
face_normals_an5 = normal_an5(mask_obj_face_an5, :);

% ------------ Origin of the base plane ------------
face_nodes_an4 = face_node_info{4, 1};
face_nodes_an5 = face_node_info{4, 2};
face_node_types_an4 = face_node_info{3, 1};
face_node_types_an5 = face_node_info{3, 2};
bp_origin = sum([face_nodes_an4; face_nodes_an5]) / size([face_nodes_an4; face_nodes_an5], 1);

% ------------ Normal of the base plane ------------
bp_normal = sum([face_normals_an4; face_normals_an5]) / size([face_normals_an4; face_normals_an5], 1);

% ------------ Get 2D points ------------
[points_2d_an4, native2d_x_an4] = project3DPointsTo2DPlane(face_nodes_an4, bp_origin, bp_normal, 'bp');
[points_2d_an5, native2d_x_an5] = project3DPointsTo2DPlane(face_nodes_an5, bp_origin, bp_normal, 'bp');


%% ------------ Fig b ------------
colors = get(gca,'colororder');

figure()
bound_cp_an4 = boundary(points_2d_an4(:,1), points_2d_an4(:,2));
pgon_cp_an4 = polyshape(points_2d_an4(bound_cp_an4, 1), points_2d_an4(bound_cp_an4, 2));
plot(pgon_cp_an4, 'FaceAlpha', 0, 'EdgeColor', colors(1,:), 'LineWidth', 2)
hold on
scatter(points_2d_an4(:,1), points_2d_an4(:,2), 'MarkerFaceColor', [1,1,1], 'MarkerEdgeColor', colors(1,:));
scatter(points_2d_an4(face_node_types_an4>=3,1), points_2d_an4(face_node_types_an4>=3,2), 'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:));

bound_cp_an5 = boundary(points_2d_an5(:,1), points_2d_an5(:,2));
pgon_cp_an5 = polyshape(points_2d_an5(bound_cp_an5, 1), points_2d_an5(bound_cp_an5, 2));
plot(pgon_cp_an5, 'FaceAlpha', 0, 'EdgeColor', colors(3,:), 'LineWidth', 2)
scatter(points_2d_an5(:,1), points_2d_an5(:,2), 'MarkerFaceColor', [1,1,1], 'MarkerEdgeColor', colors(3,:));
scatter(points_2d_an5(face_node_types_an5>=3,1), points_2d_an5(face_node_types_an5>=3,2), 'MarkerFaceColor', colors(3,:), 'MarkerEdgeColor', colors(3,:));

daspect([1 1 1])
set(gca,'visible','off')

%%
%  ------------ Fig c ------------
poly_intersect = intersect(pgon_cp_an4, pgon_cp_an5);
mask_good_nodes_an4 = inpolygon(points_2d_an4(:, 1), points_2d_an4(:, 2), ...
    poly_intersect.Vertices(:, 1), poly_intersect.Vertices(:,2));
mask_good_nodes_an5 = inpolygon(points_2d_an5(:, 1), points_2d_an5(:, 2), ...
    poly_intersect.Vertices(:, 1), poly_intersect.Vertices(:,2));

figure()
plot(pgon_cp_an4, 'FaceAlpha', 0, 'EdgeColor', colors(1, :), 'LineWidth', 2)
hold on
scatter(points_2d_an4(mask_good_nodes_an4,1),points_2d_an4(mask_good_nodes_an4,2), 'p', ...
    'MarkerFaceColor', colors(1,:), 'MarkerEdgeColor', colors(1,:))
plot(pgon_cp_an5, 'FaceAlpha', 0, 'EdgeColor', colors(3, :), 'LineWidth', 2)
scatter(points_2d_an5(mask_good_nodes_an5,1),points_2d_an5(mask_good_nodes_an5,2), 'p', ...
    'MarkerFaceColor', colors(3,:), 'MarkerEdgeColor', colors(3,:))

daspect([1 1 1])
set(gca,'visible','off')











