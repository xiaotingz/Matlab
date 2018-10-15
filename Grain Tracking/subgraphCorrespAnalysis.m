%% ############################################################## 
% Contents
% - Number of Pieces of the Subgraphs
% - Trial: Project Nodes Onto Basis Plane
% - Solve Subgraph Corresp, One Complete Plane
% ############################################################## 
%% ############################### Number of Pieces of the Subgraphs ############################### 
load('180822');
% load('180828_piecewise_face')

% ----- find the piece wise faces and their number of pieces -----
face_piecewise = [];
num_pieces = [];
for i = 1:length(num_pieces)
    idx = face_piecewise(i);
    obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
    
    mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
    mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));

    face_tri_nodeid_an4 = tri_node_an4(mask_objface_an4, :);
    face_unique_nodeid_an4 = unique(face_tri_nodeid_an4);
    face_tri_nodeid_an5 = tri_node_an5(mask_objface_an5, :);
    face_unique_nodeid_an5 = unique(face_tri_nodeid_an5);
    [subgraph_an4, subgraph_an5] = findSubgraph(face_unique_nodeid_an4, face_unique_nodeid_an5, face_tri_nodeid_an4, face_tri_nodeid_an5);
    
    m = length(unique(subgraph_1));
    n = length(unique(subgraph_2));
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



%% ############################### Trial: Project Nodes Onto Basis Plane ###############################
% ----- load data -----
facelabel_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
facelabel_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_normal_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
tri_normal_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
    % -- NOTE triNodes are indexes starting from zero --
tri_node_an4 = 1 + double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
tri_node_an5 = 1 + double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
node_coord_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
node_type_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/VertexData/NodeType'))';
node_type_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/VertexData/NodeType'))';
	% -- filter bad data --
mask_an4 = all(facelabel_an4 > 0, 2);
facelabel_an4 = facelabel_an4(mask_an4, :);
tri_node_an4 = tri_node_an4(mask_an4, :);
tri_normal_an4 = tri_normal_an4(mask_an4, :);
mask_an5 = all(facelabel_an5 > 0, 2);
facelabel_an5 = facelabel_an5(mask_an5, :);
tri_node_an5 = tri_node_an5(mask_an5, :);
tri_normal_an5 = tri_normal_an5(mask_an5, :);
color1 = [0, 0.4470, 0.7410];
color2 = [0.9290, 0.6940, 0.1250];
load('180927_tmp.mat')

%% ----- get the mesh nodes -----
idx = 1;
obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);

mask_an4_1 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
mask_an4_2 = (facelabel_an4(:,1) == obj_facelabel_an4(2) & facelabel_an4(:,2) == obj_facelabel_an4(1));
mask_an4 = (mask_an4_1 | mask_an4_2);
mask_an5_1 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));
mask_an5_2 = (facelabel_an5(:,1) == obj_facelabel_an5(2) & facelabel_an5(:,2) == obj_facelabel_an5(1));
mask_an5 = (mask_an5_1 | mask_an5_2);
node_id_an4 = unique(tri_node_an4(mask_an4, :));
node_id_an5 = unique(tri_node_an5(mask_an5, :));

% ----- calc basis plane normal -----
normal_avg_an4 = sum(tri_normal_an4(mask_an4_1, :)) - sum(tri_normal_an4(mask_an4_2, :));
normal_avg_an4 = normal_avg_an4/norm(normal_avg_an4);
normal_avg_an5 = sum(tri_normal_an5(mask_an5_1, :)) - sum(tri_normal_an5(mask_an5_2, :));
normal_avg_an5 = normal_avg_an5/norm(normal_avg_an5);
normal_bp = normal_avg_an4;

% ----- calc basis plane origion from plane centroid -----
nodes_an4 = node_coord_an4(node_id_an4, :);
nodes_an5 = node_coord_an5(node_id_an5, :);
centroid_bp = sum(nodes_an4)/length(node_id_an4);
% quiver3(centroid_bp(1), centroid_bp(2), centroid_bp(3), normal_bp(1), normal_bp(2), normal_bp(3), 5, 'color', [0.4660, 0.6740, 0.1880],  'LineWidth', 3, 'MaxHeadSize', 3);

% ----- project the mesh nodes onto a 2D plane -----
% """
% See project3DPointsTo2DPlane.m for description of native2d_x
% """
[nodes_proj2d_an4, native2d_x] = project3DPointsTo2DPlane(nodes_an4, centroid_bp, normal_bp, 'bp');
[nodes_proj2d_an5, ~] = project3DPointsTo2DPlane(nodes_an5, centroid_bp, normal_bp, native2d_x);

% ##### Boundary by Matlab & by TLNs  #####
% ----- boundary by Matlab -----
bound_an4 = boundary(nodes_proj2d_an4(:,1), nodes_proj2d_an4(:,2));
bound_an5 = boundary(nodes_proj2d_an5(:,1), nodes_proj2d_an5(:,2));

% ----- boundary by TLNs -----
tln_id_an4 = node_id_an4(node_type_an4(node_id_an4) >= 3);
tln_id_an5 = node_id_an5(node_type_an5(node_id_an5) >= 3);
tln_an4 = node_coord_an4(tln_id_an4, :);    
tln_an5 = node_coord_an5(tln_id_an5, :);

[tln_proj2d_an4, ~] = project3DPointsTo2DPlane(tln_an4, centroid_bp, normal_bp, native2d_x);
[tln_proj2d_an5, ~] = project3DPointsTo2DPlane(tln_an5, centroid_bp, normal_bp, native2d_x);

% ##### Projection Score #####
% ----- polygon and intersection -----
pgon_an4 = polyshape(nodes_proj2d_an4(bound_an4, 1), nodes_proj2d_an4(bound_an4, 2));
pgon_an5 = polyshape(nodes_proj2d_an5(bound_an5, 1), nodes_proj2d_an5(bound_an5, 2));
pgon_inter = intersect(pgon_an4,pgon_an5);

% ----- normal direction difference -----
% """
% Notice switching symmetry, max(normal_diff)=90°
% """
norm_diff = abs(atand(norm(cross(normal_avg_an4, normal_avg_an5))/dot(normal_avg_an4, normal_avg_an5)));

% ----- overlap score -----
% """
% The logic is that overlapping area and consistent normal are equally
% important. This logic maybe questionable.
% """
score = 0.5 * (area(pgon_inter)/area(pgon_an5)) + 0.5 * (1 - norm_diff/90);


% ##### Visualization of Projection #####
% ----- 3D view -----
% figure
% x_to_y = X_to_Y{idx};
% face_node_info = getSingleFaceNodes(obj_facelabel_an4, obj_facelabel_an5);
% visualizeFace(face_node_info, x_to_y)

% ----- 2D projected points -----
% figure
% scatter(nodes_proj2d_an4(:,1), nodes_proj2d_an4(:,2), 'MarkerEdgeColor', color1);
% hold on;
% scatter(nodes_proj2d_an5(:,1), nodes_proj2d_an5(:,2), 'MarkerEdgeColor', color2);
% plot(nodes_proj2d_an4(bound_an4, 1), nodes_proj2d_an4(bound_an4, 2));
% plot(nodes_proj2d_an5(bound_an5, 1), nodes_proj2d_an5(bound_an5, 2));
% scatter(tln_proj2d_an4(:,1), tln_proj2d_an4(:,2), 'MarkerEdgeColor', color1, 'MarkerFaceColor', color1);
% scatter(tln_proj2d_an5(:,1), tln_proj2d_an5(:,2), 'MarkerEdgeColor', color2, 'MarkerFaceColor', color2);
% print('pair_1_2Dproj','-dpng','-r300')

% ----- 2D polygons -----
% figure
% plot(pgon_an4, 'FaceColor', color1)
% hold on
% plot(pgon_an5, 'FaceColor', color2)
% plot(pgon_inter, 'FaceColor', [0.4660, 0.6740, 0.1880])
% legend({'face\_an4, area=124', 'face\_an5, area=68', 'intersection, area=113'}, 'FontSize',14, 'Location', 'SouthEast')
% print('pair_1_2Dpolygon','-dpng','-r300')

% ----- 3D avg norm -----
% quiver3(centroid_bp(1), centroid_bp(2), centroid_bp(3), normal_avg_an4(1), normal_avg_an4(2), normal_avg_an4(3), 5, 'LineWidth', 3, 'MaxHeadSize', 3);
% quiver3(centroid_bp(1), centroid_bp(2), centroid_bp(3), normal_avg_an5(1), normal_avg_an5(2), normal_avg_an5(3), 5, 'LineWidth', 3, 'MaxHeadSize', 3);
% daspect([1 1 1])


%% ############################### Solve Subgraph Corresp, One Complete Plane ###############################
% ----- load data -----












