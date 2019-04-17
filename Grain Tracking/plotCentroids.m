function plotCentroids(file_an4, file_an5, faces_an4, faces_an5, idx)
% ##########################################################################
% * Input
%     - faces_ = [n, 2]  
%           returned by trackFace.m or trackUniqueFace.m
%     - idx
%           index of the face of interest
% * NOTE 
%     - This file is to help debug calcFaceToGrainCentroidDist.m
% ##########################################################################


% ##### Read and Clean Data #####
face_label_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
tri_node_an4 = 1 + double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
grain_centroid_an4 = double(h5read(file_an4,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'))';
face_label_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
tri_node_an5 = 1 + double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
grain_centroid_an5 = double(h5read(file_an5,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'))';

grain_centroid_an4(1,:) = [];
tri_node_an4 = tri_node_an4(all(face_label_an4 > 0, 2), :);
face_label_an4 = face_label_an4(all(face_label_an4 > 0, 2), :);
face_label_an4 = sort(face_label_an4, 2);
grain_centroid_an5(1,:) = [];
tri_node_an5 = tri_node_an5(all(face_label_an5 > 0, 2), :);
face_label_an5 = face_label_an5(all(face_label_an5 > 0, 2), :);
face_label_an5 = sort(face_label_an5, 2);

faces_sorted_an4 = sort(faces_an4, 2);
faces_sorted_an5 = sort(faces_an5, 2);

% ##### Find face centroids #####
mask = (face_label_an4(:,1) == faces_sorted_an4(idx, 1) & face_label_an4(:,2) == faces_sorted_an4(idx, 2));
nodes = unique(tri_node_an4(mask, :));
face_centroid_an4 = sum(node_coord_an4(nodes, :))/length(nodes);
mask = (face_label_an5(:,1) == faces_sorted_an5(idx, 1) & face_label_an5(:,2) == faces_sorted_an5(idx, 2));
nodes = unique(tri_node_an5(mask, :));
face_centroid_an5 = sum(node_coord_an5(nodes, :))/length(nodes);

% ##### Find grain centroids #####
mask_left_grain_tris = any(face_label_an4==faces_sorted_an4(idx, 1), 2);
mask_right_grain_tris = any(face_label_an4==faces_sorted_an4(idx, 2), 2);
left_grain_nodes_an4 = unique(tri_node_an4(mask_left_grain_tris, :));
right_grain_nodes_an4 = unique(tri_node_an4(mask_right_grain_tris, :));
left_grain_centroid_an4 = grain_centroid_an4(faces_an4(idx, 1), :);
right_grain_centroid_an4 = grain_centroid_an4(faces_an4(idx, 2), :);

mask_left_grain_tris = any(face_label_an5==faces_sorted_an5(idx, 1), 2);
mask_right_grain_tris = any(face_label_an5==faces_sorted_an5(idx, 2), 2);
left_grain_nodes_an5 = unique(tri_node_an5(mask_left_grain_tris, :));
right_grain_nodes_an5 = unique(tri_node_an5(mask_right_grain_tris, :));
left_grain_centroid_an5 = grain_centroid_an5(faces_an5(idx, 1), :);
right_grain_centroid_an5 = grain_centroid_an5(faces_an5(idx, 2), :);


% ##### Plot #####
hold on
rotate3d on
scatter3(node_coord_an4(left_grain_nodes_an4,1), node_coord_an4(left_grain_nodes_an4,2), node_coord_an4(left_grain_nodes_an4,3), ...
    5, 'o', 'filled', 'MarkerEdgeColor', [0.8, 0.8, 0.8], 'MarkerFaceColor', [0.8, 0.8, 0.8])
scatter3(node_coord_an4(right_grain_nodes_an4,1), node_coord_an4(right_grain_nodes_an4,2), node_coord_an4(right_grain_nodes_an4,3), ...
    5, 'x', 'filled', 'MarkerEdgeColor', [0.8, 0.8, 0.8], 'MarkerFaceColor', [0.8, 0.8, 0.8])

scatter3(left_grain_centroid_an4(1), left_grain_centroid_an4(2), left_grain_centroid_an4(3), 160, 's', 'filled', 'b')
scatter3(right_grain_centroid_an4(1), right_grain_centroid_an4(2), right_grain_centroid_an4(3), 160, 's', 'filled', 'r')
scatter3(face_centroid_an4(1), face_centroid_an4(2), face_centroid_an4(3), 160, 's', 'filled', 'g')

scatter3(node_coord_an5(left_grain_nodes_an5,1), node_coord_an5(left_grain_nodes_an5,2), node_coord_an5(left_grain_nodes_an5,3), ...
    5, '.', 'MarkerEdgeColor', [0.3, 0.3, 0.3])
scatter3(node_coord_an5(right_grain_nodes_an5,1), node_coord_an5(right_grain_nodes_an5,2), node_coord_an5(right_grain_nodes_an5,3), ...
    5, '.', 'MarkerEdgeColor', [0.3, 0.3, 0.3])

scatter3(left_grain_centroid_an5(1), left_grain_centroid_an5(2), left_grain_centroid_an5(3), 160, 's', 'b')
scatter3(right_grain_centroid_an5(1), right_grain_centroid_an5(2), right_grain_centroid_an5(3), 160, 's', 'r')
scatter3(face_centroid_an5(1), face_centroid_an5(2), face_centroid_an5(3), 160, 's','g')

hold off




