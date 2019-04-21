function dists_f_g = calcFaceToGrainCentroidDist(file, faces)
% ##########################################################################
% * Input
%     - faces_ = [n, 2]  
%           returned by trackFace.m or trackUniqueFace.m
% * Output
%     - dists_f_g = [n, 2]
%           this data is aligned with faces_
% * NOTE 
%     - This function is used in: calcFaceMigSign.m
% ##########################################################################
% ----------------------- load debug data -----------------------
% file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
% % file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/181107.mat', 'tracked_uniqueface_an4');
% file = file_an5;
% faces = tracked_uniqueface_an5;
% clear tracked_uniqueface_an4
% ---------------------------------------------------------------
% ##### Read and Clean Data #####
face_label = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
tri_node = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
grain_centroid = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'))';

grain_centroid(1,:) = [];
tri_node = tri_node(all(face_label > 0, 2), :);
face_label = face_label(all(face_label > 0, 2), :);
face_label = sort(face_label, 2);

% ##### Sort grain_id in faces #####
% """
% Note in this function, grain_id in faces are not necessarily sorted becaseu we want to maintain correspondence in calcFaceMigSign.m
% Need to swap it then swap results back after calculation.
% """
mask_swap = (faces(:,1) > faces(:,2));
faces_sorted = sort(faces, 2);

% i = randi(7004, 1);
% ##### Find face_centroid #####
face_centroid = zeros(size(faces_sorted, 1), 3);
for i = 1:size(faces, 1)
    mask = (face_label(:,1) == faces_sorted(i, 1) & face_label(:,2) == faces_sorted(i, 2));
    nodes = unique(tri_node(mask, :));
    face_centroid(i, :) = sum(node_coord(nodes, :))/length(nodes);
end

% ##### Find dist(face_centroid, grain_centroid) #####
dists_f_g = zeros(size(faces_sorted, 1), 2);
for i = 1:size(faces, 1)
    for j = 1:2
        dists_f_g(i, j) = norm(face_centroid(i, :) - grain_centroid(faces_sorted(i, j), :));
    end
end

% ##### Swap Results for grain_A and grain_B in dists_f_g #####
for i = 1:size(faces, 1)
    if mask_swap(i)
        tmp = dists_f_g(i, 1);
        dists_f_g(i, 1) = dists_f_g(i, 2);
        dists_f_g(i, 2) = tmp;
    end
end

return

% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % """
% % - uncommment line39 and comment all for i = 1:size(faces, 1)
% % - plot grain nodes, and centroids of grains, and centroids of grain faces
% % to see if calculation right. 
% % """
% mask_left_grain_tris = any(face_label==faces_sorted(i, 1), 2);
% mask_right_grain_tris = any(face_label==faces_sorted(i, 2), 2);
% 
% left_grain_nodes = unique(tri_node(mask_left_grain_tris, :));
% right_grain_nodes = unique(tri_node(mask_right_grain_tris, :));
% 
% left_grain_centroid = grain_centroid(faces_sorted(i, 1), :);
% right_grain_centroid = grain_centroid(faces_sorted(i, 2), :);
% f_centroid = face_centroid(i, :);
% 
% hold on
% rotate3d on
% scatter3(node_coord(left_grain_nodes,1), node_coord(left_grain_nodes,2), node_coord(left_grain_nodes,3), ...
%     5, 'o', 'filled', 'MarkerEdgeColor', [0.8, 0.8, 0.8], 'MarkerFaceColor', [0.8, 0.8, 0.8])
% scatter3(node_coord(right_grain_nodes,1), node_coord(right_grain_nodes,2), node_coord(right_grain_nodes,3), ...
%     5, 'x', 'filled', 'MarkerEdgeColor', [0.8, 0.8, 0.8], 'MarkerFaceColor', [0.8, 0.8, 0.8])
% 
% scatter3(left_grain_centroid(1), left_grain_centroid(2), left_grain_centroid(3), 100, 's', 'filled', 'b')
% scatter3(right_grain_centroid(1), right_grain_centroid(2), right_grain_centroid(3), 100, 's', 'filled', 'r')
% scatter3(f_centroid(1), f_centroid(2), f_centroid(3), 100, 's', 'filled', 'g')
% 
% disp(dists_f_g(i, :))
% 
% f_name = ['an5_calcCentroids_pair_', num2str(i)];
% print(f_name, '-dtiff', '-r300')
% 
% hold off
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







