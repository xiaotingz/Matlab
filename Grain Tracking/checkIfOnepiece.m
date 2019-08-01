function is_one_piece = checkIfOnepiece(file, faces)
% ############################################################################
% * Input 
%     - faces = [n, 2]
%         Facelabels of the faces of interest
% * Output
%     - is_one_piece = [n, 1]
%         Indicator variable showing if the face contain multiple pieces.
% * Notes
%     - This function is written from /solve_corresp/findSubgraph.m and
%     is used in featurePrepOther.m
% ############################################################################
% ------------------ load data for debug --------------------
% file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
% file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4.dream3d';
% load('look_up_table_an4_an5crop.mat')
% [tracked_uniqueface_an4_inner, ~] = trackUniqueFace(file_an4, file_an5, look_up_table, 'use_inner_faces');
% file = file_an4;
% faces = tracked_uniqueface_an4_inner;
% clear file_an4 file_an5 tracked_uniqueface_an4_inner
% -----------------------------------------------------------

% ##### Load data #####
tri_fls = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels').';
tri_nodes = 1 + h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')';
node_types = h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType')';
mask = all(tri_fls > 0, 2);
tri_fls = tri_fls(mask, :);
tri_nodes = tri_nodes(mask, :);

tri_fls = sort(tri_fls, 2);
faces = sort(faces, 2);

is_one_piece = ones(size(faces, 1), 1);
for i = 1:size(faces, 1)
    % ##### Display calculation process #####
    if rem(i, 1000) == 0
        disp(['processed ', num2str(i), ' / ', num2str(size(faces, 1)), ' faces in checkIfOnepiece.m']);
    end
    
    % ##### Collect triangles and nodes on the face #####
    mask = (tri_fls(:,1) == faces(i, 1) & tri_fls(:,2) == faces(i, 2));
    tri_nodes_on_face = tri_nodes(mask, :);
    unique_nodes_on_face = unique(tri_nodes_on_face);

    % ##### Construct Graph of the Face #####
    % ----- convert nodes to be local, otherwise graph is too big -----
    face_map = containers.Map(unique_nodes_on_face, 1:length(unique_nodes_on_face));
    tri_nodes_local = zeros(size(tri_nodes_on_face));
    for j = 1:size(tri_nodes_on_face, 1)
        for k = 1:3
            tri_nodes_local(j, k) = face_map(tri_nodes_on_face(j, k));
        end
    end

    % ----- get the unique edges (unique connections) -----
    edges = [tri_nodes_local(:,1), tri_nodes_local(:, 2); ...
        tri_nodes_local(:,2), tri_nodes_local(:, 3); ...
        tri_nodes_local(:,3), tri_nodes_local(:, 1)];
    edges = sort(edges, 2);
    edges = unique(edges, 'rows');
    % ----- make adjacency matrix -----
    adj_mat = sparse(zeros(length(unique_nodes_on_face)));
    adj_mat_ind = sub2ind(size(adj_mat), edges(:,1), edges(:,2));
    adj_mat(adj_mat_ind) = 1;
    % ----- make graph and identify subgraph -----
    graph_face = graph(adj_mat, 'upper');
    subgraph_id = conncomp(graph_face)';
    if max(subgraph_id) > 1
        is_one_piece(i) = 0;
    end
    
end

end
% %% ################### Checks ###################
% % ----- Prepare data -----
% node_coords = h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList')';
% feature_face_label = h5read(file, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
% feature_face_label(1,:) = [];
% color1 = [0, 0.4470, 0.7410];
% ids = (1:size(is_one_piece, 1))';
% feature_face_id = (1:size(feature_face_label, 1))';
% 
% %%
% % ----- Chose a face -----
% candidates = ids(boolean(~is_one_piece));
% idx = candidates(randi(length(candidates)));
% 
% % ----- Plot the face -----
% mask = tri_fls(:, 1) == faces(idx, 1) & tri_fls(:, 2) == faces(idx, 2);
% face_nodes =  tri_nodes(mask, :);
% trisurf(face_nodes, node_coords(:,1), node_coords(:,2), node_coords(:,3),'Facecolor',color1, 'Facealpha', 0.3, 'edgealpha', 0.3);
% daspect([1 1 1])
% rotate3d on
% 
% % ----- Display D3D FeatureFaceId -----
% mask_featurefaceid = feature_face_label(:, 1) == faces(idx, 1) & feature_face_label(:, 2) == faces(idx, 2);
% disp(['facelabel = ', mat2str(faces(idx, :)), ';  FeatureFaceID = ', num2str(feature_face_id(mask_featurefaceid))]);



