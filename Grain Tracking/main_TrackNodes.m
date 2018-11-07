% % ############################################################################
% % Solve Matching
% % * Notes
% %     - In this file, the coordinates on face are not stored because the
% %     memory cost is too expensive. However, the face coordinates can be obtained by
% %     running getSingleFaceNodes.m
%       - Three strategies
%           (1) Solve corresp between all nodes directly
%           (2) First find the disconnected subgraphs, solve corresp between subgraphs then solve nodes corresp
%           (3) Sample a subset of points instead of using all points
%           The (2) and (3) strategy can be combined
% % ############################################################################
file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
load('look_up_table_an4_an5.mat')
facelabel_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
facelabel_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_normal_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
tri_normal_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
    % -- NOTE triNodes are indexes starting from zero 
tri_node_an4 = 1 + double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
tri_node_an5 = 1 + double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
node_coord_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
mask_inner_an4 = all(facelabel_an4 > 0, 2);
mask_inner_an5 = all(facelabel_an5 > 0, 2);
facelabel_an4 = facelabel_an4(mask_inner_an4, :);
facelabel_an5 = facelabel_an5(mask_inner_an5, :);
facelabel_an4 = sort(facelabel_an4, 2);
facelabel_an5 = sort(facelabel_an5, 2);
tri_node_an4 = tri_node_an4(mask_inner_an4, :);
tri_node_an5 = tri_node_an5(mask_inner_an5, :);
tri_normal_an4 = tri_normal_an4(mask_inner_an4, :);
tri_normal_an5 = tri_normal_an5(mask_inner_an5, :);

% --- tracked_unique_face = trackUniqueFace(file_an4, file_an5, look_up_table, complete) ---
[tracked_uniqueface_an4, tracked_uniqueface_an5] = trackUniqueFace(file_an4, file_an5, look_up_table, 1);
n = length(tracked_uniquean4);
X_to_Y = cell(n,1);
Y_to_X = cell(n,1);

% %% ############################################################################
% (1) Dicret Corresp
% % ############################################################################
% %% ##### (1.1) Corresp of The Small Faces #####
% large_face_id = [];
% for i = 1:length(tracked_uniqueface_an4)
%     obj_facelabel_an4 = tracked_uniqueface_an4(i, :);
%     obj_facelabel_an5 = tracked_uniqueface_an5(i, :);
%     
%     % ##### Get Objective Triangles on the Objective Face #####
%     mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
%     mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));
% 
%     % ##### Get Id and Coord of the Objective Triangles #####
%     face_tri_nodeid_an4 = tri_node_an4(mask_objface_an4, :);
%     face_unique_nodeid_an4 = unique(face_tri_nodeid_an4);
%     face_node_coord_an4 = face_node_coord_an4(face_unique_nodeid_an4,:);
%     face_tri_nodeid_an5 = tri_node_an5(mask_objface_an5, :);
%     face_unique_nodeid_an5 = unique(face_tri_nodeid_an5);
%     face_node_coord_an5 = face_node_coord_an5(face_unique_nodeid_an5,:);
%     
%     if (length(face_unique_nodeid_an4) < 1000) && (length(face_unique_nodeid_an5) < 1000)
%         [x_to_y, y_to_x] = solveNodeCorresp(face_node_coord_an4, face_node_coord_an5);
%         X_to_Y{i} = int32(x_to_y);
%         Y_to_X{i} = int32(y_to_x)';
%     else 
%         large_face_id = [large_face_id, i];
%         X_to_Y{i} = NaN;
%         Y_to_X{i} = NaN;
%     end
%     disp(['small faces: ', num2str(i), '/7004']);
% end
% save('180820_smallFaces.mat', 'X_to_Y', 'Y_to_X', 'large_id');
% 
% 
% %% ##### (1.2) Corresp of The Large Faces #####
% for i = 1:length(large_face_id)
%     idx = large_face_id(i);
%     obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
%     obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
%     
%     % ##### Get Objective Triangles on the Objective Face #####
%     mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
%     mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));
% 
%     % ##### Get Id and Coord of the Objective Triangles #####
%     face_tri_nodeid_an4 = tri_node_an4(mask_objface_an4, :);
%     face_unique_nodeid_an4 = unique(face_tri_nodeid_an4);
%     face_node_coord_an4 = face_node_coord_an4(face_unique_nodeid_an4,:);
%     face_tri_nodeid_an5 = tri_node_an5(mask_objface_an5, :);
%     face_unique_nodeid_an5 = unique(face_tri_nodeid_an5);
%     face_node_coord_an5 = face_node_coord_an5(face_unique_nodeid_an5,:);
%     
%     [x_to_y, y_to_x] = solveNodeCorresp(face_node_coord_an4, face_node_coord_an5);
%     X_to_Y_large{i} = int32(x_to_y);
%     Y_to_X_large{i} = int32(y_to_x);
%     
%     disp(['large faces: ', num2str(i), '/', num2str(large_sizee)]);
% end
% save('180820_largeFaces.mat', 'X_to_Y_large', 'Y_to_X_large');
% 

%% ############################################################################
% (2) First Find Disconnected Subgraphs 
%      - Not every subgraph will find its corresp. 
%           length(subgraph_pair) = min(length(subgraph_an4), length(subgraph_an5))
% % ############################################################################
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/180822.mat', 'tracked_uniqueface_an4', 'tracked_uniqueface_an5');
load('180822_FaceCorresp.mat', 'X_to_Y');
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/180828_piecewise_face.mat', 'face_piecewise');
file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');

% ##### Data Preparation #####
facelabel_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
facelabel_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_normal_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
tri_normal_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
tri_area_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'))';
tri_area_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'))';
    % -- NOTE triNodes are indexes starting from zero --
tri_node_an4 = 1 + double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
tri_node_an5 = 1 + double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
node_coord_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
% ---- filter bad data -----
mask_an4 = all(facelabel_an4 > 0, 2);
mask_an5 = all(facelabel_an5 > 0, 2);
facelabel_an4 = facelabel_an4(mask_an4, :);
facelabel_an5 = facelabel_an5(mask_an5, :);
tri_node_an4 = tri_node_an4(mask_an4, :);
tri_node_an5 = tri_node_an5(mask_an5, :);
tri_normal_an4 = tri_normal_an4(mask_an4, :);
tri_normal_an5 = tri_normal_an5(mask_an5, :);
tri_area_an4 = tri_area_an4(mask_an4, :);
tri_area_an5 = tri_area_an5(mask_an5, :);
%%

% face_piecewise = [];
% idx_face_piecewise = 1;
% X_to_Y = cell(n,10);
% Y_to_X = cell(n,10);
for i = 6879
    obj_facelabel_an4 = tracked_uniqueface_an4(i, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(i, :);

    % ##### Get Objective Triangles on the Objective Face #####
    mask_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2) |...
        facelabel_an4(:,1) == obj_facelabel_an4(2) & facelabel_an4(:,2) == obj_facelabel_an4(1));
    mask_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2) |...
        facelabel_an5(:,1) == obj_facelabel_an5(2) & facelabel_an5(:,2) == obj_facelabel_an5(1));
    mask_reverse_an4 = facelabel_an4(mask_an4,1) > facelabel_an4(mask_an4,2);
    mask_reverse_an5 = facelabel_an5(mask_an5,1) > facelabel_an5(mask_an5,2);
    

    % ##### Get Id and Coord of the Objective Triangles #####
    facetri_nodeid_an4 = tri_node_an4(mask_an4, :);
    facenode_id_an4 = unique(facetri_nodeid_an4);
    facenode_coord_an4 = node_coord_an4(facenode_id_an4,:);
    facetri_normal_an4 = tri_normal_an4(mask_an4, :);
    facetri_normal_an4(mask_reverse_an4, :) = - facetri_normal_an4(mask_reverse_an4, :);
    facetri_area_an4 = tri_area_an4(mask_an4, :);
    facetri_nodeid_an5 = tri_node_an5(mask_an5, :);
    facenode_id_an5 = unique(facetri_nodeid_an5);
    facenode_coord_an5 = node_coord_an5(facenode_id_an5,:);
    facetri_normal_an5 = tri_normal_an5(mask_an5, :);
    facetri_normal_an5(mask_reverse_an5, :) = - facetri_normal_an5(mask_reverse_an5, :);
    facetri_area_an5 = tri_area_an5(mask_an5, :);

    % ##### Find Disconnected Subgraphs and Solve Subgraph Correp #####
    [subgraph_id_an4, subgraph_id_an5] = findSubgraph(facenode_id_an4, facenode_id_an5, facetri_nodeid_an4, facetri_nodeid_an5);
    
    % ##### Remove A Tiny Piece Near A Large Face #####
    
    
    % ##### If One-piece #####
    if length(unique(subgraph_id_an4)) == 1 && length(unique(subgraph_id_an5)) == 1
        % ----- Solve node corresp -----
        [x_to_y, y_to_x] = solveNodeCorresp(facenode_coord_an4, facenode_coord_an5);
        X_to_Y{i, 1} = int32(x_to_y);
        Y_to_X{i, 1} = int32(y_to_x)';
    
    % ##### If Multiple-pieces  #####
    else
        % ----- Solve the corresp of subgraphs -----
        subgraph_corresp = solveSubgraphCorrespOverlap(subgraph_id_an4,subgraph_id_an5, ...
            facetri_normal_an4, facetri_normal_an5, facetri_nodeid_an4, facetri_nodeid_an5, facenode_id_an4, facenode_id_an5, ...
            facenode_coord_an4, facenode_coord_an5, facetri_area_an4, facetri_area_an5);
        
        % ----- Because there are 1-N corresps, need to first find the group of subgraphs -----
        useful_subgraph_an4 = unique(subgraph_corresp(:,1));
        useful_subgraph_an5 = unique(subgraph_corresp(:,2));
        [num_group, idx_fewerpiece] = min([length(useful_subgraph_an4), length(useful_subgraph_an5)]);
        for j = 1:num_group
            if idx_fewerpiece == 1
                group_subgraphid_an4 = useful_subgraph_an4(j);
                group_subgraphid_an5 = subgraph_corresp(subgraph_corresp(:,1)==group_subgraphid_an4, 2);
                mask_group_nodes_an4 = (subgraph_id_an4 == group_subgraphid_an4);
                mask_group_nodes_an5 = ismember(subgraph_id_an5, group_subgraphid_an5);
                subgraph_nodeid_an4 = facenode_id_an4(mask_group_nodes_an4, :);
                subgraph_nodeid_an5 = facenode_id_an5(mask_group_nodes_an5, :);
                subgraph_coord_an4 = facenode_coord_an4(mask_group_nodes_an4, :);
                subgraph_coord_an5 = facenode_coord_an5(mask_group_nodes_an5, :);
            else
                group_subgraphid_an5 = useful_subgraph_an5(j);
                group_subgraphid_an4 = subgraph_corresp(subgraph_corresp(:,2)==group_subgraphid_an5, 1);
                mask_group_nodes_an5 = (subgraph_id_an5 == group_subgraphid_an5);
                mask_group_nodes_an4 = ismember(subgraph_id_an4, group_subgraphid_an4);
                subgraph_nodeid_an4 = facenode_id_an4(mask_group_nodes_an4, :);
                subgraph_nodeid_an5 = facenode_id_an5(mask_group_nodes_an5, :);
                subgraph_coord_an4 = facenode_coord_an4(mask_group_nodes_an4, :);
                subgraph_coord_an5 = facenode_coord_an5(mask_group_nodes_an5, :);
            end
            
            % ----- Find node corresp between the subgraph groups -----
%             [x_to_y, y_to_x] = solveNodeCorresp(subgraph_coord_an4, subgraph_coord_an5);
%             X_to_Y{i, j} = int32(x_to_y);
%             Y_to_X{i, j} = int32(y_to_x)';
%             id_X_to_Y{i, j} = int32([subgraph_nodeid_an4(x_to_y), subgraph_nodeid_an5])
%             id_Y_to_X{i, j} = 
        end
        
        % ----- Record id of the piece_wise faces-----
%         face_piecewise = [face_piecewise, i, max(subgraph_id_an4), max(subgraph_id_an5)];
    end
    
end



