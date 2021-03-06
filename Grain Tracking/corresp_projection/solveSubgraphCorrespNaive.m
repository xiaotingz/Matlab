function subgraph_corresp = solveSubgraphCorrespNaive(subgraph_1, subgraph_2, face_node_coord_1, face_node_coord_2)
% ############################################################################
% * Input 
%     - face_unique_nodeid_ = [n, 1], generation see load data
%     - subgraph_ = [n, 1], containing the subgraph_id for every node in face_unique_nodeid_
%         Given by findSubgraph.m
%     - face_tri_coord_ = [n, 3], generation see load data
%  * Output
%     - subgraph_corresp = [m, 1], m = min(num_subgraph_1, num_subgraph_2)
%  * Notes
%     - Not the full corresp, but only that mapping the fewer to the more, is found.
%     - Use in main_TrackNodes.m, together with findSubgraph.m
% ############################################################################
% % ------------------ load data for debug --------------------
% load('180822_FaceCorresp');
% load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822.mat');
% idx = 5602;
% obj_facelabel_1 = tracked_uniqueface_an4(idx, :);
% obj_facelabel_2 = tracked_uniqueface_an5(idx, :);
% 
% mask_objface_1 = (facelabel_an4(:,1) == obj_facelabel_1(1) & facelabel_an4(:,2) == obj_facelabel_1(2));
% mask_objface_2 = (facelabel_an5(:,1) == obj_facelabel_2(1) & facelabel_an5(:,2) == obj_facelabel_2(2));
% 
% face_tri_nodeid_1 = tri_node_an4(mask_objface_1, :);
% face_unique_nodeid_1 = unique(face_tri_nodeid_1);
% face_node_coord_1 = node_coord_an4(face_unique_nodeid_1,:);
% face_tri_nodeid_2 = tri_node_an5(mask_objface_2, :);
% face_unique_nodeid_2 = unique(face_tri_nodeid_2);
% face_node_coord_2 = node_coord_an5(face_unique_nodeid_2,:);
% [subgraph_1, subgraph_2] = findSubgraph(face_unique_nodeid_1, face_unique_nodeid_2, face_tri_nodeid_1, face_tri_nodeid_2);
% % % -----------------------------------------------------------

num_subgraph_1 = length(unique(subgraph_1));
num_subgraph_2 = length(unique(subgraph_2));
subgraph_dist = zeros(num_subgraph_1, num_subgraph_2);
for i = 1:num_subgraph_1
    for j = 1:num_subgraph_2
        mask_subgraph_1 = (subgraph_1 == i);
        mask_subgraph_2 = (subgraph_2 == j);
        subgraph_nodes_1 = face_node_coord_1(mask_subgraph_1, :);
        subgraph_nodes_2 = face_node_coord_2(mask_subgraph_2, :);
        dist_mat = pdist2(subgraph_nodes_1,subgraph_nodes_2, 'euclidean');
        subgraph_dist(i, j) = sum(dist_mat(:)) / (size(dist_mat, 1)*size(dist_mat, 2));
    end
end

if num_subgraph_1 <= num_subgraph_2
    [~, subgraph_corresp] = min(subgraph_dist, [], 2);
elseif num_subgraph_1 > num_subgraph_2
    [~, subgraph_corresp] = min(subgraph_dist, [], 1);
end

end
