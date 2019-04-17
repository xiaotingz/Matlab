function [subgraph_1, subgraph_2] = findSubgraph(face_unique_nodeid_1, face_unique_nodeid_2, face_tri_nodeid_1, face_tri_nodeid_2)
% ############################################################################
% * Input 
%     - face_unique_nodeid_ = [n, 1], generation see load data
%     - face_tri_nodeid_ = [m, 3], generation see load data
% * Output
%     - subgraph_ = [n, 1], containing the subgraph_id for every node in face_unique_nodeid_
% * Notes
%     - use in main_TrackNodes.m
% ############################################################################
% ------------------ load data for debug --------------------
% load('180822_FaceCorresp');
% load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822.mat');
% idx = 1033;
% obj_facelabel_1 = tracked_uniqueface_an4(idx, :);
% obj_facelabel_2 = tracked_uniqueface_an5(idx, :);
% 
% mask_objface_1 = (facelabel_an4(:,1) == obj_facelabel_1(1) & facelabel_an4(:,2) == obj_facelabel_1(2));
% mask_objface_2 = (facelabel_an5(:,1) == obj_facelabel_2(1) & facelabel_an5(:,2) == obj_facelabel_2(2));
% 
% face_tri_nodeid_1 = tri_node_an4(mask_objface_1, :);
% face_unique_nodeid_1 = unique(face_tri_nodeid_1);
% face_tri_nodeid_2 = tri_node_an5(mask_objface_2, :);
% face_unique_nodeid_2 = unique(face_tri_nodeid_2);
% % -----------------------------------------------------------

% ##### Construct Graph of the Face #####
% ----- convert nodes to be local, otherwise graph is too big -----
face_map_1 = containers.Map(face_unique_nodeid_1, 1:length(face_unique_nodeid_1));
face_map_2 = containers.Map(face_unique_nodeid_2, 1:length(face_unique_nodeid_2));
tri_localnode_1 = zeros(size(face_tri_nodeid_1));
tri_localnode_2 = zeros(size(face_tri_nodeid_2));
for i = 1:size(face_tri_nodeid_1, 1)
    for j = 1:3
        tri_localnode_1(i, j) = face_map_1(face_tri_nodeid_1(i, j));
    end
end
for i = 1:size(face_tri_nodeid_2, 1)
    for j = 1:3
        tri_localnode_2(i, j) = face_map_2(face_tri_nodeid_2(i, j));
    end
end

% ----- get the unique edges (unique connections) -----
edge_1 = [tri_localnode_1(:,1), tri_localnode_1(:, 2); tri_localnode_1(:,2), tri_localnode_1(:, 3)];
edge_1 = sort(edge_1, 2);
edge_1 = unique(edge_1, 'rows');
edge_2 = [tri_localnode_2(:,1), tri_localnode_2(:, 2); tri_localnode_2(:,2), tri_localnode_2(:, 3)];
edge_2 = sort(edge_2, 2);
edge_2 = unique(edge_2, 'rows');
% ----- make adjacency matrix -----
adjmat_1 = sparse(zeros(length(face_unique_nodeid_1)));
adjmat_ind_1 = sub2ind(size(adjmat_1), edge_1(:,1), edge_1(:,2));
adjmat_1(adjmat_ind_1) = 1;
adjmat_2 = zeros(length(face_unique_nodeid_2));
adjmat_ind_2 = sub2ind(size(adjmat_2), edge_2(:,1), edge_2(:,2));
adjmat_2(adjmat_ind_2) = 1;
% ----- make graph and identify subgraph -----
graph_1 = graph(adjmat_1, 'upper');
graph_2 = graph(adjmat_2, 'upper');

% ##### Find Disconnected Subgraphs #####
subgraph_1 = conncomp(graph_1)';
subgraph_2 = conncomp(graph_2)';

end

% % ##### Check: Visualize Subgraphs #####
% % ----- need node_coord_an4 and node_coord_an4 -----
% color1 = [0, 0.4470, 0.7410];
% color2 = [0.9290, 0.6940, 0.1250];
% 
% face_node_coord_1 = node_coord_an4(face_unique_nodeid_1,:);
% face_node_coord_2 = node_coord_an5(face_unique_nodeid_2,:);
% 
% trisurf(face_tri_nodeid_1, node_coord_an4(:,1), node_coord_an4(:,2), node_coord_an4(:,3),'Facecolor',color1, 'Facealpha', 0.3, 'edgealpha', 0.3);
% hold on
% rotate3d on
% for i = 1:length(unique(subgraph_1))
%     scatter3(face_node_coord_1((subgraph_1==i), 1), face_node_coord_1((subgraph_1==i), 2), face_node_coord_1((subgraph_1==i), 3), 30, 'filled');
% end
% hold off
% % print(['pair_', num2str(idx), '_subgraph_1'], '-dpng','-r300')
% 
% figure()
% trisurf(face_tri_nodeid_2, node_coord_an5(:,1), node_coord_an5(:,2), node_coord_an5(:,3),'Facecolor',color2, 'Facealpha', 0.3, 'edgealpha', 0.3);
% hold on
% rotate3d on
% for i = 1:length(unique(subgraph_2))
%     scatter3(face_node_coord_2((subgraph_2==i), 1), face_node_coord_2((subgraph_2==i), 2), face_node_coord_2((subgraph_2==i), 3), 30, 'filled');
% end
% % print(['pair_', num2str(idx), '_subgraph_2'], '-dpng','-r300')