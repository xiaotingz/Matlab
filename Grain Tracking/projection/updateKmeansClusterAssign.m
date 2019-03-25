function cluster_id_new = updateKmeansClusterAssign(cluster_center, feature_member, tri_node_1, tri_node_2)
% ############################################################################
% * Input 
%     - feature_obj = [n, 9]
%     - feature_member = [m, 9]
%         [face_id, node_id, coordinates, normals, cluster_id]
%     - tri_node_= [k, 3]
%         nodes connection on the triangles, written from SharedTriList of D3D
% * Output
%     - cluster_id_new = [m, 1]
% * Notes
%     - This function is used in the kmeans part of calcFaceMigBySVMPlane_Proj.m
%     - There are two elements that make the distance: coordinates and normals.
%     - In order to avoid disconnection, dilate the clusters after each
%     cluster_id update. 
% ############################################################################
% ------------------ load data for debug --------------------
% feature_member = features;
% feature_obj = feature_member(1:3, :);
% clearvars -except feature_obj feature_member
% -----------------------------------------------------------

% ----- Data preparation -----
n = size(cluster_center, 1);
m = size(feature_member, 1);
k = size(cluster_center, 2);

member = repmat(feature_member, n, 1);
obj = zeros(size(member));
pos = 1;
for i = 1:n
    obj(pos:pos+m-1, :) = ones(m, k)*diag(cluster_center(i,:));
    pos = pos+m;
end

% ----- Position dist, as euclidean distance -----
dist_pos = vecnorm(member(:, 3:5) - obj(:, 3:5), 2, 2);

% ----- Normal dist, as norm of cross product -----
dist_normal = vecnorm(cross(member(:, 6:8), obj(:, 6:8)), 2, 2);

% ----- Sign dist, as the number of different sign(proj_dist) -----
% """
% ref: https://stackoverflow.com/questions/9605556/how-to-project-a-point-onto-a-plane-in-3d
%  - Projection sign, is from the following fact: within the same cluster,
% nodes on the same face should try to lie on the save side of a median
% plane. The accurate way to get the median plane is fit SVM every time
% but that's too expensive. So we just use cluster center and the
% average plane normal within the cluster. 
% - However, signs are not easy to initialize. It's not an intrinsic property
% of the cluster members thus not an intrinsic property of the cluster
% center also. 
% """

% ----- Feature dist -----
% """
% - We are going to check for connection so large weights for direction is fine. 
% - It's not clear whether a scalar weight is better or a vector weight is
% better.
% - test, vecnorm(cross([1,0,0], [cosd(20), sind(20),0]),2,2) = 0.34.
% - Suppose we hope the variation of plane normal is within +-10deg 
% - Also, suppose normal and position should be equally important, then we
% set w_normal = avg_dist_pos/0.3
% """
% w_normal = sum(dist_pos)/length(dist_pos) / 0.3;
% dist = dist_pos + w_normal*dist_normal;
w_normal = dist_pos;
dist = dist_pos + w_normal.*dist_normal;

dist = reshape(dist, m, n);
[~, cluster_id_new] = min(dist, [], 2);

% ##### Check connections and erode disconnections #####
for i = 1:2
    % ----- Get candidate nodes from boundary nodes -----
    % """
    % If three nodes all belong to the same triangle, leave them alone.
    % Don't worry about missing nodes because there are duplications and
    % the boundary nodes will show themselves in other triangles. 
    % """
    if i == 1
        tri_node = tri_node_1;
    else
        tri_node = tri_node_2;
    end
    mask_face = (feature_member(:,1) == i);
    nodeid_local = (1:size(feature_member, 1))';
    nodeid_to_clusterid = containers.Map(feature_member(mask_face, 2), cluster_id_new(mask_face));
    nodeid_to_local = containers.Map(feature_member(mask_face, 2), nodeid_local(mask_face));
    tri_node_clusterid = zeros(size(tri_node));
    tri_nodeid_local = zeros(size(tri_node));
    for j = 1:size(tri_node, 1)
        for l = 1:3
            tri_node_clusterid(j, l) = nodeid_to_clusterid(tri_node(j, l));
            tri_nodeid_local(j, l) = nodeid_to_local(tri_node(j, l));
        end
    end
    mask_candidates = (tri_node_clusterid(:,1) + tri_node_clusterid(:,2) ~= 2*tri_node_clusterid(:,3));
    cand_localid = unique(tri_nodeid_local(mask_candidates, :));
    
    % ----- Order the candidates by their distance to cluster center ----- 
    cand_coord = feature_member(cand_localid, 3:5);
    cand_center_coord = cluster_center(cluster_id_new(cand_localid), 3:5);
    cand_center_dist = vecnorm(cand_coord - cand_center_coord, 2, 2);
    [~, idx_dist_ordered] = sort(cand_center_dist, 'descend');
    
    % ----- Do the erosion ----- 
    % """
    % Within the neighborhood of the node, if no more than 2 nodes share
    % the same clusterid with it, erode its cluster id to be its neighbors'
    % """
    for j = 1:length(idx_dist_ordered)
        node_id = cand_localid(idx_dist_ordered(j));
        node_clusterid = cluster_id_new(node_id);
        mask_node_connect = any(tri_nodeid_local == node_id, 2);
        node_neigh_id = unique(tri_nodeid_local(mask_node_connect, :));
        neigh_clusterid = cluster_id_new(node_neigh_id);
        if ( sum(neigh_clusterid == node_clusterid) < 3 )
            cluster_id_new(node_id) = mode(neigh_clusterid);
        end
    end
    
end

end