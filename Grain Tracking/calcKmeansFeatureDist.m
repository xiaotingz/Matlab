function dist = calcKmeansFeatureDist(cluster_center, feature_member)
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
%     - See also updateKmeansClusterAssign.m
%     - There are two elements that make the distance: coordinates and normals.
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
w_normal = sum(dist_pos)/length(dist_pos) / 0.3;
dist = dist_pos + w_normal*dist_normal;
dist = reshape(dist, m, n);

end