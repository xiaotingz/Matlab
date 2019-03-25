function [features, cluster_id_new] = KmeansAsNormalImproved(features, num_cluster, tri_node_1, tri_node_2)
% ############################################################################
% NOTES
%   - This function is to be used in calcFaceMigBySVMPlaneProj.m
%   - This script has a dependency: updateKmeansClusterAssign.m
% Input
%   - features = [9, n], [face_id, node_id, coordinates, normals, cluster_id]
%   - num_cluster = scalar
% ############################################################################
max_kmeans_loop = 1000;

m = sum(features(:,1) == 1);
n = sum(features(:,1) == 2);

% ----- Kmeans++, Initialize cluster centers -----
% """
% ref: https://www.mathworks.com/help/stats/kmeans.html#bueftl4-1
% """
cluster_center = zeros(num_cluster, size(features, 2));
cluster_center(1, :) = features(randi(size(features, 1)), :);
dist = calcKmeansFeatureDist(cluster_center(1, :), features);
P = dist.*dist/sum(dist.*dist);
cluster_center(2, :) = features(randsample(m+n,1,true,P), :);
if num_cluster > 2
    for i = 3 : num_cluster
        dist = calcKmeansFeatureDist(cluster_center(1, :), features);
        for j = 2 : i-1
            dist_tmp = calcKmeansFeatureDist(cluster_center(j, :), features);
            dist(dist > dist_tmp) = dist_tmp(dist > dist_tmp);
        end
        P = dist.*dist/sum(dist.*dist);
        cluster_center(i, :) = features(randsample(m+n,1,true,P), :);
    end
end

% ----- Cluster nodes by Kmeans -----
kmeans_not_converge = true;
num_loop = 1;
while kmeans_not_converge 
    % ----- Update cluster assignment -----
    cluster_id_old = features(:, end);
    cluster_id_new = updateKmeansClusterAssign(cluster_center, features, tri_node_1, tri_node_2);
    features(:, end) = cluster_id_new;

    % ----- Update cluster centers -----
    % """
    % notice the first two columns of cluster center are meaningless
    % """
    for i = 1:num_cluster
        mask_cluster = (features(:,end) == i);
        cluster_center(i, :) = sum(features(mask_cluster, :))/sum(mask_cluster);
    end

    % ----- Check if converged -----
    kmeans_not_converge = (sum(cluster_id_old ~= cluster_id_new) > 0);

    % ----- Check if max num_loop exceeded -----
    num_loop = num_loop + 1;
    if num_loop > max_kmeans_loop
        warning('max_kmeans_loop reached but still not converge');
        return 
    end
end

end