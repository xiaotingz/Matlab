% function face_migration = calcFaceMigBySVMPlaneProj(coord_1, coord_2, normal_1, normal_2, x_to_y, tol)
% ############################################################################
% Input
%     - coord_1 = [m, 3], coordinates of the unique mesh nodes in an4
%     - coord_2 = [n, 3], coordinates of the unique mesh nodes in an5
%     - normal_ = [ , 3], normals of the nodes, calced as the average normal
%     of triangles sharing that node
%     - x_to_y = [m, 1], a mapping from nodes_1 to nodes_2
% Output
%     - face_migration, scalar, the average migration distance of the grain face
% Notes
%     - This script can be compared to the other two projection methods
%     also: MigrationByLocalNormProj.m and MigrationByPillarHeight.m
%     - The logic is as follows: try to find some good median planes
%     between the two   ore than one median plane would be needed if the 
%     planes are intersecting, interwining, or curved planes. Details see 
%     Migration Projection notes.
%     - If too many (>10) median planes are needed, the program will break.
%     - This script is used in calcMigrations.m, the result can be compared to 
%     calcFaceMigByLocalNormProj.m and calcFaceMigByPillarHeight.m
% ############################################################################
% ----------------------- load debug data -----------------------
load('181004_SVMcurv_pair1674')

tol = 0.05;
% ---------------------------------------------------------------
colors = get(gca,'colororder');

% ##### Prepare to Start #####
max_kmeans_loop = 10000;
coord_1_corresp = coord_2(x_to_y, :); 
dist_direct = vecnorm(coord_1 - coord_1_corresp, 2, 2);

% ##### Initial SVM trial #####
% ----- Find SVM plane -----
X = [coord_1; coord_2];
Y = [ones(length(coord_1), 1); 2*ones(length(coord_2), 1)];
SVMModel = fitcsvm(X, Y, 'KernelFunction', 'linear');
normal = SVMModel.Beta;
bias = SVMModel.Bias;

% ##### Project to see if this plane is good enough #####
% """
% - Good median plane is characterized as most projected migration distance 
% should be smaller than the direct migration distance.  
% - If the plane is not good enough, add in a cluster fit median plane within 
% the clusters.
% - Only project_dist_1 is used to check if need_another_cluster.
% - The tol parameter is discussable. 
% - Point-plane distance ref: http://mathworld.wolfram.com/Point-PlaneDistance.html
% """
dist_proj_1 = (normal(1)*coord_1(:,1) + normal(2)*coord_1(:,2) + normal(3)*coord_1(:,3) ...
    + bias*ones(size(coord_1, 1), 1)) ./ norm(normal);
dist_proj_2 = (normal(1)*coord_2(:,1) + normal(2)*coord_2(:,2) + normal(3)*coord_2(:,3) ...
    + bias*ones(size(coord_2, 1), 1)) ./ norm(normal);
need_another_cluster = sum(dist_direct - abs(dist_proj_1) < 0) / length(dist_direct) > tol; 

% ##### Kmeans with Dynamic #Centers #####
num_cluster = 1;
m = size(coord_1, 1);
n = size(coord_2, 1);
% --- features = [face_id, coordinates, normals, cluster_id] ---
features = [[ones(m, 1); 2*ones(n, 1)], [node_id_an4; node_id_an5], [coord_1; coord_2], [normal_1; normal_2], ones(m+n, 1)];

% while need_another_cluster 
num_cluster = num_cluster + 1;
if num_cluster > 10
    disp(' '); warning('Too many SVM median planes are needed in MigrationBySVMPlaneProj.m. Check the node!'); disp(' ');
    return
end






num_cluster = 

% ----- Initialize cluster center with Kmeans++ -----
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
    end
    P = dist.*dist/sum(dist.*dist);
    cluster_center(i, :) = features(randsample(m+n,1,true,P), :);
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

% ----- Fit SVM plane to every cluster -----





% ----- Project to see if num_clusters is enough -----






figure
visualizeFace(face_node_info, x_to_y)
for i = 1:length(num_cluster)
    scatter3(features(cluster_id_new==i, 3), features(cluster_id_new==i, 4), features(cluster_id_new==i, 5), 80, 'filled')
end
title('Weighted K-means')

% scale = ceil(max(length(coord_1), length(coord_2))/100);
% x_range = [min(X(:,1)), max(X(:,1))];
% y_range = [min(X(:,2)), max(X(:,2))];
% plotSVMPlane(normal, bias, x_range, y_range, scale);

% ################################## Reference: Matlab K-means ################################## 
idx = kmeans([coord_1;coord_2],2, 'MaxIter', 100);
figure
visualizeFace(face_node_info, x_to_y)
scatter3(features(idx==1, 3), features(idx==1, 4), features(idx==1, 5), 80, 'filled')
scatter3(features(idx==2, 3), features(idx==2, 4), features(idx==2, 5), 80, 'filled')
title('Naive K-means')






