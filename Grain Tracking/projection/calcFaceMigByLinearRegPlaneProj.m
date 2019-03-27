function [mig_proj_abs, mig_proj_sign, features]  = calcFaceMigByLinearRegPlaneProj(features, eps, x_to_y)
% ############################################################################
% Input
%     - features = [n+n', 9], [face_id, node_id, coordinates, normals, cluster_id] 
%         including nodes of both faces
%     - x_to_y = [n, 1], a mapping from nodes_1 to nodes_2
%         this is to calculate direct node_pair_distance, for determining good median planes.
%     - eps, scalar, the tolarence for deciting if a fitted median plane is good.
% Output
%     - dist_proj_sign & dist_proj_abs
%         If interwining shouldn't be considered as migration, use
%         mig_sign. If interwining should be consider as migration, use
%         mig_abs
%     - bad_fit
%         If most nodes on the face has been fitted.
%     - features
%         the cluster_id has been updated.
% Notes
%     - This script can be compared to the other two projection methods also: 
%     MigrationBySVMPlaneProj.m, MigrationByLocalNormProj.m, MigrationByPillarHeight.m
%     - The idea is to treat the problem as data cleaning: regression for noise cleaning. 
%     - For details, see MigrationBySVMPlaneProj.m The logic is really the same.  
%         -- First, cluster face nodes so that a median plane can be found within a cluster.
%         -- Second, fit median plane by Linear Regression within each cluster.
% ############################################################################

% ##### Prepare to Start #####
coord_1 = features(features(:,1)==1, 3:5);
coord_2 = features(features(:,1)==2, 3:5);
coord_1_corresp = coord_2(x_to_y, :); 
dist_direct = vecnorm(coord_1 - coord_1_corresp, 2, 2);
max_num_cluster = 10;


% ##### Fit Linear Regression Model #####
% """
% Plane normal direction http://mathworld.wolfram.com/Plane.html
% """
X = features(:, 3:4);
y = features(:, 5);
mdl = fitlm(X, y);
coeffs = mdl.Coefficients.Estimate;
norm_normal = norm([coeffs(2:3); -1]);
normal = [coeffs(2:3); -1]./norm_normal;
bias = coeffs(1)/norm_normal;

% ---- Project to see if this plane is good enough -----
% """ 
% - A median plane is bad if: 1) there exists many nodes for which the node
% to median plane distance is larger than the direct distance between node
% pairs. 2) normal of the median plane is far away from the average normal
% of either planes. 
% - Point-plane distance ref: http://mathworld.wolfram.com/Point-PlaneDistance.html
% """
dist_proj_1 = (sum(coord_1.*normal', 2) + bias*ones(size(coord_1, 1), 1)) ./ norm(normal);
dist_proj_2 = (sum(coord_1_corresp.*normal', 2) + bias*ones(size(coord_1_corresp, 1), 1)) ./ norm(normal);
dist_proj_sign = dist_proj_1 - dist_proj_2;
dist_proj_abs = abs(dist_proj_1) + abs(dist_proj_2);
need_another_cluster = (sum(dist_direct - dist_proj_abs < 0) / length(dist_direct)) > eps; 
bad_fit = false;

% ################################## K-means ################################## 
% """
% - if nargin == 3 ---> Naive K-means
% - if nargin == 5 ---> Normal improved K-means
% - Note in each iteration, only nodes in the imperfect clusters will be reclustered. Nodes that
% already found a good median plane will maintain their cluster_id. 
% """
num_cluster = 1;
num_good_clusters = 0;
if need_another_cluster
    mask_to_cluster = boolean(ones(size(features, 1), 1));
    dist_proj_sign = zeros(size(dist_direct));
end
while need_another_cluster
    num_cluster = num_cluster + 1;
    
    % ----- Check if reached max #cluster -----
    if num_cluster > max_num_cluster
        disp(' '); warning('Too many clusters are needed in MigrationByLinearRegPlaneProj.m. Check the node!'); disp(' ');
        mask_goodcluster = (features(:,end) <= num_good_clusters );
        if sum(mask_goodcluster) < size(features, 1) * (1 - eps)
            bad_fit = true;
        end
        break
    end
    
    % ----- Recluster nodes in the bad clusters only -----
    cluster_id_new = kmeans(features(mask_to_cluster, 3:5), num_cluster - num_good_clusters + 1);
    
    features(mask_to_cluster, end) = cluster_id_new + num_good_clusters;
    
    % ----- Fit SVM plane to every cluster -----
    mask_to_cluster = boolean(zeros(size(features, 1), 1));
    need_another_cluster = false;
    for i = 1+num_good_clusters : num_cluster
        mask_cluster_i = (features(:, end) == i);
        % --- if a cluster includes no more than 4 nodes, ignore this cluster ---
        if sum(mask_cluster_i) <= 4
            mask_to_cluster(mask_cluster_i) = true;
            continue
        else
            X = features(mask_cluster_i, 3:4);
            y = features(mask_cluster_i, 5);
            mdl = fitlm(X, y);
            coeffs = mdl.Coefficients.Estimate;
            norm_normal = norm([coeffs(2:3); -1]);
            normal = [coeffs(2:3); -1]./norm_normal;
            bias = coeffs(1);

             % ----- Check if cluster_i is good -----
            mask_cluster_i_1 = mask_cluster_i(1:size(coord_1, 1));
            dist_proj_1 = (sum(coord_1(mask_cluster_i_1, :).*normal', 2) ...
                + bias*ones(sum(mask_cluster_i_1), 1)) ./ norm(normal);
            dist_proj_2 = (sum(coord_1_corresp(mask_cluster_i_1, :).*normal', 2) ...
                + bias*ones(sum(mask_cluster_i_1), 1)) ./ norm(normal);
            dist_proj_abs(mask_cluster_i_1) = abs(dist_proj_1) + abs(dist_proj_2);
            this_cluster_bad =  (sum(dist_direct(mask_cluster_i_1) - dist_proj_abs(mask_cluster_i_1) < 0) ...
                / sum(mask_cluster_i_1)) > eps;


            if this_cluster_bad == true
                mask_to_cluster(mask_cluster_i) = true;
                need_another_cluster = true;
            else
                num_good_clusters = num_good_clusters + 1;
                dist_proj_sign(mask_cluster_i_1) = dist_proj_1 - dist_proj_2;
                features(mask_cluster_i, end) = num_good_clusters * ones(sum(mask_cluster_i), 1);
            end
        end
    end 
end


if bad_fit 
    mig_proj_abs = NaN;
    mig_proj_sign = NaN;
elseif num_good_clusters == 0
    mig_proj_abs = sum(dist_proj_abs)/length(dist_proj_abs);
    mig_proj_sign = abs(sum(dist_proj_sign))/length(dist_proj_sign);
else
    cluster_ids_face1 = features(features(:,1)==1, end);
    mask_good_dists = (cluster_ids_face1 <= num_good_clusters);
    mig_proj_abs = sum(dist_proj_abs(mask_good_dists))/sum(mask_good_dists);
    mig_proj_sign = abs(sum(dist_proj_sign(mask_good_dists)))/sum(mask_good_dists);
end


end






