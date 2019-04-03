% function [mig_proj_abs, mig_proj_sign, features, svm_plane_info]  = calcFaceMigBySVMPlaneProj(features, eps, x_to_y, tri_node_1, tri_node_2)
% 
% ############################################################################
% Input
%     - features = [n+n', 10], [face_id, node_id, coordinates, normals, curves, cluster_id] 
%         including nodes of both faces
%     - x_to_y = [n, 1], a mapping from nodes_1 to nodes_2
%         this is to calculate direct node_pair_distance, for determining good median planes.
%     - eps, scalar, the tolarence for deciting if a fitted SVM median plane is good.
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
%     - This script can be compared to the other two projection methods
%     also: MigrationByLinearRegPlaneProj.m, MigrationByLocalNormProj.m, MigrationByPillarHeight.m
%     - Dependencies: updateKmeansClusterAssign.m and KmeansAsNormalImproved.m
%     - The logic is as follows: try to find some good median planes
%     between the two faces and project distance between tracked node pairs
%     onto the mediam plane. 
%     - More than one median plane would be needed if the planes are intersecting, 
%     interwining, or curved planes. In such cases, nodes need to be
%     clustered and median planes are fitted within each cluster. 
%     - K-means clustering, break if too many (>10) median planes are needed.
%     - MATLAB svm paramters to play with 
%         !!! 'BoxConstraint', 'Standardize', 'OutlierFraction' !!!         
% ############################################################################
% ----------------------- load debug data -----------------------
% load('181004_SVMcurv_pair1674')
eps = 0.2;
naive_kmeans = true;
% ---------------------------------------------------------------

% if nargin == 3
%     naive_kmeans = true;
% elseif nargin == 5
%     naive_kmeans = false;
% end

% ##### Prepare to Start #####
max_kmeans_loop = 1000;
coord_1 = features(features(:,1)==1, 3:5);
normal_1 = features(features(:,1)==1, 6:8);
coord_2 = features(features(:,1)==2, 3:5);
normal_2 = features(features(:,1)==2, 6:8);
coord_1_corresp = coord_2(x_to_y, :); 
% normal_1_corresp = normal_2(x_to_y, :);
dist_direct = vecnorm(coord_1 - coord_1_corresp, 2, 2);
max_num_cluster = 10;

% ################################## Fit First SVM Plane ##################################
% """
% - Good median plane is characterized as most projected migration distance 
% should be smaller than the direct migration distance. 
% - If one SVM median plane is not good, needs to cluster nodes and fit
% median planes within each cluster.
% """
% ----- Find SVM plane -----
X = features(:,3:5);
Y = features(:,1);
SVMModel = fitcsvm(X, Y, 'KernelFunction', 'linear');
normal = SVMModel.Beta;
bias = SVMModel.Bias;

% ---- Project to see if this plane is good enough -----
% """ 
% - A median plane is bad if: 1) there exists many nodes for which the node
% to median plane distance is larger than the direct distance between node
% pairs. 2) normal of the median plane is far away from the average normal
% of either planes. 
% - Point-plane distance ref: http://mathworld.wolfram.com/Point-PlaneDistance.html
% - Vector projection ref: https://math.oregonstate.edu/home/programs/undergrad/CalculusQuestStudyGuides/vcalc/dotprod/dotprod.html
% """
dist_proj_1 = (sum(coord_1.*normal', 2) + bias*ones(size(coord_1, 1), 1)) ./ norm(normal);
dist_proj_2 = (sum(coord_1_corresp.*normal', 2) + bias*ones(size(coord_1_corresp, 1), 1)) ./ norm(normal);
dist_proj_sign = dist_proj_1 - dist_proj_2;
dist_proj_abs = abs(dist_proj_1) + abs(dist_proj_2);
proj_dists_toolong = (sum(dist_direct - dist_proj_abs < 0) / length(dist_direct)) > eps; 
avg_normal_1 = sum(normal_1) / size(normal_1, 1);
avg_normal_2 = sum(normal_2) / size(normal_2, 1);
avg_normal_dists = abs(atand(norm(cross(avg_normal_1,avg_normal_2)) / dot(avg_normal_1,avg_normal_2)));
median_normal_dist_1 = abs(atand(norm(cross(avg_normal_1,normal)) / dot(avg_normal_1,normal)));
median_normal_dist_2 = abs(atand(norm(cross(avg_normal_2,normal)) / dot(avg_normal_2,normal)));
plane_normal_bad = (median_normal_dist_1 > avg_normal_dists*(1+eps) && median_normal_dist_2 > avg_normal_dists*(1+eps));
need_another_cluster = proj_dists_toolong || plane_normal_bad;
bad_fit = false;

% ----- Record the SVM plane parameters -----
if ~ need_another_cluster
    svm_plane_info = [1, normal', bias];
end


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
    svm_plane_info = [];
    mask_to_cluster = boolean(ones(size(features, 1), 1));
    dist_proj_sign = zeros(size(dist_direct));
end
while need_another_cluster
    num_cluster = num_cluster + 1;
%     disp(['num_good_clusters = ', num2str(num_good_clusters), ',  num_clusters = ', num2str(num_cluster)])
    
    % ----- Check if reached max #cluster -----
    if num_cluster > max_num_cluster
        disp(' '); warning('Too many SVM median planes are needed in MigrationBySVMPlaneProj.m. Check the node!'); disp(' ');
        mask_goodcluster = (features(:,end) <= num_good_clusters );
        if sum(mask_goodcluster) < size(features, 1) * (1 - eps)
            bad_fit = true;
        end
        break
    end
    
    % ----- Recluster nodes in the bad clusters only -----
    if naive_kmeans
        cluster_id_new = kmeans(features(mask_to_cluster, 3:5), num_cluster - num_good_clusters);
    else
        [~, cluster_id_new] = KmeansAsNormalImproved(features, mask_to_cluster, num_cluster - num_good_clusters, tri_node_1, tri_node_2);
    end
    
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
            X = features(mask_cluster_i, 3:5);
            Y = features(mask_cluster_i, 1);
            SVMModel = fitcsvm(X, Y,'KernelFunction','linear', 'Standardize',false);
            normal = SVMModel.Beta;
            bias = SVMModel.Bias;

             % ----- Check if cluster_i is good -----
            mask_cluster_i_1 = mask_cluster_i(1 : size(coord_1, 1));
            mask_cluster_i_2 = mask_cluster_i(size(coord_1, 1)+1 : end);
            dist_proj_1 = (sum(coord_1(mask_cluster_i_1, :).*normal', 2) ...
                + bias*ones(sum(mask_cluster_i_1), 1)) ./ norm(normal);
            dist_proj_2 = (sum(coord_1_corresp(mask_cluster_i_1, :).*normal', 2) ...
                + bias*ones(sum(mask_cluster_i_1), 1)) ./ norm(normal);
            dist_proj_abs(mask_cluster_i_1) = abs(dist_proj_1) + abs(dist_proj_2);
            proj_dists_toolong =  (sum(dist_direct(mask_cluster_i_1) - dist_proj_abs(mask_cluster_i_1) < 0) ...
                / sum(mask_cluster_i_1)) > eps;
            avg_normal_1 = sum(normal_1(mask_cluster_i_1,:), 1) / sum(mask_cluster_i_1);
%             avg_normal_2 = sum(normal_1_corresp(mask_cluster_i_1,:), 1) / sum(mask_cluster_i_1);
            avg_normal_2 = sum(normal_2(mask_cluster_i_2,:), 1) / sum(mask_cluster_i_2);
            avg_normal_dists = abs(atand(norm(cross(avg_normal_1,avg_normal_2)) / dot(avg_normal_1,avg_normal_2)));
            median_normal_dist_1 = abs(atand(norm(cross(avg_normal_1,normal)) / dot(avg_normal_1,normal)));
            median_normal_dist_2 = abs(atand(norm(cross(avg_normal_2,normal)) / dot(avg_normal_2,normal)));
            plane_normal_bad = (median_normal_dist_1 > avg_normal_dists*(1+eps) && median_normal_dist_2 > avg_normal_dists*(1+eps));
            only_one_face = (sum(mask_cluster_i_1)==0 || sum(mask_cluster_i_2)==0);
            this_cluster_bad = proj_dists_toolong || plane_normal_bad || only_one_face;


            if this_cluster_bad == true
                mask_to_cluster(mask_cluster_i) = true;
                need_another_cluster = true;
            else
                num_good_clusters = num_good_clusters + 1;
                svm_plane_info = [svm_plane_info; num_good_clusters, normal', bias];
                dist_proj_sign(mask_cluster_i_1) = dist_proj_1 - dist_proj_2;
                features(mask_cluster_i, end) = num_good_clusters * ones(sum(mask_cluster_i), 1);
            end
        end
    end 
end

if bad_fit 
    mig_proj_abs = NaN;
    mig_proj_sign = NaN;
    mask_bad_clusters = (features(:, end) > num_good_clusters);
    features(mask_bad_clusters, end) = -1;
elseif num_good_clusters == 0
    mig_proj_abs = sum(dist_proj_abs)/length(dist_proj_abs);
    mig_proj_sign = abs(sum(dist_proj_sign))/length(dist_proj_sign);
else
    cluster_ids_face1 = features(features(:,1)==1, end);
    mask_good_dists = (cluster_ids_face1 <= num_good_clusters);
    mig_proj_abs = sum(dist_proj_abs(mask_good_dists))/sum(mask_good_dists);
    mig_proj_sign = abs(sum(dist_proj_sign(mask_good_dists)))/sum(mask_good_dists);
    mask_bad_clusters = (features(:, end) > num_good_clusters);
    features(mask_bad_clusters, end) = -1;
end

% end

% plotSVMPlane(features, face_tri_node_an4, face_tri_node_an5, x_to_y)
% title('Naive K-means')




% mahal_dist = mahal(face_tri_normal_an4, face_tri_normal_an4);
% prob = 1 - chi2cdf(mahal_dist,3);




