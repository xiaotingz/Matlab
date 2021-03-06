function [mig_proj_abs, mig_proj_sign, features, median_plane_info]  = calcFaceMigByMedianPlaneProj(features, eps, x_to_y, ...
              tri_node_1, tri_node_2, kmeans_method, classify_algorithm)
% % 
% ############################################################################
% Input
%     - features = [n+n', 10], [face_id, node_id, coordinates, normals, curves, cluster_id] 
%         including nodes of both faces
%     - eps, scalar, the tolarence for deciting if a fitted SVM median plane is good.
%     - x_to_y = [n, 1], a mapping from nodes_1 to nodes_2
%         this is to calculate direct node_pair_distance, for determining good median planes.
%     - kmeans
%         'naive' | 'improve'
%             'improve' means normal improved kmeans, may help SVM but probably won't affect LR
%     - classify_algorithm
%         'svm' | 'lr' (Linear Regression)
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
% % ----------------------- load debug data -----------------------
% eps = 0.2;
% kmeans_method = 'naive';
% classify_algorithm = 'lr';
% effective_nodes = 'all';
% tri_node_1 = facetri_nodeid_an4;
% tri_node_2 = facetri_nodeid_an5;
% % ---------------------------------------------------------------

% ##### Prepare to Start #####
% """
% The point using tracked nodes only is mainly to deal with the unrealisitc 
% large area change of twins. In such case, use track_by_nearest_node & effective_nodes = tracked. 
% """
max_num_cluster = 10;

if min(x_to_y) > 0
    all_nodes_effective = true;
else
    all_nodes_effective = false;
    mask_tracked_1 = (x_to_y > 0);
    x_to_y = (x_to_y(mask_tracked_1));
    
    % ----- if face too small -----
    if sum(mask_tracked_1) < 4
        mig_proj_abs = NaN;
        mig_proj_sign = NaN;
        median_plane_info = [];
        return
    end
        
    features_init = features;
    features_1 = features(features(:,1)==1, :);
    features_1 = features_1(mask_tracked_1, :);
    features_2 = features(features(:,1)==2, :);
    features_2 = features_2(x_to_y, :);
    features = [features_1; features_2];
end

coord_1 = features(features(:,1)==1, 3:5);
normal_1 = features(features(:,1)==1, 6:8);
coord_2 = features(features(:,1)==2, 3:5);
normal_2 = features(features(:,1)==2, 6:8);
if all_nodes_effective
    coord_1_corresp = coord_2(x_to_y, :);
else
    coord_1_corresp = coord_2;
end
mig_vec = coord_1_corresp - coord_1;
dist_direct = vecnorm(mig_vec, 2, 2);

curv_1 = features(features(:,1)==1, 9);
is_convex = sign(curv_1);
angdiff_mig_normal = abs(atan2d(norm(cross(mig_vec,normal_1, 2)), dot(mig_vec,normal_1, 2)));
mig_along_normal = sign(angdiff_mig_normal - 90);
sign_proj = is_convex .* mig_along_normal;
  
    
% ################################## Fit First Median Plane ##################################
% """
% - Good median plane is characterized as most projected migration distance 
% should be smaller than the direct migration distance. 
% - If a sinlge median plane is not good, needs to cluster nodes and fit
% median planes within each cluster.
% """
% ----- Find median plane -----
if strcmp(classify_algorithm, 'svm')
    X = features(:,3:5);
    Y = features(:,1);
    SVMModel = fitcsvm(X, Y, 'KernelFunction', 'linear','Standardize',false);
    normal = SVMModel.Beta;
    bias = SVMModel.Bias;
elseif strcmp(classify_algorithm, 'lr')
    X = features(:, 3:4);
    Y = features(:, 5);
%     mdl = fitlm(X, Y);
%     coeffs = mdl.Coefficients.Estimate;
    coeffs = robustfit(X, Y);
    norm_normal = norm([coeffs(2:3); -1]);
    normal = [coeffs(2:3); -1]./norm_normal;
    bias = coeffs(1)/norm_normal;
else 
    warning('classify_algorithm input of calcFaceMigByMedianPlaneProj.m is wrong! only surpport LR & SVM')
    return
end

% ---- Project to see if this plane is good enough -----
% """ 
% - A median plane is bad if: 1) there exists many nodes for which the node
% to median plane distance is larger than the direct distance between node
% pairs. 2) normal of the median plane is far away from the average normal
% of either planes. 
% - Note abs(dist_to_medianplane) can't be used as proj_dist. Still need to
% project dist_direct along normal of median plane. 
% - Point-plane distance 
%       ref: http://mathworld.wolfram.com/Point-PlaneDistance.html
% - Vector projection 
%       ref: https://math.oregonstate.edu/home/programs/undergrad/CalculusQuestStudyGuides/vcalc/dotprod/dotprod.html
% """
dist_proj_1 = (sum(coord_1.*normal', 2) + bias*ones(size(coord_1, 1), 1)) ./ norm(normal);
dist_proj_2 = (sum(coord_1_corresp.*normal', 2) + bias*ones(size(coord_1_corresp, 1), 1)) ./ norm(normal);
dist_to_medianplane = abs(dist_proj_1) + abs(dist_proj_2);
proj_dists_toolong = (sum(dist_direct - dist_to_medianplane < 0) / length(dist_direct)) > eps; 
avg_normal_1 = sum(normal_1) / size(normal_1, 1);
avg_normal_2 = sum(normal_2) / size(normal_2, 1);
avg_normal_dists = abs(atand(norm(cross(avg_normal_1,avg_normal_2)) / dot(avg_normal_1,avg_normal_2)));
median_normal_dist_1 = abs(atand(norm(cross(avg_normal_1,normal)) / dot(avg_normal_1,normal)));
median_normal_dist_2 = abs(atand(norm(cross(avg_normal_2,normal)) / dot(avg_normal_2,normal)));
plane_normal_bad = (median_normal_dist_1 > avg_normal_dists*(1+eps) && median_normal_dist_2 > avg_normal_dists*(1+eps));
need_another_cluster = proj_dists_toolong || plane_normal_bad;
bad_fit = false;

% ----- Record the median plane parameters -----
if ~ need_another_cluster
    median_plane_info = [1, normal', bias];
    num_nodes_1 = size(coord_1, 1);
    dist_proj_abs = abs(dot(mig_vec, repmat(normal', num_nodes_1, 1), 2)/norm(normal));
    dist_proj_sign = sign_proj .* dist_proj_abs;
else
    dist_proj_abs = zeros(size(coord_1, 1), 1);
    dist_proj_sign = zeros(size(coord_1, 1), 1);
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
    median_plane_info = [];
    mask_to_cluster = boolean(ones(size(features, 1), 1));
    dist_proj_sign = zeros(size(dist_direct));
end
while need_another_cluster
    num_cluster = num_cluster + 1;
%     disp(['num_good_clusters = ', num2str(num_good_clusters), ',  num_clusters = ', num2str(num_cluster)])
    
    % ----- Check if reached max #cluster -----
    if num_cluster > max_num_cluster
        disp(' '); warning('Too many median planes are needed in MigrationByMedianPlaneProj.m. Check the node!'); disp(' ');
        mask_goodcluster = (features(:,end) <= num_good_clusters );
        if sum(mask_goodcluster) < size(features, 1) * (1 - eps)
            bad_fit = true;
        end
        break
    end
    
    % ----- Recluster nodes in the bad clusters only -----
    if strcmp(kmeans_method, 'naive')
        cluster_id_new = kmeans(features(mask_to_cluster, 3:5), num_cluster - num_good_clusters);
    elseif strcmp(kmeans_method, 'improve')
        [~, cluster_id_new] = KmeansAsNormalImproved(features, mask_to_cluster, num_cluster - num_good_clusters, tri_node_1, tri_node_2);
    else
        warning('kmeans input of calcFaceMigByMedianPlaneProj.m is wrong! only surpport naive & improve')
        return
    end
    
    features(mask_to_cluster, end) = cluster_id_new + num_good_clusters;
    
    % ----- Fit meidan plane to every cluster -----
    mask_to_cluster = boolean(zeros(size(features, 1), 1));
    need_another_cluster = false;
    for i = 1+num_good_clusters : num_cluster
        mask_cluster_i = (features(:, end) == i);
        % --- if a cluster includes no more than 4 nodes, ignore this cluster ---
        if sum(mask_cluster_i) <= 4
            mask_to_cluster(mask_cluster_i) = true;
            continue
        % --- else, fit a median plane in this cluster ---
        else
            if strcmp(classify_algorithm, 'svm')
                X = features(mask_cluster_i, 3:5);
                Y = features(mask_cluster_i, 1);
                SVMModel = fitcsvm(X, Y, 'KernelFunction', 'linear', 'Standardize',false);
                normal = SVMModel.Beta;
                bias = SVMModel.Bias;
            else
                X = features(mask_cluster_i, 3:4);
                Y = features(mask_cluster_i, 5);
%                 mdl = fitlm(X, Y);
%                 coeffs = mdl.Coefficients.Estimate;
                coeffs = robustfit(X, Y);
                norm_normal = norm([coeffs(2:3); -1]);
                normal = [coeffs(2:3); -1]./norm_normal;
                bias = coeffs(1)/norm_normal;
            end

             % ----- Check if cluster_i is good -----
            mask_cluster_i_1 = mask_cluster_i(1 : size(coord_1, 1));
            mask_cluster_i_2 = mask_cluster_i(size(coord_1, 1)+1 : end);
            dist_proj_1 = (sum(coord_1(mask_cluster_i_1, :).*normal', 2) ...
                + bias*ones(sum(mask_cluster_i_1), 1)) ./ norm(normal);
            dist_proj_2 = (sum(coord_1_corresp(mask_cluster_i_1, :).*normal', 2) ...
                + bias*ones(sum(mask_cluster_i_1), 1)) ./ norm(normal);
            dist_to_medianplane(mask_cluster_i_1) = abs(dist_proj_1) + abs(dist_proj_2);
            proj_dists_toolong =  (sum(dist_direct(mask_cluster_i_1) - dist_to_medianplane(mask_cluster_i_1) < 0) ...
                / sum(mask_cluster_i_1)) > eps;
            avg_normal_1 = sum(normal_1(mask_cluster_i_1,:), 1);
            avg_normal_2 = sum(normal_2(mask_cluster_i_2,:), 1);
            avg_normal_dists = abs(atand(norm(cross(avg_normal_1,avg_normal_2)) / dot(avg_normal_1,avg_normal_2)));
            median_normal_dist_1 = abs(atand(norm(cross(avg_normal_1,normal)) / dot(avg_normal_1,normal)));
            median_normal_dist_2 = abs(atand(norm(cross(avg_normal_2,normal)) / dot(avg_normal_2,normal)));
            plane_normal_bad = (median_normal_dist_1 > avg_normal_dists*(1+eps) && median_normal_dist_2 > avg_normal_dists*(1+eps));
            only_one_face = (sum(mask_cluster_i_1)==0 || sum(mask_cluster_i_2)==0);
            if strcmp(classify_algorithm, 'svm')
                this_cluster_bad = proj_dists_toolong || plane_normal_bad || only_one_face;
            else
                this_cluster_bad = proj_dists_toolong || plane_normal_bad ;
            end

            if this_cluster_bad == true
                mask_to_cluster(mask_cluster_i) = true;
                need_another_cluster = true;
            else
                num_good_clusters = num_good_clusters + 1;
                median_plane_info = [median_plane_info; num_good_clusters, normal', bias];
                features(mask_cluster_i, end) = num_good_clusters * ones(sum(mask_cluster_i), 1);
                num_nodes_1 = sum(mask_cluster_i_1);
                dist_proj_abs(mask_cluster_i_1) = abs(dot(mig_vec(mask_cluster_i_1, :), repmat(normal', num_nodes_1, 1), 2) / norm(normal));
                dist_proj_sign(mask_cluster_i_1) = sign_proj(mask_cluster_i_1) .* dist_proj_abs(mask_cluster_i_1);
            end
        end
    end 
end

% --- calculate projected migration dist, use only nodes for which a median plane has been found ---
if bad_fit 
    mig_proj_abs = NaN;
    mig_proj_sign = NaN;
    mask_bad_clusters = (features(:, end) > num_good_clusters);
    features(mask_bad_clusters, end) = -1;
elseif num_good_clusters == 0
    mig_proj_abs = sum(dist_proj_abs)/length(dist_proj_abs);
    mig_proj_sign = sum(dist_proj_sign)/length(dist_proj_sign);
else
    cluster_ids_face1 = features(features(:,1)==1, end);
    mask_good_dists = (cluster_ids_face1 <= num_good_clusters);
    mig_proj_abs = sum(dist_proj_abs(mask_good_dists))/sum(mask_good_dists);
    mig_proj_sign = sum(dist_proj_sign(mask_good_dists))/sum(mask_good_dists);
    mask_bad_clusters = (features(:, end) > num_good_clusters);
    features(mask_bad_clusters, end) = -1;
end

% --- always return the initial full node feataures, but with the correct clusters (for debug purpose) ---
if ~ all_nodes_effective
    features_init_1 = features_init(features_init(:,1)==1, :);
    features_init_1(:,end) = -1;
    features_init_1(mask_tracked_1, end) = features(features(:,1)==1, end);
    features_init_2 = features_init(features_init(:,1)==2, :);
    features_init_2(:,end) = -1;
    features_init_2(x_to_y,end) = features(features(:,1)==2, end);
    features = [features_init_1; features_init_2];
end


end

% plotSVMPlane(features, face_tri_node_an4, face_tri_node_an5, x_to_y)
% title('Naive K-means')




% mahal_dist = mahal(face_tri_normal_an4, face_tri_normal_an4);
% prob = 1 - chi2cdf(mahal_dist,3);




