function dist_proj_1  = calcFaceMigBySVMPlaneProj(features, x_to_y, eps)
% ############################################################################
% Input
%     - features = [n+m, 9], [face_id, node_id, coordinates, normals, cluster_id] 
%         including nodes of both faces
%     - x_to_y = [m, 1], a mapping from nodes_1 to nodes_2
%     - eps, scalar, the tolarence for deciting if a SVM fitting plane is good.
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
%     calcFaceMigByLocalNormProj.m and calcFaceMigByPillarHeight.m.
%     - This script has a dependency: updateKmeansClusterAssign.m and
%     KmeansAsNormalImproved.m
% ############################################################################
% ----------------------- load debug data -----------------------
% load('181004_SVMcurv_pair1674')
% node_id_1 = face_node_id_an4;
% node_id_2 = face_node_id_an5;
% node_id_1 = node_id_an4;
% node_id_2 = node_id_an5;
% tri_node_1 = face_tri_node_an4;
% tri_node_2 = face_tri_node_an5;
% eps = 0.1;
% ---------------------------------------------------------------
colors = get(gca,'colororder');

% ##### Prepare to Start #####
max_kmeans_loop = 1000;
coord_1 = features(features(:,1)==1, 3:5);
coord_2 = features(features(:,1)==2, 3:5);
coord_1_corresp = coord_2(x_to_y, :); 
dist_direct = vecnorm(coord_1 - coord_1_corresp, 2, 2);

% ##### Initial SVM trial #####
% ----- Find SVM plane -----
X = features(:,3:5);
Y = features(:,1);
SVMModel = fitcsvm(X, Y, 'KernelFunction', 'linear');
normal = SVMModel.Beta;
bias = SVMModel.Bias;

% ##### Project to see if this plane is good enough #####
% """
% - Good median plane is characterized as most projected migration distance 
% should be smaller than the direct migration distance.  
% - If the plane is not good enough, add in a cluster fit median plane within 
% the clusters.
% - Point-plane distance ref: http://mathworld.wolfram.com/Point-PlaneDistance.html
% """
dist_proj_1 = (sum(coord_1.*normal', 2) + bias*ones(size(coord_1, 1), 1)) ./ norm(normal);
% dist_proj_2 = (normal(1)*coord_2(:,1) + normal(2)*coord_2(:,2) + normal(3)*coord_2(:,3) ...
%     + bias*ones(size(coord_2, 1), 1)) ./ norm(normal);
need_another_cluster = (sum(dist_direct - abs(dist_proj_1) < 0) / length(dist_direct)) > eps; 

% ##### Kmeans with Dynamic #Centers #####

% %% ################################## Normal Improved K-means ################################## 
% % """
% % Note in each iteration, all nodes are reclustered.
% % """
% while need_another_cluster
%     num_cluster = num_cluster + 1;
%     if num_cluster > 10
%         disp(' '); warning('Too many SVM median planes are needed in MigrationBySVMPlaneProj.m. Check the node!'); disp(' ');
%         return
%     end
% 
%     [features, ~] = KmeansAsNormalImproved(features, num_cluster, tri_node_1, tri_node_2);
% 
%     % ----- Fit SVM plane to every cluster -----
%     dist_proj_1 = zeros(size(dist_direct));
%     need_another_cluster = false;
%     for i = 1:max(features(:, end))
%         mask_cluster_i = (features(:, end) == i);
%         mask_cluster_i_1 = (features(:, end) == i & features(:,1) == 1);
%         X = features(mask_cluster_i, 3:5);
%         Y = features(mask_cluster_i, 1);
%         SVMModel = fitcsvm(X, Y,'KernelFunction','linear', 'Standardize',false,'ClassNames',{'1','2'});
%         normal = SVMModel.Beta;
%         bias = SVMModel.Bias;
%         dist_proj_1(mask_cluster_i_1) = (sum(features(mask_cluster_i_1, 3:5).*normal', 2) ...
%             + bias*ones(sum(mask_cluster_i_1), 1)) ./ norm(normal);
%          
%          % ----- Project to see if num_clusters is enough -----
%         this_cluster_bad =  (sum(dist_direct(mask_cluster_i_1) - abs(dist_proj_1(mask_cluster_i_1)) < 0) ...
%             / length(dist_direct(mask_cluster_i_1))) > tol;
%         need_another_cluster = (need_another_cluster || this_cluster_bad);
%     end
%     
% %      % ----- Project to see if num_clusters is enough -----
% %     need_another_cluster = (sum(dist_direct - abs(dist_proj_1) < 0) / length(dist_direct)) > tol;
% end
% 
% % % ----- Visualization -----
% visualizeFace(face_node_info, x_to_y)
% plotSVMPlane(features)
% title('Weighted K-means')
% %%

% ################################## Reference: Matlab Naive K-means ################################## 
% """
% Note in each iteration, only nodes in the imperfect clusters will be reclustered. Nodes that
% already found a good median plane will maintain their cluster_id. 
% """
num_cluster = 1;
num_good_clusters = 0;
dist_proj_1 = zeros(size(dist_direct));
mask_to_cluster = boolean(ones(size(features, 1), 1));
svm_param = [];
while need_another_cluster
    num_cluster = num_cluster + 1;
    if num_cluster > 10
        disp(' '); warning('Too many SVM median planes are needed in MigrationBySVMPlaneProj.m. Check the node!'); disp(' ');
        dist_proj_1 = NaN;
        return 
    end
    
    % ----- Increase #clusters and recluster nodes in the bad clusters only -----
    cluster_id_new = kmeans(features(mask_to_cluster, 3:5), num_cluster - num_good_clusters + 1);
%     [~, cluster_id_new] = KmeansAsNormalImproved(features, mask_to_cluster, num_cluster, tri_node_1, tri_node_2);
    
    features(mask_to_cluster, end) = cluster_id_new + num_good_clusters;
    
    % ----- Fit SVM plane to every cluster -----
    mask_to_cluster = boolean(zeros(size(features, 1), 1));
    need_another_cluster = false;
    for i = 1+num_good_clusters : max(features(:, end))
        mask_cluster_i = (features(:, end) == i);
        mask_cluster_i_1 = (features(:, end) == i & features(:,1) == 1);
        X = features(mask_cluster_i, 3:5);
        Y = features(mask_cluster_i, 1);
        SVMModel = fitcsvm(X, Y,'KernelFunction','linear', 'Standardize',false,'ClassNames',{'1','2'});
        normal = SVMModel.Beta;
        bias = SVMModel.Bias;
        dist_proj_tmp = (sum(features(mask_cluster_i_1, 3:5).*normal', 2) ...
            + bias*ones(sum(mask_cluster_i_1), 1)) ./ norm(normal);
         
         % ----- Project to see if num_clusters is enough -----
        this_cluster_bad =  (sum(dist_direct(mask_cluster_i_1) - abs(dist_proj_tmp) < 0) ...
            / length(dist_direct(mask_cluster_i_1))) > eps;
        
        if this_cluster_bad == true
            mask_to_cluster(mask_cluster_i) = true;
            need_another_cluster = true;
        else
            num_good_clusters = num_good_clusters + 1;
            svm_param = [svm_param; normal', bias];
            dist_proj_1(mask_cluster_i_1) = dist_proj_tmp;
            features(mask_cluster_i, end) = num_good_clusters * ones(sum(mask_cluster_i), 1);
        end
        
    end
    
end

end

% visualizeFace(face_node_info, x_to_y)
% plotSVMPlane(features)
% title('Naive K-means')








