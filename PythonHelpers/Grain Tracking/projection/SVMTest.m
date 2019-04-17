features_test = features_lr_all;
median_plane_test = lr_plane_info_all;

num_clusters = max(features_test(:,end));
avg_normal_1 = zeros(num_clusters, 3);
avg_normal_2 = zeros(num_clusters, 3);
norm_dists = zeros(num_clusters, 3);

for k = 1:num_clusters
    mask_cluster_i_1 = (features_test(:,end) == k) & (features_test(:,1) == 1);
    mask_cluster_i_2 = (features_test(:,end) == k) & (features_test(:,1) == 2);
    avg_normal_1(k, :) = sum(features_test(mask_cluster_i_1, 6:8), 1);
    avg_normal_2(k, :) = sum(features_test(mask_cluster_i_2, 6:8), 1);
    
    norm_dists(k, 1) = abs(atand(norm(cross(avg_normal_1(k, :),avg_normal_2(k, :))) / dot(avg_normal_1(k, :),avg_normal_2(k, :))));
    norm_dists(k, 2) = abs(atand(norm(cross(avg_normal_1(k, :),median_plane_test(k, 2:4))) / dot(avg_normal_1(k, :),median_plane_test(k, 2:4))));
    norm_dists(k, 3) = abs(atand(norm(cross(median_plane_test(k, 2:4),avg_normal_2(k, :))) / dot(median_plane_test(k, 2:4),avg_normal_2(k, :))));
end

norm_dists

%%

plotMedianPlanes(features_svm_naive, facetri_nodeid_an4, facetri_nodeid_an5, x_to_y, svm_plane_info_naive)
% figure
% plotMedianPlanes(features_svm_improve, facetri_nodeid_an4, facetri_nodeid_an5, x_to_y, svm_plane_info_improve)
figure
plotMedianPlanes(features_lr, facetri_nodeid_an4, facetri_nodeid_an5, x_to_y, lr_plane_info)

%%
OptionZ.FrameRate=20; OptionZ.Duration=8; OptionZ.Periodic=true;
file_name = ['SVM_nearest_pair_', num2str(i)];
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], file_name ,OptionZ)

%%
OptionZ.FrameRate=20; OptionZ.Duration=8; OptionZ.Periodic=true;
file_name = ['LR_nearest_pair_', num2str(i)];
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], file_name ,OptionZ)

%%
figure
plotMedianPlanes(features_lr_tracked, facetri_nodeid_an4, facetri_nodeid_an5, x_to_y, lr_plane_info_tracked)

%%
OptionZ.FrameRate=20; OptionZ.Duration=8; OptionZ.Periodic=true;
file_name = ['pair_', num2str(i)];
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], file_name ,OptionZ)


%% ##################### Test Sign #####################
% 
% [mig_lr_abs, mig_lr_sign, features_lr, lr_plane_info] = calcFaceMigByLinearRegPlaneProj(features, 0.2, x_to_y);
% [mig_lr_abs_all, mig_lr_sign_all, features_lr_all, lr_plane_info_all] = calcFaceMigByMedianPlaneProj(features, 0.2, x_to_y, ...
%             facetri_nodeid_an4, facetri_nodeid_an5, 'naive', 'lr', 'all');


plotSingleFaceWithNormal(file_an4, tracked_uniqueface_an4(1500, :), 0)

colors = get(gca,'colororder');
plotMedianPlanes(features_lr_all, facetri_nodeid_an4, facetri_nodeid_an5, x_to_y, lr_plane_info_all)
% figure
node_curvs = features(features(:,1) == 1, 9);
node_curv_pos = node_curvs > 0;
pos_curv_coord = facenode_coord_an4(node_curv_pos, :);
neg_curv_coord = facenode_coord_an4(~node_curv_pos, :);

scatter3(pos_curv_coord(:,1), pos_curv_coord(:,2), pos_curv_coord(:,3), ...
            80, 'filled', 'MarkerFaceColor',colors(2,:), 'MarkerEdgeColor',colors(2,:))
        
scatter3(neg_curv_coord(:,1), neg_curv_coord(:,2), neg_curv_coord(:,3), ...
            80, 'filled', 'MarkerFaceColor',colors(1,:), 'MarkerEdgeColor',colors(1,:))
        
rotate3d on


% plotMedianPlanes(features_lr, facetri_nodeid_an4, facetri_nodeid_an5, x_to_y, lr_plane_info)



%%
OptionZ.FrameRate=20; OptionZ.Duration=8; OptionZ.Periodic=true;
file_name = ['SVMimprove_pair_', num2str(i)];
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], file_name ,OptionZ)
