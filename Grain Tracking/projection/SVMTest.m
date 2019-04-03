num_clusters = max(features_svm(:,end));
avg_normal_1 = zeros(num_clusters, 3);
avg_normal_2 = zeros(num_clusters, 3);
norm_dists = zeros(num_clusters, 3);

for k = 1:num_clusters
    mask_cluster_i_1 = (features_svm(:,end) == k) & (features_svm(:,1) == 1);
    mask_cluster_i_2 = (features_svm(:,end) == k) & (features_svm(:,1) == 2);
    avg_normal_1(k, :) = sum(features_svm(mask_cluster_i_1, 6:8))/sum(mask_cluster_i_1);
    avg_normal_1(k, :) = avg_normal_1(k, :)/norm(avg_normal_1(k, :));
    avg_normal_2(k, :) = sum(features_svm(mask_cluster_i_2, 6:8))/sum(mask_cluster_i_2);
    avg_normal_2(k, :) = avg_normal_2(k, :)/norm(avg_normal_2(k, :));
    
    norm_dists(k, 1) = abs(atand(norm(cross(avg_normal_1(k, :),avg_normal_2(k, :))) / dot(avg_normal_1(k, :),avg_normal_2(k, :))));
    norm_dists(k, 2) = abs(atand(norm(cross(avg_normal_1(k, :),svm_plane_info(k, 2:4))) / dot(avg_normal_1(k, :),svm_plane_info(k, 2:4))));
    norm_dists(k, 3) = abs(atand(norm(cross(svm_plane_info(k, 2:4),avg_normal_2(k, :))) / dot(svm_plane_info(k, 2:4),avg_normal_2(k, :))));
end

norm_dists



plotMedianPlanes(features_svm_naive, facetri_nodeid_an4, facetri_nodeid_an5, x_to_y, svm_plane_info_naive)
figure
plotMedianPlanes(features_svm, facetri_nodeid_an4, facetri_nodeid_an5, x_to_y, svm_plane_info)
figure
plotMedianPlanes(features_lr, facetri_nodeid_an4, facetri_nodeid_an5, x_to_y, lr_plane_info)



%%
OptionZ.FrameRate=20; OptionZ.Duration=8; OptionZ.Periodic=true;
file_name = ['SVM_pair_', num2str(i)];
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], file_name ,OptionZ)




%%
OptionZ.FrameRate=20; OptionZ.Duration=8; OptionZ.Periodic=true;
file_name = ['LR_pair_', num2str(i)];
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], file_name ,OptionZ)



%%
OptionZ.FrameRate=20; OptionZ.Duration=8; OptionZ.Periodic=true;
file_name = ['SVMimprove_pair_', num2str(i)];
CaptureFigVid([-20,10;-110,10;-190,80;-290,10;-380,10], file_name ,OptionZ)
