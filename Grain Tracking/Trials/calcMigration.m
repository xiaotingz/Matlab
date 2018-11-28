% % function calcMigration(file_an4, file_an4, X_to_Y)
% % ############################################################################
% % Notes
% %     - This function is doing data preparation.
% %     - Migration of each grain face is calced by calling
% %     calcFaceMigBySVMPlaneProj.m | calcFaceMigByPillarHeight.m |
% %     calcFaceMigByLocalNormProj.m
% % ############################################################################
% % ----------------------- load debug data -----------------------
% % load('180822', 'file_an4', 'file_an5');
% file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
% file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh.dream3d';
% load('180822_FaceCorresp', 'tracked_uniqueface_an4', 'tracked_uniqueface_an5', 'X_to_Y');
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/180828_piecewise_face.mat', 'face_piecewise')
% % ---------------------------------------------------------------
% 
% % ##### Load Data & Basic Cleans #####
% facelabel_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
% facelabel_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
% tri_normal_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
% tri_normal_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
% tri_node_an4 = 1 + double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
% tri_node_an5 = 1 + double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
% node_coord_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
% node_coord_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
% mask_an4 = all(facelabel_an4 > 0, 2);
% facelabel_an4 = facelabel_an4(mask_an4, :);
% tri_node_an4 = tri_node_an4(mask_an4, :);
% tri_normal_an4 = tri_normal_an4(mask_an4, :);
% mask_an5 = all(facelabel_an5 > 0, 2);
% facelabel_an5 = facelabel_an5(mask_an5, :);
% tri_node_an5 = tri_node_an5(mask_an5, :);
% tri_normal_an5 = tri_normal_an5(mask_an5, :);
% % tri_centroid_an4 = h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceCentroids')';
% % tri_centroid_an5 = h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceCentroids')';
% tri_centroid_an4 = tri_centroid_an4(mask_an4, :);
% tri_centroid_an5 = tri_centroid_an5(mask_an5, :);
% 
% % ##### ID of one-piece faces #####
% face_to_calc = (1:length(tracked_uniqueface_an4))';
% face_to_calc(face_piecewise) = [];
% migration = zeros(length(tracked_uniqueface_an4), 1);


% ##### Data for Current Face #####
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/181025_mig_input.mat');
% X_to_Y_onepiece = X_to_Y(face_to_calc);
% trackedface_an4_onepiece = tracked_uniqueface_an4(face_to_calc, :);
% trackedface_an5_onepiece = tracked_uniqueface_an5(face_to_calc, :);
% mig_svm_proj = zeros(length(face_to_calc), 3);
% mig_normal_proj = zeros(length(face_to_calc), 4);

eps = 0.2;
% parfor i = 1:length(mig_svm_proj)
% [idx, ~] = find(face_to_calc == 2945);
idx = 1;

%     x_to_y = X_to_Y_onepiece{i};
%     obj_facelabel_an4 = trackedface_an4_onepiece(i, :);
%     obj_facelabel_an5 = trackedface_an5_onepiece(i, :);
    obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);

    % """
    % Note it's not legal to apply mask_half1 and mask_half2 seperately,
    % otherwise the data order would be broken.
    % """
    mask_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2)) ...
         | (facelabel_an4(:,1) == obj_facelabel_an4(2) & facelabel_an4(:,2) == obj_facelabel_an4(1));
    mask_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2)) ... 
         | (facelabel_an5(:,1) == obj_facelabel_an5(2) & facelabel_an5(:,2) == obj_facelabel_an5(1));

    facetri_node_an4 = tri_node_an4(mask_an4, :);
    facetri_node_an5 = tri_node_an5(mask_an5, :);
    facetri_normal_an4 = tri_normal_an4(mask_an4, :);
    facetri_normal_an5 = tri_normal_an5(mask_an5, :);
    facenode_id_an4 = unique(facetri_node_an4);
    facenode_id_an5 = unique(facetri_node_an5);
    facenode_coord_an4 = node_coord_an4(facenode_id_an4, :);
    facenode_coord_an5 = node_coord_an5(facenode_id_an5, :);
    % ----- note triangle normal direction has to be consistent -----
    mask_reverse_an4 = facelabel_an4(mask_an4,1) > facelabel_an4(mask_an4,2);
    mask_reverse_an5 = facelabel_an5(mask_an5,1) > facelabel_an5(mask_an5,2);
    facetri_normal_an4(mask_reverse_an4, :) = - facetri_normal_an4(mask_reverse_an4, :);
    facetri_normal_an5(mask_reverse_an5, :) = - facetri_normal_an5(mask_reverse_an5, :);

    coord_1 = facenode_coord_an4;
    coord_2 = facenode_coord_an5;

    % ##### calculate normal of a node as the average of all its resident triangles' nomral #####
    % """
    % Here it's done every time for a face. It's possible to calc for all nodes in one
    % time. However, in that case, TLNs and QNs can't have a normal direction.
    % """
    normal_an4 = zeros(length(coord_1), 3);
    normal_an5 = zeros(length(coord_2), 3);
    for j = 1:length(facenode_id_an4)
        mask = any(facetri_node_an4 == facenode_id_an4(j), 2);
        normal_tmp = facetri_normal_an4(mask, :);
        if size(normal_tmp,1) > 1
            normal_tmp = sum(normal_tmp);
        end
        normal_an4(j,:) = normal_tmp/norm(normal_tmp);
    end
    for j = 1:length(facenode_id_an5)
        mask = any(facetri_node_an5 == facenode_id_an5(j), 2);
        normal_tmp = facetri_normal_an5(mask, :);
        if size(normal_tmp,1) > 1
            normal_tmp = sum(normal_tmp);
        end
        normal_an5(j,:) = normal_tmp/norm(normal_tmp);
    end
    
    % ##### x_to_y  #####
%     m = size(coord_1, 1);
%     n = size(coord_2, 1);
%     features = [[ones(m, 1); 2*ones(n, 1)], [face_node_id_an4; face_node_id_an5], ...
%         [face_node_coord_an4; face_node_coord_an5], [normal_an4; normal_an5], ones(m+n, 1)];
% %     mig_svm_proj(i, :) = calcFaceMigBySVMPlaneProj(features, x_to_y, eps)
% %     mig_normal_proj(i, :) = calcFaceMigByLocalNormProj(features, x_to_y)
    % ##### corresp_localid #####
    % """
    % See main_TrackNodes.m for the making of corresp_localid. Basically,
    % nodes_an4 are kept the way they are, nodes_an5 are formulated as
    % correspondence of nodes_an4.
    % """
    id_use_an4 = corresp_localid{idx, 1}(:,1);
    id_use_an5 = corresp_localid{idx, 1}(:,2);
    m = length(id_use_an4);
    n = length(id_use_an5);
    features = [[ones(m, 1); 2*ones(n, 1)], [facenode_id_an4(id_use_an4); facenode_id_an5(id_use_an5)], ...
        [facenode_coord_an4(id_use_an4, :); facenode_coord_an5(id_use_an5, :)], [normal_an4(id_use_an4, :); normal_an5(id_use_an5, :)], ones(m+n, 1)];
    mig_svm_proj = calcFaceMigBySVMPlaneProj(features, eps);
    mig_normal_proj = calcFaceMigByLocalNormProj(features);


%     disp(i);
%     [mig_svm, features] = calcFaceMigBySVMPlaneProj(features, x_to_y, eps);

end



% ############################# Visualization #############################
% figure
% plotSVMPlane(features, face_tri_node_an4, face_tri_node_an5, x_to_y);

% face_node_info = getSingleFaceNodes(tracked_uniqueface_an4(face_to_calc(idx),:), tracked_uniqueface_an5(face_to_calc(idx),:));
% visualizeFace(face_node_info, x_to_y)
% unique(model.BoxConstraints)
% plotSingleFaceWithNormal(file_an4, obj_facelabel_an4, 0);
% hold on
% quiver3(coord_1(:,1),coord_1(:,2),coord_1(:,3), ...
%      normal_1(:,1),normal_1(:,2),normal_1(:,3),2,'color','g');
