% % ############################################################################
% * Key data
%     - mig_svm_proj = [n, 3]; 
%         [mig_svm_sign_avg, mig_svm_abs_avg, bad_fit];
%     - mig_normal_proj = [n, 4];
%         [mig_1_sign_avg, mig_2_sign_avg, mig_1_abs_avg, mig_2_abs_avg];
%     - mig_
% % ############################################################################

load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/180822.mat', 'tracked_uniqueface_an4', 'tracked_uniqueface_an5');
load('/Volumes/XIAOTING/Ni/190328/180822_FaceCorresp.mat', 'X_to_Y', 'Y_to_X')
load('/Volumes/XIAOTING/Ni/190328/181107_mig_piececorresp_comb.mat', 'face_onepiece')
file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');

% ###################################### Data Preparation ######################################
facelabel_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
facelabel_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_normal_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
tri_normal_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
tri_curv_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'))';
tri_curv_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'))';
tri_area_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'))';
tri_area_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'))';
    % -- NOTE triNodes are indexes starting from zero --
tri_node_an4 = 1 + double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
tri_node_an5 = 1 + double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
node_coord_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
% ---- filter bad data -----
mask_an4 = all(facelabel_an4 > 0, 2);
mask_an5 = all(facelabel_an5 > 0, 2);
facelabel_an4 = facelabel_an4(mask_an4, :);
facelabel_an5 = facelabel_an5(mask_an5, :);
tri_node_an4 = tri_node_an4(mask_an4, :);
tri_node_an5 = tri_node_an5(mask_an5, :);
tri_normal_an4 = tri_normal_an4(mask_an4, :);
tri_normal_an5 = tri_normal_an5(mask_an5, :);
tri_curv_an4 = tri_curv_an4(mask_an4, :);
tri_curv_an5 = tri_curv_an5(mask_an5, :);
tri_area_an4 = tri_area_an4(mask_an4, :);
tri_area_an5 = tri_area_an5(mask_an5, :);

n = length(tracked_uniqueface_an4);
corresp_localid = cell(n,1);
corresp_globalid = cell(n,1);

clear mask_an4 mask_an5
%%
% load('/Volumes/XIAOTING/Ni/190328/181107.mat')
% % X_to_Y = X_to_Y_nearest; clear X_to_Y_nearest;
% load('/Volumes/XIAOTING/Ni/190328/180822_FaceCorresp.mat', 'X_to_Y', 'Y_to_X')
% load('/Volumes/XIAOTING/Ni/190328/181107_mig_piececorresp_comb.mat', 'face_onepiece')
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/181107.mat');
% load('/Users/xiaotingzhong/Documents/Matlab/Grain Tracking/data/180822_FaceCorresp.mat', 'X_to_Y', 'Y_to_X');
% load('/Users/xiaotingzhong/Documents/Matlab/Grain Tracking/data/181107_mig_piececorresp_comb.mat', 'face_onepiece');
eps = 0.2;

% tracked_uniqueface_an4 = tracked_uniqueface_an4(face_onepiece,:);
% tracked_uniqueface_an5 = tracked_uniqueface_an5(face_onepiece,:);
% X_to_Y = X_to_Y(face_onepiece);
% Y_to_X = Y_to_X(face_onepiece);


n = length(tracked_uniqueface_an4);
mig_localnorm_proj = zeros(n, 4);
mig_svm_proj_naivekmeans = zeros(n, 2);
mig_svm_proj_improve = zeros(n, 2);
mig_lr_proj = zeros(n, 2);
mig_pillar_proj = zeros(n, 3);

% ###################################### Distance Projection ######################################
% parfor i = 1:size(tracked_uniqueface_an4,1)
    i = 886;
    disp(['working on pair ', num2str(i)]);
    obj_facelabel_an4 = tracked_uniqueface_an4(i, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(i, :);

    % ##### Get Objective Triangles on the Objective Face #####
    mask_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2) |...
        facelabel_an4(:,1) == obj_facelabel_an4(2) & facelabel_an4(:,2) == obj_facelabel_an4(1));
    mask_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2) |...
        facelabel_an5(:,1) == obj_facelabel_an5(2) & facelabel_an5(:,2) == obj_facelabel_an5(1));
    mask_reverse_an4 = facelabel_an4(mask_an4,1) > facelabel_an4(mask_an4,2);
    mask_reverse_an5 = facelabel_an5(mask_an5,1) > facelabel_an5(mask_an5,2);
    

    % ##### Get Id and Coord of the Objective Triangles #####
    facetri_nodeid_an4 = tri_node_an4(mask_an4, :);
    facetri_normal_an4 = tri_normal_an4(mask_an4, :);
    facetri_normal_an4(mask_reverse_an4, :) = - facetri_normal_an4(mask_reverse_an4, :);
    facetri_curv_an4 = tri_curv_an4(mask_an4, :);
    facetri_curv_an4(mask_reverse_an4, :) = - facetri_curv_an4(mask_reverse_an4, :);
    facetri_area_an4 = tri_area_an4(mask_an4, :);
    facenode_id_an4 = unique(facetri_nodeid_an4);
    facenode_coord_an4 = node_coord_an4(facenode_id_an4,:);
    facetri_nodeid_an5 = tri_node_an5(mask_an5, :);
    facetri_normal_an5 = tri_normal_an5(mask_an5, :);
    facetri_normal_an5(mask_reverse_an5, :) = - facetri_normal_an5(mask_reverse_an5, :);
    facetri_curv_an5 = tri_curv_an5(mask_an5, :);
    facetri_curv_an5(mask_reverse_an5, :) = - facetri_curv_an5(mask_reverse_an5, :);
    facetri_area_an5 = tri_area_an5(mask_an5, :);
    facenode_id_an5 = unique(facetri_nodeid_an5);
    facenode_coord_an5 = node_coord_an5(facenode_id_an5,:);

    % ---------------------------------------- If One-piece ---------------------------------------- 
    % """
    % Here it's done every time for a face. It's possible to calc for all nodes in one
    % time. However, in that case, TLNs and QNs can't have a normal direction.
    % """
%     if length(unique(subgraph_id_an4)) == 1 && length(unique(subgraph_id_an5)) == 1
        % ----- Solve node corresp -----
        x_to_y = X_to_Y{i};
        y_to_x = Y_to_X{i};
%         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        subgraph_nodeid_local_an4 = (1:length(facenode_id_an4))';
        subgraph_nodeid_local_an5 = (1:length(facenode_id_an5))';
%         !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        corresp_localid{i} = int32([subgraph_nodeid_local_an4; subgraph_nodeid_local_an5]);
        corresp_globalid{i} = int32([facenode_id_an4; facenode_id_an5]);
        
        % ----- Calculate normal of a node as the average of all its resident triangles' nomral -----
        normal_an4 = zeros(length(facenode_coord_an4), 3);
        curv_an4 =  zeros(length(facenode_coord_an4), 1);
        normal_an5 = zeros(length(facenode_coord_an5), 3);
        curv_an5 =  zeros(length(facenode_coord_an5), 1);
        for j = 1:length(facenode_id_an4)
            mask = any(facetri_nodeid_an4 == facenode_id_an4(j), 2);
            normal_tmp = facetri_normal_an4(mask, :);
            if size(normal_tmp,1) > 1
                normal_tmp = sum(normal_tmp);
            end
            normal_an4(j,:) = normal_tmp/norm(normal_tmp);
            curv_an4(j) = sum(facetri_curv_an4(mask, :)) / sum(mask);
        end
        for j = 1:length(facenode_id_an5)
            mask = any(facetri_nodeid_an5 == facenode_id_an5(j), 2);
            normal_tmp = facetri_normal_an5(mask, :);
            if size(normal_tmp,1) > 1
                normal_tmp = sum(normal_tmp);
            end
            normal_an5(j,:) = normal_tmp/norm(normal_tmp);
            curv_an5(j) = sum(facetri_curv_an5(mask, :)) / sum(mask);
        end
        
        % ----- Calculate migration distance -----
        m = length(facenode_id_an4);
        n = length(facenode_id_an5);
        features = [[ones(m, 1); 2*ones(n, 1)], [facenode_id_an4; facenode_id_an5], ...
            [facenode_coord_an4; facenode_coord_an5], [normal_an4; normal_an5], [curv_an4; curv_an5], ones(m+n, 1)];
        %%
        [mig_1_abs, mig_2_abs, mig_1_sign, mig_2_sign] = calcFaceMigByLocalNormProj(features, x_to_y);
        mig_localnorm_proj(i, :) = [mig_1_abs, mig_2_abs, mig_1_sign, mig_2_sign];
        [mig_pillar_abs, mig_pillar_sign, fit_ratio] = calcFaceMigByPillarHeight(features, 10, x_to_y, facetri_nodeid_an4, y_to_x, facetri_nodeid_an5);
        mig_pillar_proj(i, :) = [mig_pillar_abs, mig_pillar_sign, fit_ratio];
        [mig_svm_abs_naive, mig_svm_sign_naive, features_svm_naive, svm_plane_info_naive] = calcFaceMigBySVMPlaneProj(features, 0.2, x_to_y);
        mig_svm_proj_naivekmeans(i, :) =  [mig_svm_abs_naive, mig_svm_sign_naive];
        [mig_svm_abs, mig_svm_sign, features_svm, svm_plane_info] = calcFaceMigBySVMPlaneProj(features, 0.2, x_to_y, facetri_nodeid_an4, facetri_nodeid_an5);
        mig_svm_proj_improve(i, :) =  [mig_svm_abs, mig_svm_sign];
        [mig_lr_abs, mig_lr_sign, features_lr, lr_plane_info] = calcFaceMigByLinearRegPlaneProj(features, 0.2, x_to_y);
        mig_lr_proj(i, :) = [mig_lr_abs, mig_lr_sign];
        
        disp('migration_abs: ');
        disp(['norm_1 = ', num2str(mig_1_abs), ', norm_2 = ',  num2str(mig_2_abs), ', svm_naive = ', num2str(mig_svm_abs_naive), ...
            ', svm_improve = ', num2str(mig_svm_abs), ', lr = ', num2str(mig_lr_abs), ', pillar = ', num2str(mig_pillar_abs)]);
        disp('migration_sign: ');
        disp(['norm_1 = ', num2str(mig_1_sign), ', norm_2 = ',  num2str(mig_2_sign), ', svm_naive = ', num2str(mig_svm_sign_naive), ...
            ', svm_improve = ', num2str(mig_svm_sign), ', lr = ', num2str(mig_lr_sign), ', pillar = ', num2str(mig_pillar_sign)]);
        
%     % ---------------------------------------- If Multiple-pieces ---------------------------------------- 
%     else
%         % ----- Solve the corresp of subgraphs -----
%         subgraph_corresp = solveSubgraphCorrespOverlap(subgraph_id_an4,subgraph_id_an5, ...
%             facetri_normal_an4, facetri_normal_an5, facetri_nodeid_an4, facetri_nodeid_an5, facenode_id_an4, facenode_id_an5, ...
%             facenode_coord_an4, facenode_coord_an5, facetri_area_an4, facetri_area_an5);
% 
%          % ----- calculate normal of a node to prepare migration calculation -----
%         normal_an4 = zeros(length(facenode_coord_an4), 3);
%         normal_an5 = zeros(length(facenode_coord_an5), 3);
%         for j = 1:length(facenode_id_an4)
%             mask = any(facetri_nodeid_an4 == facenode_id_an4(j), 2);
%             normal_tmp = facetri_normal_an4(mask, :);
%             if size(normal_tmp,1) > 1
%                 normal_tmp = sum(normal_tmp);
%             end
%             normal_an4(j,:) = normal_tmp/norm(normal_tmp);
%         end
%         for j = 1:length(facenode_id_an5)
%             mask = any(facetri_nodeid_an5 == facenode_id_an5(j), 2);
%             normal_tmp = facetri_normal_an5(mask, :);
%             if size(normal_tmp,1) > 1
%                 normal_tmp = sum(normal_tmp);
%             end
%             normal_an5(j,:) = normal_tmp/norm(normal_tmp);
%         end
%         % ----- prepare vector holders -----
%         mig_svm_sign = [];  mig_svm_abs = [];  bad_fit = false;
%         mig_1_sign = [];  mig_2_sign =[];  mig_1_abs = [];  mig_2_abs = [];
%         
%         % ----- Because there are 1-N corresps, need to first find the group of subgraphs -----
%         useful_subgraph_an4 = unique(subgraph_corresp(:,1));
%         useful_subgraph_an5 = unique(subgraph_corresp(:,2));
%         [num_group, idx_fewerpiece] = min([length(useful_subgraph_an4), length(useful_subgraph_an5)]);
%         X_to_Y_tmp = cell(num_group, 1);
%         corresp_localid_tmp = cell(num_group, 1);
%         corresp_globalid_tmp = cell(num_group, 1);
%         for j = 1:num_group
%             if idx_fewerpiece == 1
%                 group_subgraphid_an4 = useful_subgraph_an4(j);
%                 group_subgraphid_an5 = subgraph_corresp(subgraph_corresp(:,1)==group_subgraphid_an4, 2);
%                 mask_group_nodes_an4 = (subgraph_id_an4 == group_subgraphid_an4);
%                 mask_group_nodes_an5 = ismember(subgraph_id_an5, group_subgraphid_an5);
%                 subgraph_nodeid_locali_an4 = subgraph_nodeid_local_an4(mask_group_nodes_an4);
%                 subgraph_nodeid_locali_an5 = subgraph_nodeid_local_an5(mask_group_nodes_an5);
%                 subgraph_nodeid_an4 = facenode_id_an4(mask_group_nodes_an4, :);
%                 subgraph_nodeid_an5 = facenode_id_an5(mask_group_nodes_an5, :);
%                 subgraph_coord_an4 = facenode_coord_an4(mask_group_nodes_an4, :);
%                 subgraph_coord_an5 = facenode_coord_an5(mask_group_nodes_an5, :);
%             else
%                 group_subgraphid_an5 = useful_subgraph_an5(j);
%                 group_subgraphid_an4 = subgraph_corresp(subgraph_corresp(:,2)==group_subgraphid_an5, 1);
%                 mask_group_nodes_an5 = (subgraph_id_an5 == group_subgraphid_an5);
%                 mask_group_nodes_an4 = ismember(subgraph_id_an4, group_subgraphid_an4);
%                 subgraph_nodeid_locali_an4 = subgraph_nodeid_local_an4(mask_group_nodes_an4);
%                 subgraph_nodeid_locali_an5 = subgraph_nodeid_local_an5(mask_group_nodes_an5);
%                 subgraph_nodeid_an4 = facenode_id_an4(mask_group_nodes_an4, :);
%                 subgraph_nodeid_an5 = facenode_id_an5(mask_group_nodes_an5, :);
%                 subgraph_coord_an4 = facenode_coord_an4(mask_group_nodes_an4, :);
%                 subgraph_coord_an5 = facenode_coord_an5(mask_group_nodes_an5, :);
%             end
%             
%             % ----- Find node corresp between the subgraph groups -----
%             [x_to_y, y_to_x] = solveNodeCorresp(subgraph_coord_an4, subgraph_coord_an5);
%             X_to_Y_tmp{j} = x_to_y;
%             corresp_localid_tmp{j} = int32([subgraph_nodeid_locali_an4, subgraph_nodeid_locali_an5(x_to_y)]);
%             corresp_globalid_tmp{j} = int32([subgraph_nodeid_an4, subgraph_nodeid_an5(x_to_y)]);
%             
%             % ----- Calculate migration distance -----
%             id_use_an4 = subgraph_nodeid_locali_an4;
%             id_use_an5 = subgraph_nodeid_locali_an5(x_to_y);
%             m = length(id_use_an4);
%             n = length(id_use_an5);
%             features = [[ones(m, 1); 2*ones(n, 1)], [facenode_id_an4(id_use_an4); facenode_id_an5(id_use_an5)], ...
%                 [facenode_coord_an4(id_use_an4, :); facenode_coord_an5(id_use_an5, :)], [normal_an4(id_use_an4, :); normal_an5(id_use_an5, :)], ones(m+n, 1)];
%             [mig_svm_sign_i, mig_svm_abs_i, bad_fit_i, ~] = calcFaceMigBySVMPlaneProj(features, eps);
%             [mig_1_sign_i, mig_2_sign_i, mig_1_abs_i, mig_2_abs_i] = calcFaceMigByLocalNormProj(features);
% 
%             mig_svm_sign = [mig_svm_sign; mig_svm_sign_i];
%             mig_svm_abs = [mig_svm_abs; mig_svm_abs_i];
%             bad_fit = (bad_fit || bad_fit_i);
%             mig_1_sign = [mig_1_sign; mig_1_sign_i];
%             mig_2_sign = [mig_2_sign; mig_2_sign_i];
%             mig_1_abs = [mig_1_abs; mig_1_abs_i];
%             mig_2_abs = [mig_2_abs; mig_2_abs_i];
%         end
%         
%         % ----- average migration for the piece-wise face -----
%         X_to_Y{i} = X_to_Y_tmp;
%         corresp_localid{i} = corresp_localid_tmp;
%         corresp_globalid{i} = corresp_globalid_tmp;
%         
%         % ----- average migration for the piece-wise face -----
%         mig_svm_sign_avg = sum(mig_svm_sign)/length(mig_svm_sign);
%         mig_svm_abs_avg = sum(mig_svm_abs)/length(mig_svm_abs);
%         mig_1_sign_avg = sum(mig_1_sign)/length(mig_1_sign);
%         mig_2_sign_avg = sum(mig_2_sign)/length(mig_2_sign);
%         mig_1_abs_avg = sum(abs(mig_1_abs))/length(mig_1_abs);
%         mig_2_abs_avg = sum(abs(mig_2_abs))/length(mig_2_abs);
%         mig_svm_proj(i, :) = [mig_svm_sign_avg, mig_svm_abs_avg, bad_fit];
%         mig_normal_proj(i, :)  = [mig_1_sign_avg, mig_2_sign_avg, mig_1_abs_avg, mig_2_abs_avg];
% 
%         % ----- Record id of the piece_wise faces-----
% %         face_piecewise = [face_piecewise, i, max(subgraph_id_an4), max(subgraph_id_an5)];
%     end
    
% end

% save('190316.mat', 'mig_svm', 'mig_lr');
