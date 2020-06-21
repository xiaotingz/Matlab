function [num_corners, num_edges] = getFaceCharacter(faces, TLs, QNs, FCNs)
% ############################################################################
% * Inputs
%   - faces = [n, 2]
%   - TLs, given by findTripleLines.m
%   - QNs and FCNs, given by findQuadNodes.m result{1} and result{2}
% * Note
%   - This function takes only FiveCoordSuperQN but not other superQNs.
%         - The assumption behind every FiveCoordNode is counted as two QN
%           is that the formation of five grain set is likely due to identifiability, or resolution problem. 
%         - In that case, there should be three grains with two corners and
%           two grains with one corner, (E.g. [1,2,3,4,5] --> [1,2,3,4] & [2,3,4,5])
%         - We don't know which are the true two QN pairs, so we just assign every grain two QNs.
% ############################################################################
% ----------------------- load debug data -----------------------
% A = [277, 277, 277]';
% B = [8, 69, 100]';
% labels = [277, 8];
% ---------------------------------------------------------------
num_edges = zeros(size(faces, 1), 1);
num_corners = zeros(size(faces, 1), 1);

% ##### if only QuadPoints #####
if nargin == 3
    for i = 1:size(faces, 1)
        A = faces(i, 1);
        B = faces(i, 2);

        mask_TL_A = (TLs == A);
        mask_TL_B = (TLs == B);
        mask_TL_AB = ((sum(mask_TL_A, 2) + sum(mask_TL_B, 2)) == 2);
        num_edges(i) = sum(mask_TL_AB);

        mask_QN_A = (QNs == A);
        mask_QN_B = (QNs == B);
        mask_QN_AB = ((sum(mask_QN_A, 2) + sum(mask_QN_B, 2)) == 2);
        num_corners(i) = sum(mask_QN_AB);
    end
    
% ##### QuadPoints and FCNoints #####
elseif nargin == 4
    for i = 1:size(faces, 1)
        A = faces(i, 1);
        B = faces(i, 2);

        mask_TL_A = (TLs == A);
        mask_TL_B = (TLs == B);
        mask_TL_AB = ((sum(mask_TL_A, 2) + sum(mask_TL_B, 2)) == 2);
        num_edges(i) = sum(mask_TL_AB);

        mask_QN_A = (QNs == A);
        mask_QN_B = (QNs == B);
        mask_QN_AB = ((sum(mask_QN_A, 2) + sum(mask_QN_B, 2)) == 2);
        
        mask_FCN_A = (FCNs == A);
        mask_FCN_B = (FCNs == B);
        mask_FCN_AB = ((sum(mask_FCN_A, 2) + sum(mask_FCN_B, 2)) == 2);
        
        num_corners(i) = sum(mask_QN_AB) + sum(mask_FCN_AB)*2;
    end
end

% mask_TL = zeros(length(TLs),2);
% for i = 1:length(TLs)
%     mask_TL(i, :) = ismember([A, B], TLs(i,:));
% end
% numEdges = sum(sum(mask_TL,2) == 2)
% mask_QN = zeros(length(generalQNs),2);
% for i = 1:length(generalQNs)
%     mask_QN(i, :) = ismember([A, B], generalQNs(i,:));
% end
% numCorners = sum(sum(mask_QN,2) == 2)
    
% end


end


% ################################## Checks ##################################
% % ------------------- Info for paraview -------------------
% % """
% % num_corners & num_nnface_avgcorner
% % """
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190730_Lsmoth_Topologies.mat', ...
%     'file_an4', 'file_an5', 'faces_an4', 'faces_an5', 'triple_line_full_an4', 'triple_line_full_an5', ...
%     'num_corners_an4', 'num_corners_an5', 'num_edges_an4', 'num_edges_an5');
% tls_an4 = triple_line_full_an4;
% tls_an5 = triple_line_full_an5;
% face_feature_label_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
% face_feature_label_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
% face_feature_label_an4(1,:) = [];
% face_feature_label_an5(1,:) = [];
% 
% % ----- Prepare dictionary for facefeatureid -----
% keys_face_feature_label = num2str(face_feature_label_an4);
% keys_face_feature_label = mat2cell(keys_face_feature_label, ones(1, size(keys_face_feature_label, 1)), size(keys_face_feature_label, 2));
% idx_face_feature_label = (1:length(face_feature_label_an4))';
% faceid_dict_an4 = containers.Map(keys_face_feature_label, idx_face_feature_label);
% keys_faces_an4 = num2str(faces_an4);
% keys_face_feature_label = num2str(face_feature_label_an5);
% keys_face_feature_label = mat2cell(keys_face_feature_label, ones(1, size(keys_face_feature_label, 1)), size(keys_face_feature_label, 2));
% idx_face_feature_label = (1:length(face_feature_label_an5))';
% faceid_dict_an5 = containers.Map(keys_face_feature_label, idx_face_feature_label);
% keys_faces_an5 = num2str(faces_an5);
% 
% 
% %% ----- Display info -----
% idx = randi(length(faces_an4));
% 
% disp('------------')
% disp(['label_an4 = ', num2str(faces_an4(idx, :)), ';   face_feature_id_an4 = ', num2str(faceid_dict_an4(keys_faces_an4(idx, :)))])
% disp(['num_corners_an4 = ', num2str(num_corners_an4(idx)), ',  num_edges_an4 = ', num2str(num_edges_an4(idx))])
% % disp(['label_an5 = ', num2str(faces_an5(idx, :)), ';   face_feature_id_an5 = ', num2str(faceid_dict_an5(keys_faces_an5(idx, :)))])
% % disp(['num_corners_an5 = ', num2str(num_corners_an5(idx)), ',  num_edges_an5 = ', num2str(num_edges_an5(idx))])
% 
% mask_tl = (sum(ismember(tls_an4, faces_an4(idx, :)), 2) == 2);
% resident_tl = tls_an4(mask_tl, :);
% disp('resident triple lines:')
% disp(resident_tl)
% 
% 
