function nn_faces = getLeftRightFaceConnections(tracked_faces, tls, file, faces_in_use)
% ##########################################################################
% * Input
%     - unique_faces = [n1,2]  
%           returned by trackFace.m or trackUniqueFace.m, usually named tracked_uniqueface_ in other functions
%     - tls = [m, 3]
%           returned by Topologies/findTripleLines.m
%     - faces_in_use
%           'use_complete_faces'
%           'use_all_true_faces'
% * Output
%     - nn_faces = {{k1, 3}*n, {k2, 3}*n}
%           [face_label_1, face_label_2, feature_face_id]
%           k1 is #left_connected_faces, k2 is #right_connected_faces, of the ith grain face. 
%           There are in total n grain faces.
% * NOTE
%     - This function is to calculate left/right connection of grain 
%       faces. The idea is related to that of high/low energy (big/small size) grains, 
%       but applied on the level of grain faces. The objective is to predict migration sign.
%     - Related functions: calcFaceToGrainCentroidDist.m
% ##########################################################################
% % ----------------------- load debug data -----------------------
% % file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
% % file_an5 = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4_Hsmooth.dream3d';
% % tracked_faces = faces;
% % tls = triple_line_full_an4;
% % faces_in_use = 'use_inner_faces';
% % --- tls are returned by findTripleLines.m in /Topology
% % """ don't threshold trianlge quality for finding face connection, threshold quality only for diheral angle """
% % [triple_line_full_an4, tl_info_full_an4] = findTripleLines(file_an4, 1e6, 1e6, 0);
% tls = triple_line_full_an5;
% % clear tracked_uniqueface_an4 tl_an4
% % ---------------------------------------------------------------
feature_face_label = h5read(file, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
feature_face_id = (1:size(feature_face_label, 1))' - 1;
if strcmp(faces_in_use, 'use_inner_faces')
    inner_faces = feature_face_label(all(feature_face_label>0, 2), :);
    inner_ids = feature_face_id(all(feature_face_label>0, 2));
elseif strcmp(faces_in_use, 'use_all_true_faces')
    inner_faces = feature_face_label(any(feature_face_label>0, 2), :);
    inner_ids = feature_face_id(any(feature_face_label>0, 2));
else
    warning('Please specify faces of interest: ''use_inner_faces'' or ''use_all_true_faces''')
end

nn_faces = cell(size(tracked_faces, 1), 2);
for i = 1:size(tracked_faces)
    left_grain = tracked_faces(i, 1);
    right_grain = tracked_faces(i, 2);
    mask_tls = (sum((tls == left_grain | tls == right_grain), 2) == 2);
    candidate_tl = tls(mask_tls, :);
    candidate_grains = setdiff(unique(candidate_tl), [left_grain, right_grain]);
    
    label_l_tmp = [left_grain * int32(ones(size(candidate_grains))), candidate_grains];
    label_l_tmp = sort(label_l_tmp, 2);
    mask = ismember(inner_faces, label_l_tmp, 'rows');
    nn_faces{i, 1} = [inner_faces(mask, :), inner_ids(mask)];
    
    label_r_tmp = [right_grain * int32(ones(size(candidate_grains))), candidate_grains];
    label_r_tmp = sort(label_r_tmp, 2);
    mask = ismember(inner_faces, label_r_tmp, 'rows');
    nn_faces{i, 2} = [inner_faces(mask, :), inner_ids(mask)];
    
end

end
%%
% nn_faces_an5 = nn_faces;
% load('../190425.mat')
% rng('shuffle');
% idx = randi(length(tracked_faces));
% mask = (inner_faces(:, 1)==min(tracked_faces(idx, :)) & inner_faces(:, 2)==max(tracked_faces(idx, :)));
% 
% disp('###################################################')
% disp(['pair_', num2str(idx)])
% disp(['face_label = ', num2str(tracked_faces(idx, :)), ',     face_id = ', num2str(inner_ids(mask))]);
% disp('previous left:       ')
% disp(nn_faces_an4{idx, 1})
% disp('current left:       ')
% disp(nn_faces{idx, 1})
% disp('------')
% disp('previous right:       ')
% disp(nn_faces_an4{idx, 2})
% disp('current right:       ')
% disp(nn_faces{idx, 2})







