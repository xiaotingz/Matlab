file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixedOrigin_smooth.dream3d';
file_an5 = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4.dream3d';
load('/Volumes/XIAOTING/Ni/working/Grain Tracking/look_up_table_an4_an5crop.mat')
% file_an4 = '/Volumes/XIAOTING/Ni/simu_An4_clean_seg.dream3d';
% file_an5 = '/Volumes/XIAOTING/Ni/simu_An5_clean_seg.dream3d';

eps_curv = 1;
eps_area = 7;
eps_min_ang = 10;

faces_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
faces_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
faces_an4 = faces_an4(any(faces_an4>0, 2), :);
faces_an5 = faces_an5(any(faces_an5>0, 2), :);
faces_an4 = sortrows(sort(faces_an4, 2));
faces_an5 = sortrows(sort(faces_an5, 2));


% [result_an4, ~, ~, ~] = findQuadNodes(file_an4);
% [triple_line_an4, tl_info_an4] = findTripleLines(file_an4, eps_area, eps_curv, eps_min_ang);
% [triple_line_full_an4, tl_info_full_an4] = findTripleLines(file_an4, 10e6, 10e6, 0);
[num_corners_an4, num_edges_an4] = getFaceCharacter(faces_an4, triple_line_full_an4, result_an4{1}, result_an4{2});
[num_nnface_avgcorner_an4] = findFaceNNAvgCorner(file_an4, faces_an4, num_corners_an4, triple_line_full_an4);

% [result_an5, ~, ~, ~] = findQuadNodes(file_an5);
% [triple_line_an5, tl_info_an5] = findTripleLines(file_an5, eps_area, eps_curv, eps_min_ang);
% [triple_line_full_an5, tl_info_full_an5] = findTripleLines(file_an5, 10e6, 10e6, 0);
[num_corners_an5, num_edges_an5] = getFaceCharacter(faces_an5, triple_line_full_an5, result_an5{1}, result_an5{2});
[num_nnface_avgcorner_an5] = findFaceNNAvgCorner(file_an5, faces_an5, num_corners_an5, triple_line_full_an5);

[avg_nng_diff, max_nng_diff] = findFaceLocalTopologyChange(file_an4, file_an5, faces_an4, look_up_table, triple_line_full_an4);

% 
% [da_len_weighted_an4, da_num_weighted_an4] = calcGrainFaceDAs(faces_an4, triple_line_an4, tl_info_an4);
% [da_len_weighted_an5, da_num_weighted_an5] = calcGrainFaceDAs(faces_an5, triple_line_an5, tl_info_an5);

%% ##################################### Checks #####################################
% """
% num_corners & num_nnface_avgcorner
% """
file = file_an4;
num_corners = num_corners_an4;
num_edges = num_edges_an4;
faces = faces_an4;
num_nnface_avgcorner = num_nnface_avgcorner_an4;
tl = triple_line_full_an4;
face_feature_label = h5read(file, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';

[~, num_neigh_face, neigh_list_faceid] = findFaceConnection(file, tl);
keys_face_feature_label = num2str(face_feature_label);
keys_face_feature_label = mat2cell(keys_face_feature_label, ones(1, size(keys_face_feature_label, 1)), size(keys_face_feature_label, 2));
idx_face_feature_label = (1:length(face_feature_label))' - 1;
faceid_dict = containers.Map(keys_face_feature_label, idx_face_feature_label);
keys_faces = num2str(faces);

idx = randi(length(faces));
disp(['idx_full_faces = ', num2str(idx), ';   idx_face_feature_id = ', num2str(faceid_dict(keys_faces(idx, :)))])
disp(['num_corners = ', num2str(num_corners(idx)), ',  num_edges = ', num2str(num_edges(idx))])

neighbors = getNeighList(idx, num_neigh_face, neigh_list_faceid);
neighbor_corners = num_corners(neighbors);
disp(['num_nnface_avgcorner_calc = ', num2str(num_nnface_avgcorner(idx)), ...
        ';   check = ', num2str(sum(neighbor_corners)/length(neighbor_corners))]);

mask_tl = (sum(ismember(tl, faces(idx, :)), 2) == 2);
resident_tl = tl(mask_tl, :);
disp('resident triple lines:')
disp(resident_tl)






