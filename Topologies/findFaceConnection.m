function [faces, num_neigh_face, neigh_list_faceid] = findFaceConnection(file, tl)
% ############################################################################
% Input 
%     - tl = [m, 3]
%         returned by findTripleLines.m
% Output
%     - faces = [n, 2]
%         face_labels
%     - num_neigh_face = [n, 1]
%     - neigh_list_faceid = [sum(num_neigh), 1]
%         Id of the neighbors, given with respect to faces
% Notes
%     - Neighbor of a face is qualified as sharing a common edge
%           num_neigh_face and neigh_list_face are similar to that of grains. 
%           both the tracked and untracked neighbors are listed. 
%     - USE FULL TL for finding connections
%     - Similar function: GrainTracking/getLeftRightFaceConnections.m
% ############################################################################
% % ----------------------- load debug data -----------------------
% file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
% file_an5 = '/Volumes/XIAOTING/Ni/An5new6_cropToAn4_Hsmooth.dream3d';
% load('/Volumes/XIAOTING/Ni/working/190425_Hsmooth_geo_topo_an5crop.mat', 'triple_line_full_an4', 'triple_line_full_an5');
% file = file_an4;
% tl = triple_line_full_an4;
% % ---------------------------------------------------------------

%  ##### Get Labels of the Inner Grain Faces #####
faces = h5read(file, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
faces = faces(all(faces>0, 2), :);
faces = sortrows(faces);

% ##### Calculate num_neigh_face and neigh_list_face #####
num_neigh_face = zeros(length(faces), 1);
neigh_list_facelabel = [];
for i = 1:length(faces)
    mask_tl_contain_face = (sum(ismember(tl, faces(i,:)), 2) == 2);
    candidate_grain = setdiff(unique(tl(mask_tl_contain_face, :)), faces(i, :));
    tmp = int32(ones(length(candidate_grain), 1));
    candidate_neigh = [faces(i, 1)*tmp, candidate_grain; tmp*faces(i, 2), candidate_grain];
    candidate_neigh = sort(candidate_neigh, 2);
    mask_neigh = ismember(candidate_neigh, faces, 'rows');
    num_neigh_face(i) = sum(mask_neigh);
    neigh_list_facelabel = [neigh_list_facelabel; candidate_neigh(mask_neigh, :)];
end


% ##### Convert Elements of neigh_list_face, From Face Labels Into ids in faces  #####
% ----- Make a map by numerical hash -----
% faces_hash = int64(1e5 * faces(:,1)) +  int64(faces(:,2));
% neigh_list_face_hash = int64(1e5 * neigh_list_facelabel(:,1)) + int64(neigh_list_facelabel(:,2));
% 
% faceid_dict = containers.Map(faces_hash, (1:length(faces_hash))');
% 
% neigh_list_faceid = zeros(size(neigh_list_face_hash));
% for i = 1:length(neigh_list_facelabel)
%     neigh_list_faceid(i) = faceid_dict(neigh_list_face_hash(i));
% end

% ----- Make a map by convert key to string -----
keys_faces = num2str(faces);
keys_faces = mat2cell(keys_faces, ones(1, size(keys_faces, 1)), size(keys_faces, 2));
faceid_dict = containers.Map(keys_faces, (1:length(faces))');
num_digits = length(num2str(max(faces(:))));
if length(keys_faces{1,1}) ~= num_digits*2 + 2
    warning('key (from face_labels) format length wrong in findFaceConnection!')
end
key_format = ['%', num2str(num_digits), 'd  %', num2str(num_digits),'d'];
neigh_list_faceid = zeros(size(neigh_list_facelabel, 1), 1);
for i = 1:length(neigh_list_facelabel)
    key = sprintf(key_format, neigh_list_facelabel(i, 1), neigh_list_facelabel(i, 2));
    neigh_list_faceid(i) = faceid_dict(key);
end

end

% %% ####################################### Check #######################################
% % """ See 190425_190712_connected_faces 190712_info.txt for detail records """
% % ---------- Load Data ----------
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/records/190425_190712_connected_faces/190712_data_an4an5.mat')
% % """
% % file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
% % file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4.dream3d';
% % load('look_up_table_an4_an5crop.mat')
% % load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190624_Ni_crop_TopologyResult_uniqiue.mat', ...
% %     'triple_line_full_an4')
% % tls = triple_line_full_an4;
% % [tracked_uniqueface_an4, tracked_uniqueface_an5] = trackUniqueFace(file_an4, file_an5, look_up_table, 'none');
% % [faces, num_neigh_face, neigh_list_faceid] = findFaceConnection(file_an4, tls);
% % nn_faces = getLeftRightFaceConnections(tracked_uniqueface_an4, tls, file_an4);
% % """
% 
% % ---------- Add faces_feature_face_id_d3d for faces returned by findFaceConnection ----------
% feature_face_label = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
% feature_face_label(1, :) = [];
% feature_face_id = (1:size(feature_face_label, 1))';
% mask = all(feature_face_label>0, 2);
% feature_face_id = feature_face_id(mask);
% feature_face_label = feature_face_label(mask, :);
% feature_face_label = sort(feature_face_label, 2);
% dict_facefeatureid = containers.Map('KeyType','char','ValueType','int32'); 
% for i = 1:size(feature_face_label, 1)
%     key_tmp = mat2str(feature_face_label(i, :));
%     dict_facefeatureid(key_tmp) = feature_face_id(i);
% end
% faces_featurefaceid = zeros(size(faces, 1), 1);
% for i = 1:size(faces, 1)
%     tmp = dict_facefeatureid(mat2str(faces(i, :)));
%     faces_featurefaceid(i) = tmp;
% end
% 
% % ---------- Make a dictionray of connected faces for the tracked faces ----------
% n = size(tracked_uniqueface_an4, 1);
% dict_neighbor_leftright = containers.Map('KeyType','char','ValueType','any'); 
% for i = 1:n
%     key_faceid = mat2str(tracked_uniqueface_an4(i, :));
%     dict_neighbor_leftright(key_faceid) = [nn_faces{i, 1}; nn_faces{i, 2}];
% end
% 
% %%
% for i = 1:20
%     % ---------- Pick an Id from the full faces of an4 ----------
%     rng('shuffle')
%     idx = randi(size(faces, 1));
%     label = mat2str(faces(idx, :));
% 
%     % ----- If tracked, cross Check with /Grain Tracking/getLeftRightFaceConnections -----
%     disp('-------------------------------------')
%     disp(['face = ', label, '; face_feature_id = ', num2str(dict_facefeatureid(label))]);
%     id_neigh_face = getNeighList(idx, num_neigh_face, neigh_list_faceid);
%     res_neigh_list = [faces(id_neigh_face, :), faces_featurefaceid(id_neigh_face)];
%     if isKey(dict_neighbor_leftright, label)
%         res_leftright = dict_neighbor_leftright(label);
%         res_leftright = sortrows(res_leftright);
%         res_neigh_list = sortrows(res_neigh_list);
%         if ~ all(all(res_leftright == res_neigh_list))
%             warning('Neighbor results are different!')
%         else
%             disp('tracked and results matched, neighbors = ');
%             disp(res_neigh_list)
%         end
%         disp(' ')
%     % ----- If untracked, display info and check in D3D -----
%     else
%         disp('untracked, neighbors = ')
%         disp(res_neigh_list)
%     end
% 
% end



