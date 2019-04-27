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
%     - num_neigh_face and neigh_list_face are similar to that of grains. 
%     - USE FULL TL for finding connections
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

% 
% %% ####################################### Check #######################################
% % """
% % Note the index in nn_faces_ (featurefacelabel) is different from what is being used here (faces).
% % As long as the returned neighbor_ids are the same, the result is correct. 
% % """
% load('/Volumes/XIAOTING/Ni/working/190425_Hsmooth_energygrad_an5crop.mat', ...
%     'nn_faces_an4', 'nn_faces_an5', 'tracked_uniqueface_an4', 'tracked_uniqueface_an5');
% tracked_unique_faces = tracked_uniqueface_an4;
% nn_faces = nn_faces_an4;
% idx = 4215;
% 
% [r, ~] = find(faces(:,1) == tracked_unique_faces(idx, 1) & faces(:,2) == tracked_unique_faces(idx, 2));
% neighbors = getNeighList(r, num_neigh_face, neigh_list_faceid);
% 
% connection_fl = [nn_faces{idx, 1}(:, 1:2); nn_faces{idx, 2}(:, 1:2)]
% faces(neighbors, :)




