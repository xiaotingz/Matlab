function [faces, num_neigh_face, neigh_list_faceid] = findGrainFaceConnection(file)
% ############################################################################
% Output
%     - faces = [n, 2]
%         face_labels
%     - num_neigh_face = [n, 1]
%     - neigh_list_faceid = [sum(num_neigh), 1]
%         Id of the neighbors, given with respect to faces
% Notes
%     - neighbor of a face is qualified as sharing a common edge
%     - num_neigh_face and neigh_list_face are similar to that of grains. 
% ############################################################################
% ----------------------- load debug data -----------------------
% file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
% % '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d'
% ---------------------------------------------------------------

%  ##### Get Labels of the Inner Grain Faces #####
faces = h5read(file, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
faces = faces(all(faces>0, 2), :);
faces = sortrows(faces);


%  ##### Get Triple Lines #####
tl = findTripleLines(file);

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
% ----- Convert the two-number facelabels into one single number -----
faces_hash = int64(1e5 * faces(:,1)) +  int64(faces(:,2));
neigh_list_face_hash = int64(1e5 * neigh_list_facelabel(:,1)) + int64(neigh_list_facelabel(:,2));

faceid_dict = containers.Map(faces_hash, (1:length(faces_hash))');

neigh_list_faceid = zeros(size(neigh_list_face_hash));
for i = 1:length(neigh_list_facelabel)
    neigh_list_faceid(i) = faceid_dict(neigh_list_face_hash(i));
end


% % ----- Test if neigh_list_faceid is correct -----
% test = faces(neigh_list_faceid, :);
% sum(neigh_list_facelabel == test);

end