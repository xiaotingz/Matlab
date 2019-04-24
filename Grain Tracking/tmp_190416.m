%% ###################### Calculate face centroid of interest ##########################
load('/Volumes/XIAOTING/Ni/working/Grain Tracking/tmp_190416.mat', ...
    'file_an4', 'file_an5', 'tracked_uniqueface_an4_sort', 'tracked_uniqueface_an5');
% ##### Read and Clean Data #####
fl_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
tri_node_an4 = 1 + double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
fl_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
tri_node_an5 = 1 + double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';

tri_node_an4 = tri_node_an4(all(fl_an4 > 0, 2), :);
fl_an4 = fl_an4(all(fl_an4 > 0, 2), :);
fl_an4 = sort(fl_an4, 2);
tri_node_an5 = tri_node_an5(all(fl_an5 > 0, 2), :);
fl_an5 = fl_an5(all(fl_an5 > 0, 2), :);
fl_an5 = sort(fl_an5, 2);

% ##### Id for faces of interest #####
obj_idx = [4484, 5623];
obj_face_an4 = tracked_uniqueface_an4_sort(obj_idx, :);
obj_face_an5 = tracked_uniqueface_an5(obj_idx, :);

% ##### Find face_centroid #####
face_centroid_an4 = zeros(size(obj_face_an4, 1), 3);
face_centroid_an5 = zeros(size(obj_face_an4, 1), 3);
for i = 1:size(obj_face_an4)
    mask_an4 = (fl_an4(:,1) == obj_face_an4(i, 1) & fl_an4(:,2) == obj_face_an4(i, 2) ...
             | fl_an4(:,1) == obj_face_an4(i, 2) & fl_an4(:,2) == obj_face_an4(i, 1));
    nodes_an4 = unique(tri_node_an4(mask_an4, :));
    face_centroid_an4(i, :) = sum(node_coord_an4(nodes_an4, :))/length(nodes_an4);
    mask_an5 = (fl_an5(:,1) == obj_face_an5(i, 1) & fl_an5(:,2) == obj_face_an5(i, 2) ...
             | fl_an5(:,1) == obj_face_an5(i, 2) & fl_an5(:,2) == obj_face_an5(i, 1));
    nodes_an5 = unique(tri_node_an5(mask_an5, :));
    face_centroid_an5(i, :) = sum(node_coord_an5(nodes_an5, :))/length(nodes_an5);
end

diff = face_centroid_an5 - face_centroid_an4;

rfvec_twin = ([1,1,1]/norm([1,1,1]))*tand(60/2);
eps = 0.05;
mask = vecnorm(rfvecs_an4 - rfvec_twin, 2, 2)<eps & vecnorm(rfvecs_an5 - rfvec_twin, 2, 2)<eps;

idx_grains = (1:size(dist_ctwin, 1))';
idx_twins = idx_grains(mask);


%% ###################### Plot Twins ##########################
rng('shuffle');
i = randi(length(idx_twins));
idx = idx_twins(i);


face_node_info = getSingleFaceNodes(tracked_uniqueface_an4(idx,:), tracked_uniqueface_an5(idx,:));

figure
visualizeFace(face_node_info)
disp(['pair ', num2str(idx), ',   dist_ctwin=', num2str(dist_ctwin(idx))]);


f_name = ['twins_pair_', num2str(idx)];
% f_name = ['twins_pair_', num2str(idx_twins(idx)), '_2'];
print(f_name, '-dtiff', '-r300')















