% % ############################################################################
% * Find faces with large local motion
%     - See if large local motion always happen near triple lines
%     - The definition of large local motion is currently by standard deviation
% % ############################################################################
dist = cell(length(small_face_id), 1);
dist_std = zeros(length(dist), 1);
for i = 1:length(small_face_id)
    idx = small_face_id(i);

    % ##### get triangles on the objective face #####
    obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
    corresp_1to2 = X_to_Y{idx};

    obj_facelabel_an4 = sort(obj_facelabel_an4, 2);
    facelabel_1 = sort(facelabel_1, 2);
    mask_objface_1 = (facelabel_1(:,1) == obj_facelabel_an4(1) & facelabel_1(:,2) == obj_facelabel_an4(2));
    obj_facelabel_an5 = sort(obj_facelabel_an5, 2);
    facelabel_2 = sort(facelabel_2, 2);
    mask_objface_2 = (facelabel_2(:,1) == obj_facelabel_an5(1) & facelabel_2(:,2) == obj_facelabel_an5(2));

    % ##### get id and coord of the objective nodes from the triangles #####
    node_id_face1 = tri_node_1(mask_objface_1, :);
    node_id_face1 = unique(node_id_face1);
    node_coord_face1 = node_coord_1(node_id_face1,:);
    node_id_face2 = tri_node_2(mask_objface_2, :);
    node_id_face2 = unique(node_id_face2);
    node_coord_face2 = node_coord_2(node_id_face2,:);
    node_coord_2for1 = node_coord_face2(corresp_1to2, :);

    dist{i} = vecnorm(node_coord_face1 - node_coord_2for1, 2, 2);
    dist_std(i) = std(dist{i});
    disp(num2str(i));
end

face_largelocalmotion_1 = [];
for i = 1:length(small_face_id)
    if dist_std(i) > 5 && length(X_to_Y{i}) > 30 && length(Y_to_X{i}) > 30
        face_largelocalmotion_1 = [face_largelocalmotion_1; i];
    end
end



%%
face_largelocalmotion_2 = [];
for i = 1:length(dist)
    face_dist = dist{i};
    if sum(face_dist < 5) > length(face_dist)/2
        avg_smalldist = sum(face_dist(face_dist < 5)) / sum(face_dist < 5);
        if sum(face_dist > 3*avg_smalldist) > 10
            face_largelocalmotion_2 = [face_largelocalmotion_2; i];
        end
    end
end

        
%%
% % ############################################################################
% Get Mobile Face Info for Paraview
% % ############################################################################
%% ##### Load Data ##### 
file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
% file_an4 = '/Volumes/XIAOTING/Ni_180823/An4new6_fixedOrigin_smooth_forParaview.dream3d';
% file_an5 = '/Volumes/XIAOTING/Ni_180823/An5new6_smooth_forParaview.dream3d';
num_neigh_an4 = h5read(file_an4,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors').';
num_neigh_an5 = h5read(file_an5,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors').';
num_neigh_an4(1) = [];
num_neigh_an5(1) = [];
neigh_list_an4 = h5read(file_an4,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList');
neigh_list_an5 = h5read(file_an5,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList');
face_id_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FeatureFaceId')';
face_label_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
face_id_an4 = face_id_an4(all(face_label_an4 > 0, 2), :);
face_label_an4 = face_label_an4(all(face_label_an4 > 0, 2), :);
face_id_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceData/FeatureFaceId')';
face_label_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
face_id_an5 = face_id_an5(all(face_label_an5 > 0, 2), :);
face_label_an5 = face_label_an5(all(face_label_an5 > 0, 2), :);
face_label_an4 = sort(face_label_an4, 2);
face_label_an5 = sort(face_label_an5, 2);

load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/180822.mat', 'tracked_uniqueface_an4', 'tracked_uniqueface_an5');
% obj_face = [2723, 654, 1500, 2158, 3102, 3832, 2945]';

%% ##### Individual Face Info ##### 
clc
idx = 6710;

obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);

neighbor_an4_1 = getNeighList(obj_facelabel_an4(1), num_neigh_an4, neigh_list_an4);
neighbor_an4_2 = getNeighList(obj_facelabel_an4(2), num_neigh_an4, neigh_list_an4);
neighbor_an5_1 = getNeighList(obj_facelabel_an5(1), num_neigh_an5, neigh_list_an5);
neighbor_an5_2 = getNeighList(obj_facelabel_an5(2), num_neigh_an5, neigh_list_an5);

neighbor_an4 = intersect(neighbor_an4_1, neighbor_an4_2);
neighbor_an5 = intersect(neighbor_an5_1, neighbor_an5_2);

mask_an4 = ((face_label_an4(:,1) == obj_facelabel_an4(1) & face_label_an4(:,2) == obj_facelabel_an4(2)));
faceid_an4 = unique(face_id_an4(mask_an4));
mask_an5 = ((face_label_an5(:,1) == obj_facelabel_an5(1) & face_label_an5(:,2) == obj_facelabel_an5(2)));
faceid_an5 = unique(face_id_an5(mask_an5));

disp(['face id in an4 =  ', num2str(faceid_an4)]);
disp(['face label in an4: [', num2str(obj_facelabel_an4), ']']);
disp(['common neighbors in an4: [', num2str(neighbor_an4'), ']']);
disp(['face id in an5 =  ', num2str(faceid_an5)]);
disp(['face label in an5: [', num2str(obj_facelabel_an5), ']']);
disp(['common neighbors in an5: [', num2str(neighbor_an5'), ']']);
disp(' ');
