% ############################################################################
% Solve Matching
% * Notes
%     - In this file, the coordinates on face are not stored because the
%     memory cost is too expensive. However, the face coordinates can be obtained by
%     running getSingleFaceNodes.m
% ############################################################################
facelabel_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
facelabel_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
    % -- NOTE triNodes are indexes starting from zero 
tri_node_an4 = 1 + double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
tri_node_an5 = 1 + double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
node_coord_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';

% --- tracked_unique_face = trackUniqueFace(file_an4, file_an5, look_up_table, complete) ---
[tracked_uniqueface_an4, tracked_uniqueface_an5] = trackUniqueFace(file_an4, file_an5, look_up_table, 1);
n = length(tracked_uniqueface_an4);

%%
% for i = 1:length(tracked_uniqueface_an4)
X_to_Y = cell(n,1);
Y_to_X = cell(n,1);
for i = 1:n
    obj_facelabel_an4 = tracked_uniqueface_1(i, :);
    obj_facelabel_an5 = tracked_uniqueface_2(i, :);
    
    % ##### get objective triangles on the objective face #####
    obj_facelabel_an4 = sort(obj_facelabel_an4, 2);
    facelabel_an4 = sort(facelabel_an4, 2);
    mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
    obj_facelabel_an5 = sort(obj_facelabel_an5, 2);
    facelabel_an5 = sort(facelabel_an5, 2);
    mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));

    % ##### get id and coord of the objective triangles #####
    face_node_id_an4 = tri_node_an4(mask_objface_an4, :);
    face_node_id_an4 = reshape(face_node_id_an4', 1, [])';
    face_node_coord_an4 = node_coord_an4(face_node_id_an4,:);
    face_node_id_an5 = tri_node_an5(mask_objface_an5, :);
    face_node_id_an5 = reshape(face_node_id_an5', 1, [])';
    face_node_coord_an5 = node_coord_an5(face_node_id_an5,:);
    
    [x_to_y, y_to_x] = solveNodeCorresp(face_node_coord_an4, face_node_coord_an5);
    X_to_Y{i} = int32(x_to_y);
    Y_to_X{i} = int32(y_to_x);
end


%% ##### Visualize Face Correspondences #####
idx = small_face(randi(length(small_face)));
obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);

face_node_info = getSingleFaceNodes(file_an4, obj_facelabel_an4, file_an5, obj_facelabel_an5);
figure(1)
visualizeFace(face_node_info, face_node_info{4,1}, face_node_info{4,2}, X_to_Y{idx})
