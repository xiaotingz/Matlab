% % ############################################################################
% % Solve Matching
% % * Notes
% %     - In this file, the coordinates on face are not stored because the
% %     memory cost is too expensive. However, the face coordinates can be obtained by
% %     running getSingleFaceNodes.m
%       - 
% % ############################################################################
% facelabel_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
% facelabel_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
%     % -- NOTE triNodes are indexes starting from zero 
% tri_node_an4 = 1 + double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
% tri_node_an5 = 1 + double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
% node_coord_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
% node_coord_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
% 
% % --- tracked_unique_face = trackUniqueFace(file_an4, file_an5, look_up_table, complete) ---
% [tracked_uniqueface_an4, tracked_uniqueface_an5] = trackUniqueFace(file_an4, file_an5, look_up_table, 1);
% n = length(tracked_uniqueface_an4);

%%
facelabel_an4 = sort(facelabel_an4, 2);
facelabel_an5 = sort(facelabel_an5, 2);

X_to_Y = cell(n,1);
Y_to_X = cell(n,1);
large_face_id = [];
for i = 1:length(tracked_uniqueface_an4)
    obj_facelabel_an4 = tracked_uniqueface_an4(i, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(i, :);
    
    % ##### get objective triangles on the objective face #####
    mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
    mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));

    % ##### get id and coord of the objective triangles #####
    face_node_id_an4 = tri_node_an4(mask_objface_an4, :);
    face_node_id_an4 = unique(face_node_id_an4);
    face_node_coord_an4 = node_coord_an4(face_node_id_an4,:);
    face_node_id_an5 = tri_node_an5(mask_objface_an5, :);
    face_node_id_an5 = unique(face_node_id_an5);
    face_node_coord_an5 = node_coord_an5(face_node_id_an5,:);
    
    if (length(face_node_id_an4) < 1000) && (length(face_node_id_an5) < 1000)
        [x_to_y, y_to_x] = solveNodeCorresp(face_node_coord_an4, face_node_coord_an5);
        X_to_Y{i} = int32(x_to_y);
        Y_to_X{i} = int32(y_to_x)';
    else 
        large_face_id = [large_face_id, i];
        X_to_Y{i} = NaN;
        Y_to_X{i} = NaN;
    end
    disp(['small faces: ', num2str(i), '/7004']);
end
save('180820_smallFaces.mat', 'X_to_Y', 'Y_to_X', 'large_face_id');


huge_face_id = [];
large_face_size = length(large_face_id);
X_to_Y_large = cell(large_face_size,1);
Y_to_X_large= cell(large_face_size,1);
parfor i = 1:large_face_size
    idx = large_face_id(i);
    obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
    
    % ##### get objective triangles on the objective face #####
    mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
    mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));

    % ##### get id and coord of the objective triangles #####
    face_node_id_an4 = tri_node_an4(mask_objface_an4, :);
    face_node_id_an4 = unique(face_node_id_an4);
    face_node_coord_an4 = node_coord_an4(face_node_id_an4,:);
    face_node_id_an5 = tri_node_an5(mask_objface_an5, :);
    face_node_id_an5 = unique(face_node_id_an5);
    face_node_coord_an5 = node_coord_an5(face_node_id_an5,:);
    
    if length(face_node_coord_an4) < 10000 && length(face_node_coord_an5) < 10000
        [x_to_y, y_to_x] = solveNodeCorresp(face_node_coord_an4, face_node_coord_an5);
        X_to_Y_large{i} = int32(x_to_y);
        Y_to_X_large{i} = int32(y_to_x);
    else 
        huge_face_id = [huge_face_id, idx];
        X_to_Y_large{i} = NaN;
        Y_to_X_large{i} = NaN;
    end
    
    disp(['large faces: ', num2str(i), '/', num2str(large_face_sizee)]);
end
save('180820_largeFaces.mat', 'X_to_Y_large', 'Y_to_X_large', 'huge_face_id');



%% ##### Visualize Face Correspondences #####
idx = small_face_id(randi(length(small_face_id)))
% idx = 3102;

% ----- get the object face triangles and nodes -----
obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
face_node_info = getSingleFaceNodes(file_an4, obj_facelabel_an4, file_an5, obj_facelabel_an5);

% ----- plot with alpha shape-----
figure(1)
visualizeFace(face_node_info, X_to_Y{idx})

% ----- alpha shape -----
% shp = alphaShape([face_node_info{4,1}; face_node_info{4,2}]);
% shp.Alpha = 16.4938;
% plot(shp)
% alpha(.5)
% rotate3d on
% volume = shp.volume;
