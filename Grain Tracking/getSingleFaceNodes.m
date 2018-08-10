function face_node_info = getSingleFaceNodes(file_1, obj_facelabel_1, file_2, obj_facelabel_2)
% ##########################################################################
% * Input
%     - facelabe = [1,2], order ([A, B] or [B, A]) doesn't matter 
% * Output
%     - face_nodes = {3, 1} 
%         1st dimension: facelabel
%         2nd dimension = [n*3, 1]: id of the nodes on the face, row first order ([tri1_3nodes, tri_3nodes, ..]) 
%         3rd dimension = [n*3, 3]: coordinate of the nodes on the face
% * NOTE
%     - make sure samples are aligned by changing ORIGIN.
% ##########################################################################
% ----------------------- load debug data -----------------------
% file = file_an4;
% obj_facelabel = [620,  752];
% ---------------------------------------------------------------
facelabel_1 = double(h5read(file_1,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
facelabel_2 = double(h5read(file_2,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
    % -- NOTE triNodes are indexes starting from zero 
tri_node_1 = 1 + double(h5read(file_1,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
tri_node_2 = 1 + double(h5read(file_2,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_1 = double(h5read(file_1,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
node_coord_2 = double(h5read(file_2,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
node_type_1 = double(h5read(file_1,'/DataContainers/TriangleDataContainer/VertexData/NodeType'))';
node_type_2 = double(h5read(file_2,'/DataContainers/TriangleDataContainer/VertexData/NodeType'))';

% ##### filter bad data #####
mask_1 = all(facelabel_1 > 0, 2);
facelabel_1 = facelabel_1(mask_1,:);
tri_node_1 = tri_node_1(mask_1,:);
mask_2 = all(facelabel_2 > 0, 2);
facelabel_2 = facelabel_2(mask_2,:);
tri_node_2 = tri_node_2(mask_2,:);

% ##### get objective triangles on the objective face #####
obj_facelabel_1 = sort(obj_facelabel_1, 2);
facelabel_1 = sort(facelabel_1, 2);
mask_objface_1 = (facelabel_1(:,1) == obj_facelabel_1(1) & facelabel_1(:,2) == obj_facelabel_1(2));
obj_facelabel_2 = sort(obj_facelabel_2, 2);
facelabel_2 = sort(facelabel_2, 2);
mask_objface_2 = (facelabel_2(:,1) == obj_facelabel_2(1) & facelabel_2(:,2) == obj_facelabel_2(2));

% ##### get id and coord of the objective triangles #####
face_node_id_1 = tri_node_1(mask_objface_1, :);
face_node_id_1 = reshape(face_node_id_1', 1, [])';
face_node_coord_1 = node_coord_1(face_node_id_1,:);
face_node_id_2 = tri_node_2(mask_objface_2, :);
face_node_id_2 = reshape(face_node_id_2', 1, [])';
face_node_coord_2 = node_coord_2(face_node_id_2,:);

% ##### get node types of the objective triangles #####
face_node_type_1 = node_type_1(face_node_id_1, :);
face_node_type_2 = node_type_1(face_node_id_2, :);

% ##### write data #####
face_node_info = cell(4,2);
face_node_info{1,1} = obj_facelabel_1;
face_node_info{2,1} = face_node_id_1;
face_node_info{3,1} = face_node_type_1;
face_node_info{4,1} = face_node_coord_1;
face_node_info{1,2} = obj_facelabel_2;
face_node_info{2,2} = face_node_id_2;
face_node_info{3,2} = face_node_type_2;
face_node_info{4,2} = face_node_coord_2;
end