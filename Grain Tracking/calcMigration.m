% function calcMigration(file_an4, file_an4, X_to_Y)
% ############################################################################
% Notes
%     - This function is doing data preparation.
%     - Migration of each grain face is calced by calling
%     calcFaceMigBySVMPlaneProj.m | calcFaceMigByPillarHeight.m |
%     calcFaceMigByLocalNormProj.m
% ############################################################################
% ----------------------- load debug data -----------------------
% load('180822', 'file_an4', 'file_an5');
file_an4 = '/Volumes/XIAOTING/Ni/An4new6_fixedOrigin_smooth.dream3d';
file_an5 = '/Volumes/XIAOTING/Ni/An5new6_smooth.dream3d';
load('180822_FaceCorresp', 'tracked_uniqueface_an4', 'tracked_uniqueface_an5', 'X_to_Y');
% ---------------------------------------------------------------

% ##### Load Data & Basic Cleans #####
facelabel_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
facelabel_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_normal_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
tri_normal_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
tri_node_an4 = 1 + double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
tri_node_an5 = 1 + double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
node_coord_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
mask_an4 = all(facelabel_an4 > 0, 2);
facelabel_an4 = facelabel_an4(mask_an4, :);
tri_node_an4 = tri_node_an4(mask_an4, :);
tri_normal_an4 = tri_normal_an4(mask_an4, :);
mask_an5 = all(facelabel_an5 > 0, 2);
facelabel_an5 = facelabel_an5(mask_an5, :);
tri_node_an5 = tri_node_an5(mask_an5, :);
tri_normal_an5 = tri_normal_an5(mask_an5, :);
tri_centroid_an4 = h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceCentroids')';
tri_centroid_an5 = h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceCentroids')';
tri_centroid_an4 = tri_centroid_an4(mask_an4, :);
tri_centroid_an5 = tri_centroid_an5(mask_an5, :);


% ##### Data for Current Face #####
idx = 1674;
x_to_y = X_to_Y{idx};
obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);

% """
% Note it's not legal to apply mask_half1 and mask_half2 seperately,
% otherwise the data order would be broken.
% """
mask_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2)) ...
     | (facelabel_an4(:,1) == obj_facelabel_an4(2) & facelabel_an4(:,2) == obj_facelabel_an4(1));
mask_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2)) ... 
     | (facelabel_an5(:,1) == obj_facelabel_an5(2) & facelabel_an5(:,2) == obj_facelabel_an5(1));

face_tri_node_an4 = tri_node_an4(mask_an4, :);
face_tri_node_an5 = tri_node_an5(mask_an5, :);
face_node_id_an4 = unique(face_tri_node_an4);
face_node_id_an5 = unique(face_tri_node_an5);
face_node_coord_an4 = node_coord_an4(face_node_id_an4, :);
face_node_coord_an5 = node_coord_an5(face_node_id_an5, :);
face_tri_normal_an4 = tri_normal_an4(mask_an4, :);
face_tri_normal_an5 = tri_normal_an5(mask_an5, :);
face_tri_centroid_an4 = tri_centroid_an4(mask_an4, :);
face_tri_centroid_an5 = tri_centroid_an5(mask_an5, :);
% ----- note triangle normal direction has to be consistent -----
mask_reverse_an4 = facelabel_an4(mask_an4,1) > facelabel_an4(mask_an4,2);
mask_reverse_an5 = facelabel_an5(mask_an5,1) > facelabel_an5(mask_an5,2);
face_tri_normal_an4(mask_reverse_an4, :) = - face_tri_normal_an4(mask_reverse_an4, :);
face_tri_normal_an5(mask_reverse_an5, :) = - face_tri_normal_an5(mask_reverse_an5, :);

coord_1 = face_node_coord_an4;
coord_2 = face_node_coord_an5;
% ##### calculate normal of a node as the average of all its resident triangles' nomral #####
% """
% When calculating normal directions associated with the nodes, here it's
% done every time for a face. It's possible to calc for all nodes in one
% time. However, in that case, TLNs and QNs can't have a normal direction.
% """
normal_1 = zeros(length(coord_1), 3);
normal_2 = zeros(length(coord_2), 3);
for i = 1:length(face_node_id_an4)
    mask = any(face_tri_node_an4 == face_node_id_an4(i), 2);
    normal_tmp = face_tri_normal_an4(mask, :);
    if size(normal_tmp,1) > 1
        normal_tmp = sum(normal_tmp);
    end
    normal_1(i,:) = normal_tmp/norm(normal_tmp);
end
for i = 1:length(face_node_id_an5)
    mask = any(face_tri_node_an5 == face_node_id_an5(i), 2);
    normal_tmp = face_tri_normal_an5(mask, :);
    if size(normal_tmp,1) > 1
        normal_tmp = sum(normal_tmp);
    end
    normal_2(i,:) = normal_tmp/norm(normal_tmp);
end








% ############################# Visualization #############################
face_node_info = getSingleFaceNodes(obj_facelabel_an4, obj_facelabel_an5);
figure
visualizeFace(face_node_info, x_to_y)
% plotSingleFaceWithNormal(file_an4, obj_facelabel_an4, 0);
% hold on
% quiver3(coord_1(:,1),coord_1(:,2),coord_1(:,3), ...
%      normal_1(:,1),normal_1(:,2),normal_1(:,3),2,'color','g');
