function plotSingleFaceWithNormal(file, obj_facelabel, reverse_winding)
% ###########################################################################
% * Inputs
%     - reverse_winding: a bool variable determining the triangle normal
%     direction
% * Notes
%     - can be compared with plotSVMPlane.m
% ###########################################################################
% ------------------ load data for debug --------------------
% file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
% obj_facelabel = [332, 2275];
% reverse_winding = 1;
% -----------------------------------------------------------
facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
normal = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals')).';
centroid = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceCentroids')).';
tri_node = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';

mask_half1 = (facelabel(:,1) == obj_facelabel(1) & facelabel(:,2) == obj_facelabel(2));
mask_half2 = (facelabel(:,1) == obj_facelabel(2) & facelabel(:,2) == obj_facelabel(1));
mask_full = (mask_half1 | mask_half2);

obj_node = tri_node(mask_full, :);
obj_centroid = centroid(mask_full, :);
obj_normal = normal(mask_full, :);

% --- adjust the winding of triangle normal in display ---
obj_label = facelabel(mask_full, :);
mask_adjust = obj_label(:,1) > obj_label(:,2);
obj_normal(mask_adjust, :) = - obj_normal(mask_adjust, :);
if reverse_winding == 1
    obj_normal = - obj_normal;
end

trimesh(obj_node, node_coord(:,1), node_coord(:,2), node_coord(:,3));
rotate3d
hold on
quiver3(obj_centroid(:,1),obj_centroid(:,2),obj_centroid(:,3), ...
     obj_normal(:,1),obj_normal(:,2),obj_normal(:,3),1.5,'color','r');
hold off 


%% ##### Calc Triangle Normal by MATLAB #####
% mesh_matlab = triangulation(obj_node, [node_coord(:,1), node_coord(:,2), node_coord(:,3)]);
% tri_center_matlab = incenter(mesh_matlab);
% tri_normal_matlab = faceNormal(mesh_matlab);
% tri_normal_matlab(mask_adjust, :) = - tri_normal_matlab(mask_adjust, :);
% if reverse_winding == 1
%     tri_normal_matlab = - tri_normal_matlab;
% end
% 
% figure(2)
% trimesh(mesh_matlab)
% rotate3d on
% hold on  
% quiver3(tri_center_matlab(:,1),tri_center_matlab(:,2),tri_center_matlab(:,3), ...
%      tri_normal_matlab(:,1),tri_normal_matlab(:,2),tri_normal_matlab(:,3),0.5,'color','r');
% hold off
end






