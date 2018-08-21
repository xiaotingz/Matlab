function visualizeFace(obj_face, x_to_y, obj_point_smooth_an4, obj_point_smooth_an5)
% ############################################################################
% Input
%     - obj_face = {4, 2}, wrote by getSingleFaceNodes.m
%     - x_to_y = [m, 1], wrote by solveNodeCorresp.m
%     - obj_point_... = [n, 3], surface mesh node coordinates from D3D
% NOTES
%     - change inputs according to what is desired. 
% ############################################################################
% ----------------------- load debug data -----------------------
% obj_face = face_node_info;
% obj_point_smooth_an4 = face_node_info{4,1};
% obj_point_smooth_an5 = face_node_info{4,2};
% x_to_y = X_to_Y{idx};
% alpha_value = 0.3;
% ---------------------------------------------------------------
load('node_coord');
color1 = [0, 0.4470, 0.7410];
color2 = [0.9290, 0.6940, 0.1250];

% figure(1)
% trisurf(reshape(obj_face{2,1}, 3, [])', node_coordmesh_an4(:,1), node_coordmesh_an4(:,2), node_coordmesh_an4(:,3),'Facecolor',color1, 'Facealpha', 0.3, 'edgealpha', 0.3 );
% hold on
% trisurf(reshape(obj_face{2,2}, 3, [])', node_coordmesh_an5(:,1), node_coordmesh_an5(:,2), node_coordmesh_an5(:,3),'Facecolor',color2, 'Facealpha', 0.3, 'edgealpha', 0.3);
% rotate3d on

% figure(2)
trisurf(reshape(obj_face{2,1}, 3, [])', node_coordsmooth_an4(:,1), node_coordsmooth_an4(:,2), node_coordsmooth_an4(:,3),'Facecolor',color1, 'Facealpha', alpha_value, 'edgealpha', 0.3);
hold on
trisurf(reshape(obj_face{2,2}, 3, [])', node_coordsmooth_an5(:,1), node_coordsmooth_an5(:,2), node_coordsmooth_an5(:,3),'Facecolor',color2, 'Facealpha', alpha_value, 'edgealpha', 0.3);
rotate3d on

obj_point_smooth_an4 = unique(obj_point_smooth_an4, 'rows');
obj_point_smooth_an5 = unique(obj_point_smooth_an5, 'rows');

scatter3(obj_point_smooth_an4(:,1), obj_point_smooth_an4(:,2), obj_point_smooth_an4(:,3), 20, 'filled', 'MarkerFaceColor',color1, 'MarkerEdgeColor',color1);
scatter3(obj_point_smooth_an5(:,1), obj_point_smooth_an5(:,2), obj_point_smooth_an5(:,3), 20, 'filled', 'MarkerFaceColor',color2, 'MarkerEdgeColor',color2);
if nargin == 2
    for i = 1:length(obj_point_smooth_an4)
        plot3([obj_point_smooth_an4(i,1), obj_point_smooth_an5(x_to_y(i),1)], [obj_point_smooth_an4(i,2), obj_point_smooth_an5(x_to_y(i),2)], [obj_point_smooth_an4(i,3), obj_point_smooth_an5(x_to_y(i),3)], 'k', 'LineWidth', 1);
%         text(obj_point_smooth_an4(i,1), obj_point_smooth_an4(i,2), obj_point_smooth_an4(i,3),[' ', num2str(int32(i))],'FontSize',16, 'Color', color1);
    end
elseif nargin == 4
    
end

end








