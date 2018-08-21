function visualizeFace(obj_face, obj_point_smooth_an4, obj_point_smooth_an5, x_to_y, alpha_value)
% ############################################################################
% Input
%     - obj_face = {4, 2}, wrote by getSingleFaceNodes.m
%     - x_to_y = [m, 1], wrote by solveNodeCorresp.m
%     - obj_point_... = [n, 3], surface mesh node coordinates from D3D
% NOTES
%     - change inputs according to what is desired. 
% ############################################################################
% ----------------------- load debug data -----------------------
% file_coordsmooth_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
% file_coordsmooth_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
% file_1 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_forParaview.dream3d');
% file_2 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh_forParaview.dream3d');

% obj_face = face_node_info_3;
% obj_point_id_an4 = obj_face{2,1}([1, 64, 472]);
% obj_point_smooth_an4 = node_coordsmooth_an4(obj_point_id_an4, :);
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
trisurf(reshape(obj_face{2,1}, 3, [])', node_coordsmooth_an4(:,1), node_coordsmooth_an4(:,2), node_coordsmooth_an4(:,3),'Facecolor',color1, 'Facealpha', alpha_value, 'edgealpha', alpha_value);
hold on
trisurf(reshape(obj_face{2,2}, 3, [])', node_coordsmooth_an5(:,1), node_coordsmooth_an5(:,2), node_coordsmooth_an5(:,3),'Facecolor',color2, 'Facealpha', alpha_value, 'edgealpha', alpha_value);
rotate3d on
scatter3(obj_point_smooth_an4(:,1), obj_point_smooth_an4(:,2), obj_point_smooth_an4(:,3), 20, 'filled', 'MarkerFaceColor',color1, 'MarkerEdgeColor',color1);
scatter3(obj_point_smooth_an5(:,1), obj_point_smooth_an5(:,2), obj_point_smooth_an5(:,3), 20, 'filled', 'MarkerFaceColor',color2, 'MarkerEdgeColor',color2);

if nargin == 5
    for i = 1:length(x_to_y)
        plot3([obj_point_smooth_an4(i,1), obj_point_smooth_an5(x_to_y(i),1)], [obj_point_smooth_an4(i,2), obj_point_smooth_an5(x_to_y(i),2)], [obj_point_smooth_an4(i,3), obj_point_smooth_an5(x_to_y(i),3)], 'k', 'LineWidth', 1);
%         text(obj_point_smooth_an4(i,1), obj_point_smooth_an4(i,2), obj_point_smooth_an4(i,3),[' ', num2str(int32(i))],'FontSize',16, 'Color', color1);
    end
end



end








