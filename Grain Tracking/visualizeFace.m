function visualizeFace(obj_face, x_to_y, obj_node_an4, obj_node_an5, obj)
% ############################################################################
% Input
%     - obj_face = {4, 2}, wrote by getSingleFaceNodes.m
%     - x_to_y = [m, 1], wrote by solveNodeCorresp.m
%     - obj_point_... = [n, 3], surface mesh node coordinates from D3D
%     - obj = 'distort_tri' || 'corresp' 
%           use with getDistortedMatchTriangles.m
% NOTES
%     - change inputs according to what is desired. 
%     - default is to plot the correspondence between all nodes, but
%     obj_node can be specified to a set of nodes
% ############################################################################
% ----------------------- load debug data -----------------------
% obj_face = face_node_info;
% obj_point_smooth_an4 = face_node_info{4,1};
% obj_point_smooth_an5 = face_node_info{4,2};
% x_to_y = X_to_Y{idx};
% ---------------------------------------------------------------
load('node_coord');
color1 = [0, 0.4470, 0.7410];
color2 = [0.9290, 0.6940, 0.1250];


% ##### Plot Mesh #####
trisurf(reshape(obj_face{2,1}, 3, [])', node_coordsmooth_an4(:,1), node_coordsmooth_an4(:,2), node_coordsmooth_an4(:,3),'Facecolor',color1, 'Facealpha', 0.3, 'edgealpha', 0.3);
hold on
trisurf(reshape(obj_face{2,2}, 3, [])', node_coordsmooth_an5(:,1), node_coordsmooth_an5(:,2), node_coordsmooth_an5(:,3),'Facecolor',color2, 'Facealpha', 0.3, 'edgealpha', 0.3);
rotate3d on


% ##### If No obj_node Specified, Default Will Use All Nodes On Surface #####
if nargin == 2
    node_id_an4 = unique(obj_face{2,1});
    node_id_an5 = unique(obj_face{2,2});
    obj_node_an4 = node_coordsmooth_an4(node_id_an4, :);
    obj_node_an5 = node_coordsmooth_an5(node_id_an5, :);
    scatter3(obj_node_an4(:,1), obj_node_an4(:,2), obj_node_an4(:,3), 20, 'filled', 'MarkerFaceColor',color1, 'MarkerEdgeColor',color1);
    scatter3(obj_node_an5(:,1), obj_node_an5(:,2), obj_node_an5(:,3), 20, 'filled', 'MarkerFaceColor',color2, 'MarkerEdgeColor',color2);
    for i = 1:length(x_to_y)
        plot3([obj_node_an4(i,1), obj_node_an5(x_to_y(i),1)], [obj_node_an4(i,2), obj_node_an5(x_to_y(i),2)], [obj_node_an4(i,3), obj_node_an5(x_to_y(i),3)], 'k', 'LineWidth', 1);
%         text(obj_point_smooth_an4(i,1), obj_point_smooth_an4(i,2), obj_point_smooth_an4(i,3),[' ', num2str(int32(i))],'FontSize',16, 'Color', color1);
    end
    
% ##### Specified obj_node, corresp #####
elseif nargin == 5 && strcmp(obj, 'corresp')
    scatter3(obj_node_an4(:,1), obj_node_an4(:,2), obj_node_an4(:,3), 20, 'filled', 'MarkerFaceColor',color1, 'MarkerEdgeColor',color1);
    scatter3(obj_node_an5(:,1), obj_node_an5(:,2), obj_node_an5(:,3), 20, 'filled', 'MarkerFaceColor',color2, 'MarkerEdgeColor',color2);
    for i = 1:length(obj_node_an4)
        plot3([obj_node_an4(i,1), obj_node_an5(x_to_y(i),1)], [obj_node_an4(i,2), obj_node_an5(x_to_y(i),2)], [obj_node_an4(i,3), obj_node_an5(x_to_y(i),3)], 'k', 'LineWidth', 1);
    end
% ##### Specified obj_node, distort_tri #####
elseif nargin == 5 && strcmp(obj, 'distort_tri')
    scatter3(obj_node_an4(:,1), obj_node_an4(:,2), obj_node_an4(:,3), 20, 'filled', 'MarkerFaceColor',color1, 'MarkerEdgeColor',color1);
    scatter3(obj_node_an5(:,1), obj_node_an5(:,2), obj_node_an5(:,3), 20, 'filled', 'MarkerFaceColor',color2, 'MarkerEdgeColor',color2);
    for i = 1:length(obj_node_an5)
        plot3([obj_node_an4(i,1), obj_node_an5(i,1)], [obj_node_an4(i,2), obj_node_an5(i,2)], [obj_node_an4(i,3), obj_node_an5(i,3)], 'k', 'LineWidth', 1);
    end
    num_obj_tri = length(obj_node_an5)/3;
    for i = 1:num_obj_tri
        plot3([obj_node_an4(i,1), obj_node_an4(i+num_obj_tri,1)], [obj_node_an4(i,2), obj_node_an4(i+num_obj_tri,2)], [obj_node_an4(i,3), obj_node_an4(i+num_obj_tri,3)], 'Color', color1, 'LineWidth', 2);
        plot3([obj_node_an4(i+num_obj_tri,1), obj_node_an4(i+2*num_obj_tri,1)], [obj_node_an4(i+num_obj_tri,2), obj_node_an4(i+2*num_obj_tri,2)], [obj_node_an4(i+num_obj_tri,3), obj_node_an4(i+2*num_obj_tri,3)], 'Color', color1, 'LineWidth', 2);
        plot3([obj_node_an4(i,1), obj_node_an4(i+2*num_obj_tri,1)], [obj_node_an4(i,2), obj_node_an4(i+2*num_obj_tri,2)], [obj_node_an4(i,3), obj_node_an4(i+2*num_obj_tri,3)], 'Color', color1, 'LineWidth', 2);
        plot3([obj_node_an5(i,1), obj_node_an5(i+num_obj_tri,1)], [obj_node_an5(i,2), obj_node_an5(i+num_obj_tri,2)], [obj_node_an5(i,3), obj_node_an5(i+num_obj_tri,3)], 'Color', color2, 'LineWidth', 2);
        plot3([obj_node_an5(i+num_obj_tri,1), obj_node_an5(i+2*num_obj_tri,1)], [obj_node_an5(i+num_obj_tri,2), obj_node_an5(i+2*num_obj_tri,2)], [obj_node_an5(i+num_obj_tri,3), obj_node_an5(i+2*num_obj_tri,3)], 'Color', color2, 'LineWidth', 2);
        plot3([obj_node_an5(i,1), obj_node_an5(i+2*num_obj_tri,1)], [obj_node_an5(i,2), obj_node_an5(i+2*num_obj_tri,2)], [obj_node_an5(i,3), obj_node_an5(i+2*num_obj_tri,3)], 'Color', color2, 'LineWidth', 2);
    end
end

end








