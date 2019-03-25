function visualizeFace(obj_face, x_to_y, corresp, obj_node_an4, obj_node_an5, obj)
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
% x_to_y = X_to_Y{idx};
% ---------------------------------------------------------------
color1 = [0, 0.4470, 0.7410];
color2 = [0.9290, 0.6940, 0.1250];

% ##### Convert Triangle Connection to be Local #####
face_map_an4 = containers.Map(unique(obj_face{5,1}), (1:length(unique(obj_face{5,1})))');
face_map_an5 = containers.Map(unique(obj_face{5,2}), (1:length(unique(obj_face{5,2})))');
tri_connect_1 = zeros(size(obj_face{5,1}, 1), 3);
tri_connect_2 = zeros(size(obj_face{5,2}, 1), 3);
for i = 1:size(tri_connect_1, 1)
    for j = 1:3
        tri_connect_1(i, j) = face_map_an4(obj_face{5,1}(i,j));
    end
end
for i = 1:size(tri_connect_2, 1)
    for j = 1:3
        tri_connect_2(i, j) = face_map_an5(obj_face{5,2}(i,j));
    end
end


% ##### Plot Mesh #####
trisurf(tri_connect_1, obj_face{4,1}(:,1), obj_face{4,1}(:,2), obj_face{4,1}(:,3),'Facecolor',color1, 'Facealpha', 0.3, 'edgealpha', 0.3);
hold on
trisurf(tri_connect_2, obj_face{4,2}(:,1), obj_face{4,2}(:,2), obj_face{4,2}(:,3),'Facecolor',color2, 'Facealpha', 0.3, 'edgealpha', 0.3);
rotate3d on

% ##### If No obj_node Specified, Default Will Use All Nodes On Surface #####
if nargin == 2
    obj_node_an4 = obj_face{4,1};
    obj_node_an5 = obj_face{4,2};
    if any(isnan(x_to_y)) || ismember(0, x_to_y)
        obj_node_an4_use = obj_node_an4(~isnan(x_to_y), :);
        obj_node_an4_use = obj_node_an4_use((x_to_y ~= 0), :);
        node_id_an5 = setdiff(unique(x_to_y), 0);
        node_id_an5 = node_id_an5(~isnan(node_id_an5));
        obj_node_an5_use = obj_node_an5(node_id_an5, :);
    else
        obj_node_an4_use = obj_node_an4;
        obj_node_an5_use = obj_node_an5;
    end
    scatter3(obj_node_an4_use(:,1), obj_node_an4_use(:,2), obj_node_an4_use(:,3), 20, 'filled', 'MarkerFaceColor',color1, 'MarkerEdgeColor',color1);
    scatter3(obj_node_an5_use(:,1), obj_node_an5_use(:,2), obj_node_an5_use(:,3), 20, 'filled', 'MarkerFaceColor',color2, 'MarkerEdgeColor',color2);
    for i = 1:length(x_to_y)
        if ~isnan(x_to_y(i)) && x_to_y(i)~= 0 
            plot3([obj_node_an4(i,1), obj_node_an5(x_to_y(i),1)], [obj_node_an4(i,2), obj_node_an5(x_to_y(i),2)], [obj_node_an4(i,3), obj_node_an5(x_to_y(i),3)], 'k', 'LineWidth', 1);
%         text(obj_point_smooth_an4(i,1), obj_point_smooth_an4(i,2), obj_point_smooth_an4(i,3),[' ', num2str(int32(i))],'FontSize',16, 'Color', color1);
        end
    end

% ##### In the case of node corresp by min_dist, only some nodes will have corresps #####
elseif nargin == 3
    obj_node_an4 = obj_face{4,1};
    obj_node_an5 = obj_face{4,2};
    use_an4 = corresp(:,1);
    use_an5 = corresp(:,2);
    obj_node_an4_use = obj_node_an4(use_an4, :);
    obj_node_an5_use = obj_node_an5(use_an5, :);
    scatter3(obj_node_an4_use(:,1), obj_node_an4_use(:,2), obj_node_an4_use(:,3), 20, 'filled', 'MarkerFaceColor',color1, 'MarkerEdgeColor',color1);
    scatter3(obj_node_an5_use(:,1), obj_node_an5_use(:,2), obj_node_an5_use(:,3), 20, 'filled', 'MarkerFaceColor',color2, 'MarkerEdgeColor',color2);
    
    for i = 1:length(x_to_y)
        if ~isnan(x_to_y(i))
            plot3([obj_node_an4(i,1), obj_node_an5(x_to_y(i),1)], [obj_node_an4(i,2), obj_node_an5(x_to_y(i),2)], [obj_node_an4(i,3), obj_node_an5(x_to_y(i),3)], 'k', 'LineWidth', 1);
%         text(obj_point_smooth_an4(i,1), obj_point_smooth_an4(i,2), obj_point_smooth_an4(i,3),[' ', num2str(int32(i))],'FontSize',16, 'Color', color1);
        end
    end
    
% ##### Specified obj_node, corresp #####
elseif nargin == 5 && strcmp(obj, 'corresp')
    scatter3(obj_node_an4(:,1), obj_node_an4(:,2), obj_node_an4(:,3), 20, 'filled', 'MarkerFaceColor',color1, 'MarkerEdgeColor',color1);
    scatter3(obj_node_an5(:,1), obj_node_an5(:,2), obj_node_an5(:,3), 20, 'filled', 'MarkerFaceColor',color2, 'MarkerEdgeColor',color2);
    for i = 1:length(obj_node_an4)
        if ~isnan(x_to_y(i))
            plot3([obj_node_an4(i,1), obj_node_an5(x_to_y(i),1)], [obj_node_an4(i,2), obj_node_an5(x_to_y(i),2)], [obj_node_an4(i,3), obj_node_an5(x_to_y(i),3)], 'k', 'LineWidth', 1);
        end
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


daspect([1 1 1])

end








