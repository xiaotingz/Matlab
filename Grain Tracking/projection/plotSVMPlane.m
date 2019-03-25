function plotSVMPlane(features, face_tri_node_1, face_tri_node_2, x_to_y, normal, bias, x_range, y_range, scale, color)
% ############################################################################
% NOTES
%   - use together with visualizeFace.m, plotTriNormals.m, calcFaceMigBySVMPlaneProj.m
%   - This script can plot one or mutiple SVM planes.
%         To plot one SVM plane, given input = [normal, bias, x_range, y_range, scale]
%         To plot several SVM planes, given input = [features]
% Inputs
%   Either
%   - features = [n, 9], [face_id, node_id, coordinates, normals, cluster_id]
%       see calcFaceMigBySVMPlaneProj.m
%   - face_tri_node_ = [m, 3], node_ids on face triangles
%   - x_to_y = [m, 1], returned by solveNodeCorresp.m
%   Or
%   - normal = [3, 1], the plane normal direction
%   - bias, scalar, the bias term
%   - x_range, [2, 1], [x_min, x_max]
%   - y_range, [2, 1], [y_min, y_max]
%   - scale, the scale for the plane normal arrow
% ############################################################################
% ------------------ load data for debug --------------------
% normal = SVMModel.Beta;
% bias = SVMModel.Bias;
% x_range = [min(X(:,1)), max(X(:,1))];
% y_range = [min(X(:,2)), max(X(:,2))];
% visualizeFace(face_node_info, x_to_y)
% -----------------------------------------------------------

% ############################## Plotting multiple SVM planes for calcFaceMigBySVMPlaneProj.m ############################## 
if nargin == 4
    colors = get(gca,'colororder');
    colors = repmat(colors, 3, 1);
    hold on
    
    % ##### Plot Mesh #####
    % """
    % same as in visualizeFace.m
    % """
    face_map_an4 = containers.Map(unique(face_tri_node_1), (1:length(unique(face_tri_node_1)))');
    face_map_an5 = containers.Map(unique(face_tri_node_2), (1:length(unique(face_tri_node_2)))');
    tri_connect_1 = zeros(size(face_tri_node_1, 1), 3);
    tri_connect_2 = zeros(size(face_tri_node_2, 1), 3);
    for i = 1:size(tri_connect_1, 1)
        for j = 1:3
            tri_connect_1(i, j) = face_map_an4(face_tri_node_1(i,j));
        end
    end
    for i = 1:size(tri_connect_2, 1)
        for j = 1:3
            tri_connect_2(i, j) = face_map_an5(face_tri_node_2(i,j));
        end
    end

    node_coords_1 = features(features(:,1)==1, 3:5);
    node_coords_2 = features(features(:,1)==2, 3:5);
    trisurf(tri_connect_1, node_coords_1(:,1), node_coords_1(:,2), node_coords_1(:,3),'Facecolor',colors(1,:), 'Facealpha', 0.3, 'edgealpha', 0.3);
    trisurf(tri_connect_2, node_coords_2(:,1), node_coords_2(:,2), node_coords_2(:,3),'Facecolor',colors(3,:), 'Facealpha', 0.3, 'edgealpha', 0.3);
    rotate3d on
    for i = 1:length(x_to_y)
        plot3([node_coords_1(i,1), node_coords_2(x_to_y(i),1)], [node_coords_1(i,2), node_coords_2(x_to_y(i),2)], [node_coords_1(i,3), node_coords_2(x_to_y(i),3)], 'k', 'LineWidth', 1);
    end

    % ##### Plot SVM plane #####
    % --- limite the plot range --- 
    x_limscale = 0.05*(max(features(:,3)) - min(features(:,3)));
    y_limscale = 0.05*(max(features(:,4)) - min(features(:,4)));
    z_limscale = 0.05*(max(features(:,5)) - min(features(:,5)));
    xlim([min(features(:,3)) - x_limscale, max(features(:,3)) + x_limscale]);
    ylim([min(features(:,4)) - y_limscale, max(features(:,4)) + y_limscale]);
    zlim([min(features(:,5)) - z_limscale, max(features(:,5)) + z_limscale]);
    
    % --- plot all clusters --- 
    for i = 1:max(features(:,end))
%       for i = 5
        color = colors(3+i, :);
        
        % --- color nodes according to its cluster --- 
        scatter3(features(features(:,end)==i, 3), features(features(:,end)==i, 4), features(features(:,end)==i, 5), ...
            80, 'filled', 'MarkerFaceColor',color, 'MarkerEdgeColor',color)
        
        % --- fit SVM plane within cluster --- 
        mask_cluster_i = (features(:, end) == i);
        X = features(mask_cluster_i, 3:5);
        Y = features(mask_cluster_i, 1);
        % """
        % paramters: 'BoxConstraint', 'Standardize', 'OutlierFraction'       
        % """        
        svm_model = fitcsvm(X, Y,'KernelFunction','linear');
        normal = svm_model.Beta;
        bias = svm_model.Bias;
        
        % --- prepare to plot SVM plane  --- 
        x_range = [min(X(:,1)), max(X(:,1))];
        y_range = [min(X(:,2)), max(X(:,2))];
        [xx,yy]=ndgrid(x_range(1) : (x_range(2) - x_range(1))/2 : x_range(2), y_range(1) : (y_range(2) - y_range(1))/2 : y_range(2));
            % --- calculate corresponding z: Ax + By + Cz + D = 0  --- 
        z = (-normal(1)*xx - normal(2)*yy - bias)/normal(3);
        scale = ceil(max(sum(features(:,1)==1), sum(features(:,1)==2))/100);
        scale = min(scale, 6);
        
        % --- plot the SVM seperation plane and plane normal--- 
        surf(xx,yy,z, 'FaceColor', color, 'EdgeColor', color, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
        idx = 2;
        quiver3(xx(idx)+0.1, yy(idx), z(idx), normal(1), normal(2), normal(3), scale, 'color', color,  'LineWidth', 3, 'MaxHeadSize', 3);
        quiver3(xx(idx)+0.1, yy(idx), z(idx), -normal(1), -normal(2), -normal(3), scale, 'color', color, 'LineWidth', 3,'MaxHeadSize', 3);
    end
    

    
    
% ############################## Plotting A Single SVM Plane ############################## 
elseif nargin == 5
    color = [0.4660, 0.6740, 0.1880];
    [xx,yy]=ndgrid(x_range(1) : (x_range(2) - x_range(1))/2 : x_range(2), y_range(1) : (y_range(2) - y_range(1))/2 : y_range(2));
    % --- calculate corresponding z: Ax + By + Cz + D = 0  --- 
    z = (-normal(1)*xx - normal(2)*yy - bias)/normal(3);
    scale = min(scale, 6);

    % --- plot the seperation plane and plane normal--- 
    surf(xx,yy,z, 'FaceColor', color, 'EdgeColor', color, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
    idx = 2;
    quiver3(xx(idx)+0.1, yy(idx), z(idx), normal(1), normal(2), normal(3), scale, 'color', color,  'LineWidth', 3, 'MaxHeadSize', 3);
    quiver3(xx(idx)+0.1, yy(idx), z(idx), -normal(1), -normal(2), -normal(3), scale, 'color', color, 'LineWidth', 3,'MaxHeadSize', 3);
    % --- plot support vector points  --- 
    % plot3(X(SVMModel.IsSupportVector,1),X(SVMModel.IsSupportVector,2),X(SVMModel.IsSupportVector,3), 'ko');
    % --- set axis ratio to be realistic ---
    
elseif nargin == 6
    [xx,yy]=ndgrid(x_range(1) : (x_range(2) - x_range(1))/2 : x_range(2), y_range(1) : (y_range(2) - y_range(1))/2 : y_range(2));
    % --- calculate corresponding z: Ax + By + Cz + D = 0  --- 
    z = (-normal(1)*xx - normal(2)*yy - bias)/normal(3);
    scale = min(scale, 6);

    % --- plot the seperation plane and plane normal--- 
    surf(xx,yy,z, 'FaceColor', color, 'EdgeColor', color, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
    idx = 2;
    quiver3(xx(idx)+0.1, yy(idx), z(idx), normal(1), normal(2), normal(3), scale, 'color', color,  'LineWidth', 3, 'MaxHeadSize', 3);
    quiver3(xx(idx)+0.1, yy(idx), z(idx), -normal(1), -normal(2), -normal(3), scale, 'color', color, 'LineWidth', 3,'MaxHeadSize', 3);
    % --- plot support vector points  --- 
    % plot3(X(SVMModel.IsSupportVector,1),X(SVMModel.IsSupportVector,2),X(SVMModel.IsSupportVector,3), 'ko');
    % --- set axis ratio to be realistic ---
end

daspect([1 1 1])
end


