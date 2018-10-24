function plotSVMPlane(features, normal, bias, x_range, y_range, scale, color)
% ############################################################################
% NOTES
%   - use together with visualizeFace.m, plotTriNormals.m, calcFaceMigBySVMPlaneProj.m
%   - This script can plot one or mutiple SVM planes.
%         To plot one SVM plane, given input = [normal, bias, x_range, y_range, scale]
%         To plot several SVM planes, given input = [features]
% Inputs
%   Either
%   - features = [9, n], [face_id, node_id, coordinates, normals, cluster_id]
%       see calcFaceMigBySVMPlaneProj.m
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
if nargin == 1
    colors = get(gca,'colororder');
    hold on
    for i = 1:max(features(:,end))
        if 3+i <= length(colors)
            color = colors(3+i, :);
        else
            color = colors(3+i-length(colors), :);
        end
        scatter3(features(features(:,end)==i, 3), features(features(:,end)==i, 4), features(features(:,end)==i, 5), ...
            80, 'filled', 'MarkerFaceColor',color, 'MarkerEdgeColor',color)
        mask_cluster_i = (features(:, end) == i);
        X = features(mask_cluster_i, 3:5);
        Y = features(mask_cluster_i, 1);
        SVMModel = fitcsvm(X, Y,'KernelFunction','linear', 'Standardize',false,'ClassNames',{'1','2'});
        normal = SVMModel.Beta;
        bias = SVMModel.Bias;
        x_range = [min(X(:,1)), max(X(:,1))];
        y_range = [min(X(:,2)), max(X(:,2))];
        [xx,yy]=ndgrid(x_range(1) : (x_range(2) - x_range(1))/2 : x_range(2), y_range(1) : (y_range(2) - y_range(1))/2 : y_range(2));
        % --- calculate corresponding z: Ax + By + Cz + D = 0  --- 
        z = (-normal(1)*xx - normal(2)*yy - bias)/normal(3);
        scale = ceil(max(sum(features(:,1)==1), sum(features(:,1)==2))/100);
        scale = min(scale, 6);
        % --- plot the seperation plane and plane normal--- 
        surf(xx,yy,z, 'FaceColor', color, 'EdgeColor', color, 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
        idx = 2;
        quiver3(xx(idx)+0.1, yy(idx), z(idx), normal(1), normal(2), normal(3), scale, 'color', color,  'LineWidth', 3, 'MaxHeadSize', 3);
        quiver3(xx(idx)+0.1, yy(idx), z(idx), -normal(1), -normal(2), -normal(3), scale, 'color', color, 'LineWidth', 3,'MaxHeadSize', 3);
    end
    % --- limite the plot range --- 
    x_limscale = 0.05*(max(features(:,3)) - min(features(:,3)));
    y_limscale = 0.05*(max(features(:,4)) - min(features(:,4)));
    z_limscale = 0.05*(max(features(:,5)) - min(features(:,5)));
    xlim([min(features(:,3)) - x_limscale, max(features(:,3)) + x_limscale]);
    ylim([min(features(:,4)) - y_limscale, max(features(:,4)) + y_limscale]);
    zlim([min(features(:,5)) - z_limscale, max(features(:,5)) + z_limscale]);

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


