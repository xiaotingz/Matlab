function plotSVMPlane(normal, bias, x_range, y_range, scale)
% ############################################################################
% Inputs
%   - normal = [3, 1], the plane normal direction
%   - bias, scalar, the bias term
%   - x_range, [2, 1], [x_min, x_max]
%   - y_range, [2, 1], [y_min, y_max]
%   - scale, the scale for the plane normal arrow
% Notes
%   - use together with visualizeFace.m and visualizeSingleFaceWithNormal.m
% ############################################################################
% ------------------ load data for debug --------------------
% normal = SVMModel.Beta;
% bias = SVMModel.Bias;
% x_range = [min(X(:,1)), max(X(:,1))];
% y_range = [min(X(:,2)), max(X(:,2))];
% visualizeFace(face_node_info, x_to_y)
% -----------------------------------------------------------
[xx,yy]=ndgrid(x_range(1) : (x_range(2) - x_range(1))/2 : x_range(2), y_range(1) : (y_range(2) - y_range(1))/2 : y_range(2));
% --- calculate corresponding z: Ax + By + Cz + D = 0  --- 
z = (-normal(1)*xx - normal(2)*yy - bias)/normal(3);
scale = min(scale, 6);

% --- plot the seperation plane and plane normal--- 
surf(xx,yy,z, 'FaceColor', [0.4660, 0.6740, 0.1880], 'EdgeColor', [0.4660, 0.6740, 0.1880], 'FaceAlpha', 0.3, 'EdgeAlpha', 0.3);
idx = 2;
quiver3(xx(idx)+0.1, yy(idx), z(idx), normal(1), normal(2), normal(3), scale, 'color', [0.4660, 0.6740, 0.1880],  'LineWidth', 3, 'MaxHeadSize', 3);
quiver3(xx(idx)+0.1, yy(idx), z(idx), -normal(1), -normal(2), -normal(3), scale, 'color', [0.4660, 0.6740, 0.1880], 'LineWidth', 3,'MaxHeadSize', 3);
% --- plot support vector points  --- 
% plot3(X(SVMModel.IsSupportVector,1),X(SVMModel.IsSupportVector,2),X(SVMModel.IsSupportVector,3), 'ko');
% --- set axis ratio to be realistic ---
daspect([1 1 1])

end


