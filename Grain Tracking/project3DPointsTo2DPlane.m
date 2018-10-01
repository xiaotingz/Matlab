function [points_proj_2d, native2d_x] = project3DPointsTo2DPlane(points, origin, normal, cond)
% ############################################################################
% * Input 
%     - points = [n, 3], coordinates of the nodes to be projected onto the basis plane
%     - origin = [3, 1], origin of the basis plane
%     - normal = [3, 1], normal of the basis plane
%     - cond, 'bp' || a [3, 1] vector, specify local coordinate frame
%          If cond == 'bp', the points are on the basis plane, the native2d_x should
%             be calced. The choice native2d_x is not unique, can be chose randomly as long as it's 
%             consistent for all points projected onto this basis plane. 
%          If the native_x is specified, then use it.
% * Output
%     - points_proj = [n, 2], coordinates of the projected points
%     - native2d_x = [3, 1], x-axis of the local coordinate frame being used.
% * Notes
%     - use in subgraphCorrespAnalysis.m
%     - ref: https://stackoverflow.com/questions/9605556/how-to-project-a-point-onto-a-plane-in-3d
% ############################################################################
% ------------------ load data for debug --------------------
% points = nodes_an4;
% origin = centroid_bp;
% normal = normal_bp;
% cond = 'bp';
% -----------------------------------------------------------
% ##### Project the Points Onto a 3D Plane #####
v = points - repmat(origin, size(points,1), 1);
dists = v*normal';

points_proj_3d = points - dists*normal;


% ##### Project Points On the 3D Plane Onto a 2D Plane #####
if strcmp(cond, 'bp')
% --- This choice of native_x is rather random ---
    native2d_x = points_proj_3d(1,:) - origin;
    native2d_x = native2d_x/norm(native2d_x);
elseif isa(class(native2d_x), 'double') && length(native2d_x) == 3
    native2d_x = native2d_x/norm(native2d_x);
else
    warning('There is a problem with the cond parameter of projectPointsTo2DPlane.');
    warning('cond has to be a string bp or a [3,1] vector');
end
% --- native2d_y can be derived from native2d_x ---
native2d_y = cross(normal, native2d_x);
native2d_y = native2d_y/norm(native2d_y);

points_proj_2d(:,1) = (points_proj_3d - origin) * native2d_x';
points_proj_2d(:,2) = (points_proj_3d - origin) * native2d_y';


% % ##### Checks #####
% % ----- Visual check -----
% scatter3(points_proj_3d(:,1), points_proj_3d(:,2), points_proj_3d(:,3), 'filled');
% hold on
% rotate3d on
% % scatter3(points(:,1), points(:,2), points(:,3), 'filled');
% quiver3(origin(1,1), origin(2), origin(3), native2d_x(1), native2d_x(2), native2d_x(3), 5, 'color', 'r',  'LineWidth', 3, 'MaxHeadSize', 3);
% quiver3(origin(1,1), origin(2), origin(3), native2d_y(1), native2d_y(2), native2d_y(3), 5, 'color', 'g',  'LineWidth', 3, 'MaxHeadSize', 3);
% quiver3(origin(1,1), origin(2), origin(3), normal(1), normal(2), normal(3), 5, 'color', 'k',  'LineWidth', 3, 'MaxHeadSize', 3);
% daspect([1 1 1])
% 
% figure
% scatter(points_proj_2d(:,1), points_proj_2d(:,2), 'filled')
% daspect([1 1 1])
% 
% % ----- Distance invariant check -----
% norm(points_proj_2d(10, :) - points_proj_2d(15, :))
% norm(points_proj_3d(10, :) - points_proj_3d(15, :))

end
