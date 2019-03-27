function min_angle = calcTriangleMinAngle(tri_coords)
% ############################################################################
% Input
%     - tri_coords = [3, 3], each row is [x, y, z] of one node, each
%     column is a node.
% Output
%     - min_angle is the smallest inner angle of the triangle
% Notes
%     - Can be used to check triangle shape
% ############################################################################
% ----------------------- load debug data -----------------------
% tri_coords = [0,0,0; 1,0,0; 0,sqrt(3),0];
% ---------------------------------------------------------------
edge_21 = tri_coords(2, :) - tri_coords(1, :);
edge_31 = tri_coords(3, :) - tri_coords(1, :);
edge_12 = tri_coords(1, :) - tri_coords(2, :);
edge_32 = tri_coords(3, :) - tri_coords(2, :);
edge_13 = tri_coords(1, :) - tri_coords(3, :);
edge_23 = tri_coords(2, :) - tri_coords(3, :);

angle_1 = atan2d(norm(cross(edge_21,edge_31)),dot(edge_21,edge_31));
angle_2 = atan2d(norm(cross(edge_12,edge_32)),dot(edge_12,edge_32));
angle_3 = atan2d(norm(cross(edge_13,edge_23)),dot(edge_13,edge_23));

min_angle = min([angle_1, angle_2, angle_3]);
end
