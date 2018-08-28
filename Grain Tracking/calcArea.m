function area = calcArea(node_coords)
% ###########################################################################
% * Input 
%     - node_coords = [3, 3], each row contain [x, y, z] of a node
% ###########################################################################
% --- The length of the cross product of two vectors is equal to the area of the parallelogram determined by the two vectors ---
    area = 1/2*norm(cross(node_coords(1,:) - node_coords(2,:), node_coords(3,:) - node_coords(1,:)));
end

