function h = calcHeightTruncatedConeModel(coords_1, coords_2)
% ############################################################################
% * Inputs
%     - coords_ = [3,3]
%           Coordinates of the three nodes on the top and bottom surface.
% * Notes
%     - Given two surface areas and the volume, calculate heigt as if the
%       geometry is truncated cone
%     - Refs, pillar
%           http://mathworld.wolfram.com/PyramidalFrustum.html 
%           https://math.stackexchange.com/questions/1966507/the-volume-for-truncated-pyramid-with-irregular-base
%     - Refs, triangle area. There are two ways to calculate triangle
%     areas, the difference seems to be just numerical error. 
%           https://www.mathworks.com/matlabcentral/answers/14928-area-of-triangle
%           http://demonstrations.wolfram.com/Signed2DTriangleAreaFromTheCrossProductOfEdgeVectors/
% ############################################################################
% ----------------------- load debug data -----------------------
% area1 = 16;
% area2 = 23.309;
% v = 25.9274;
% coords_1 = [0,0,0;1,0,0;0,1,0];
% coords_2 = [0,0,1;1,0,1;0,1,1];
% ---------------------------------------------------------------
%     ons = [1 1 1];

%     x_1 = coords_1(:,1)';
%     y_1 = coords_1(:,2)';
%     z_1 = coords_1(:,3)';
%     area_1 = 0.5 * sqrt(det([x_1;y_1;ons])^2 + det([y_1;z_1;ons])^2 + det([z_1;x_1;ons])^2);
    area_1 = 1/2 * norm(cross(coords_1(1,:) - coords_1(2,:), coords_1(3,:) - coords_1(1,:)));

%     x_2 = coords_2(:,1)';
%     y_2 = coords_2(:,2)';
%     z_2 = coords_2(:,3)';
%     area_2 = 0.5 * sqrt(det([x_2;y_2;ons])^2 + det([y_2;z_2;ons])^2 + det([z_2;x_2;ons])^2);
    area_2 = 1/2 * norm(cross(coords_2(1,:) - coords_2(2,:), coords_2(3,:) - coords_2(1,:)));

    coords = [coords_1; coords_2];
    [~,v] = convhull(coords(:,1), coords(:,2), coords(:,3));
    
    h = 3*v/(area_1 + area_2 + sqrt(area_1*area_2));
end


