%%  ##### MATLAB Simple Examples #####
mu = [1 -1];

SIGMA = [.9 .4; .4 .3];
[x_corner,y_corner] = meshgrid(linspace(-1,3,25)',linspace(-3,1,25)');
z_corner = 4 * mvncdf([x_corner(:) y_corner(:)],mu,SIGMA);
tri_corner = delaunay(x_corner(:), y_corner(:));

mesh_corner = triangulation(tri_corner, [x_corner(:) y_corner(:), z_corner(:)]);
tri_center_corner = incenter(mesh_corner);
tri_normal_corner = faceNormal(mesh_corner);

% ----- visualize mesh -----
% trimesh(mesh_corner)
% rotate3d on
% hold on  
% quiver3(tri_center_corner(:,1),tri_center_corner(:,2),tri_center_corner(:,3), ...
%      tri_normal_corner(:,1),tri_normal_corner(:,2),tri_normal_corner(:,3),0.5,'color','r');


% ----- solve correpondence -----
X = [x_corner(:), y_corner(:), z_corner(:)];
Y = [x_corner(:), y_corner(:), zeros(length(x_corner(:)), 1)];
[x_to_y, y_to_x] = solveNodeCorresp(X, Y);


% SIGMA = [.6 .2; .2 .3];
% [x_ribbon, y_ribbon] = meshgrid(linspace(-3,3,30)',linspace(0,2,10)');
% z_ribbon = mvncdf([x_ribbon(:) y_ribbon(:)],mu,SIGMA);
% tri_ribbon = delaunay(x_ribbon(:), y_ribbon(:));
% trisurf(tri_ribbon, x_ribbon(:), y_ribbon(:), z_ribbon(:));
% rotate3d on



%%  ##### Surface from D3D #####

% ----- visualizeSingleFaceWithNormal(file, obj_facelabel, reverse_winding) -----
% visualizeSingleFaceWithNormal(file_an4, obj_facelabel, 1)



