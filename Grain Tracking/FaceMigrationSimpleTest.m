mu = [1 -1];

SIGMA = [.9 .4; .4 .3];
[x_corner,y_corner] = meshgrid(linspace(-1,3,25)',linspace(-3,1,25)');
z_corner = mvncdf([x_corner(:) y_corner(:)],mu,SIGMA);
tri_corner = delaunay(x_corner(:), y_corner(:));
trisurf(tri_corner, x_corner(:), y_corner(:), z_corner(:));
rotate3d on
F = faceNormal(tri_corner);


% SIGMA = [.6 .2; .2 .3];
% [x_ribbon, y_ribbon] = meshgrid(linspace(-3,3,30)',linspace(0,2,10)');
% z_ribbon = mvncdf([x_ribbon(:) y_ribbon(:)],mu,SIGMA);
% tri_ribbon = delaunay(x_ribbon(:), y_ribbon(:));
% trisurf(tri_ribbon, x_ribbon(:), y_ribbon(:), z_ribbon(:));
% rotate3d on


