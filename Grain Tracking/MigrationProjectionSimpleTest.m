%%  ##### MATLAB Simple Examples #####
mu = [1 -1];

% ##### A Corner-Shape Surface #####
% SIGMA = [.9 .4; .4 .3];
% [x,y] = meshgrid(linspace(-1,3,25)',linspace(-3,1,25)');
% x = x(:);
% y = y(:);
% z = 4 * mvncdf([x y],mu,SIGMA);
% tri = delaunay(x, y);
% ##### A Ribbon-Shape Surface #####
% SIGMA = [.6 .2; .2 .3];
% [x, y] = meshgrid(linspace(-3,3,30)',linspace(0,2,10)');
% z = mvncdf([x y],mu,SIGMA);
% tri = delaunay(x, y);
% ##### circle X Y #####
SIGMA = [.9 .4; .4 .3];
[T,R] = meshgrid(linspace(0,2*pi,64),linspace(0,2,16));
x = R.*cos(T);
y = R.*sin(T);
x = x(:);
y = y(:);
z = 20 * mvncdf([x y],mu,SIGMA);
tri = delaunay(x, y);

tri_mesh = triangulation(tri, [x y, z(:)]);
tri_center = incenter(tri_mesh);
tri_normal = faceNormal(tri_mesh);

% ----- the migration of a triangle is the average of its three nodes -----
tri_migration = sum(z(tri),2)/3 * [0,0,1];
tri_migration_project = dot(tri_migration, tri_normal, 2);

avg_migration = sum(tri_migration(:,3))/length(tri_migration);
% ----- projection = dot(a,b)/|b|, |b|=1 in this case -----
avg_migration_project = sum(tri_migration_project)/length(tri_migration);
avg_migration_diff = avg_migration - avg_migration_project;


% ##### Alpha Shape #####
% x_full = [x; x];
% y_full = [y; y];
% z_full = [z; zeros(length(x), 1)];
% shp = alphaShape(x_full, y_full, z_full);
% shp.Alpha = 0.5;
% plot(shp)
% rotate3d on
% volume = shp.volume;

% ##### Migration by individual volume/surface_area #####
node_coord_an4 = [x, y, z]; 
node_coord_an5 = [x, y, zeros(length(x),1)];

avg_migration_trunctedcone_sep = zeros(size(tri, 1), 1);
area1 = 0;
area2 = 0;
v = 0;
for i = 1:length(avg_migration_trunctedcone_sep)
    tri1_nodes = node_coord_an4(tri(i, :), :);
    tri2_nodes = node_coord_an5(tri(i, :), :);
    DT = delaunayTriangulation([tri1_nodes; tri2_nodes]);
    [C, v_tmp] = convexHull(DT);
    area1_tmp = calcArea(tri1_nodes);
    area2_tmp = calcArea(tri2_nodes);
    avg_migration_trunctedcone_sep(i) = heightUsingTruncatedConeModel(area1_tmp, area2_tmp, v_tmp);
    v= v + v_tmp;
    area1 = area1 + area1_tmp;
    area2 = area2 + area2_tmp;
end
migration = sum(avg_migration_trunctedcone_sep)/length(avg_migration_trunctedcone_sep);


% ##### Migration by entire volume/surface_area #####
avg_migration_trunctedcone_itg = heightUsingTruncatedConeModel(area1, area2, v);

% ----- the volume by convex hull is actually larger than the real volume -----
% DT = delaunayTriangulation([x_full, y_full, z_full]);
% [C, v] = convexHull(DT);
% trisurf(C,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3), 'Facecolor', [0.4660, 0.6740, 0.1880], 'edgecolor', 'k');



% ##### Visualize Mesh #####
color1 = [0, 0.4470, 0.7410];
color2 = [0.9290, 0.6940, 0.1250];
trisurf(tri, x, y, z(:), 'Facealpha', 0, 'edgecolor', color1)
rotate3d on
hold on  
quiver3(tri_center(:,1),tri_center(:,2),tri_center(:,3), ...
     tri_normal(:,1),tri_normal(:,2),tri_normal(:,3),0.5,'color','r', ...
     'AutoScaleFactor', 3);
trisurf(tri, x, y, zeros(length(x), 1), 'Facealpha', 0, 'edgecolor', color2);
% plot3([x(end), x(end)], [y(end), y(end)], [0, z(end)], 'k', 'LineWidth', 1);
% idx = 573;
% plot3([x(idx), x(idx)], [y(idx), y(idx)], [0, z(idx)], 'k', 'LineWidth', 1);


% ----- solve correpondence -----
% X = [x, y, z(:)];
% Y = [x, y, zeros(length(x), 1)];
% [x_to_y, y_to_x] = solveNodeCorresp(X, Y);




%% ##### volume of three simple prisms ##### 
figure()
points_1 = [0,0,0; 2,0,0; 0,1,0; 0,0,1; 4,0,1; 0,2,1];  
DT_1 = delaunayTriangulation(points_1);
[C_1, v_1] = convexHull(DT_1);
trisurf(C_1,DT_1.Points(:,1),DT_1.Points(:,2),DT_1.Points(:,3));
rotate3d on

figure()
points_2 = [0,0,0; 1,0,0; 0,1,0; 1,1,0; 0,0,1; 2,0,1; 0,2,1; 2,2,1];  
DT_2 = delaunayTriangulation(points_2);
[C_2, v_2] = convexHull(DT_2);
trisurf(C_2,DT_2.Points(:,1),DT_2.Points(:,2),DT_2.Points(:,3));
rotate3d on

figure()
points_3 = [0,0,0; 1,0,0; 0,1,0; 1,1,0; 0,0,1; 4,0,1; 0,1,1; 4,1,1];  
DT_3 = delaunayTriangulation(points_3);
[C_3, v_3] = convexHull(DT_3);
trisurf(C_3,DT_3.Points(:,1),DT_3.Points(:,2),DT_3.Points(:,3));
rotate3d on


