% ############################################################################
% Simple Examples 1). Direct, Projection and Height=Volume/Area
% ############################################################################
mu = [1 -1];

% ##### A Corner-Shape Surface #####
% SIGMA = [.9 .4; .4 .3];
% [x,y] = meshgrid(linspace(-1,3,25)',linspace(-3,1,25)');
% x = x(:);
% y = y(:);
% z = 50 * mvncdf([x y],mu,SIGMA);
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

migration_direct = sum(tri_migration(:,3))/length(tri_migration);
% ----- projection = dot(a,b)/|b|, |b|=1 in this case -----
migration_project = sum(tri_migration_project)/length(tri_migration);

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
node_coord_an5 = [x, y, -z];

migration_trunctedcone_seps = zeros(size(tri, 1), 1);
area1 = 0;
area2 = 0;
v = 0;
for i = 1:length(migration_trunctedcone_seps)
    tri_top = node_coord_an4(tri(i, :), :);
    tri_bottom = node_coord_an5(tri(i, :), :);
    DT = delaunayTriangulation([tri_top; tri_bottom]);
    [C, v_tmp] = convexHull(DT);
    area1_tmp = calcArea(tri_top);
    area2_tmp = calcArea(tri_bottom);
    migration_trunctedcone_seps(i) = heightUsingTruncatedConeModel(area1_tmp, area2_tmp, v_tmp);
    v= v + v_tmp;
    area1 = area1 + area1_tmp;
    area2 = area2 + area2_tmp;
end
migration_trunctedcone_sep = sum(migration_trunctedcone_seps)/length(migration_trunctedcone_seps);


% ##### Migration by entire volume/surface_area #####
migration_trunctedcone_itg = heightUsingTruncatedConeModel(area1, area2, v);

% ----- the volume by convex hull is actually larger than the real volume -----
% DT = delaunayTriangulation([x_full, y_full, z_full]);
% [C, v] = convexHull(DT);
% trisurf(C,DT.Points(:,1),DT.Points(:,2),DT.Points(:,3), 'Facecolor', [0.4660, 0.6740, 0.1880], 'edgecolor', 'k');



% ##### Visualize Mesh #####
color1 = [0, 0.4470, 0.7410];
color2 = [0.9290, 0.6940, 0.1250];
trisurf(tri, node_coord_an4(:,1), node_coord_an4(:,2), node_coord_an4(:,3), 'Facealpha', 0, 'edgecolor', color1)
rotate3d on
hold on  
% quiver3(tri_center(:,1),tri_center(:,2),tri_center(:,3), ...
%      tri_normal(:,1),tri_normal(:,2),tri_normal(:,3),0.5,'color','r', ...
%      'AutoScaleFactor', 3);
trisurf(tri, node_coord_an5(:,1), node_coord_an5(:,2), node_coord_an5(:,3), 'Facealpha', 0, 'edgecolor', color2);
% plot3([x(end), x(end)], [y(end), y(end)], [0, z(end)], 'k', 'LineWidth', 1);
% idx = 573;
% plot3([x(idx), x(idx)], [y(idx), y(idx)], [0, z(idx)], 'k', 'LineWidth', 1);


%% ############################################################################
% Simple Examples 2). Shrinkage and Expansion
% ############################################################################
color1 = [0, 0.4470, 0.7410];
color2 = [0.9290, 0.6940, 0.1250];

mu = [0 0];
Sigma = [.25 .3; .3 1];

xrange = 1.6;
yrange = 3;
interval = 0.2;


% ##### Generate Surface Coordinates #####
[x_1,y_1] = meshgrid(-xrange:interval:xrange, -yrange:interval:yrange);
x_1 = x_1(:);
y_1 = y_1(:);
z_1_bell = 3 * mvnpdf([x_1, y_1],mu,Sigma);
z_1_ref = - 0.5*ones(size(x_1));
coords_1_bell = [x_1, y_1, z_1_bell];
coords_1_ref = [x_1, y_1, z_1_ref];

[x_2,y_2] = meshgrid(-xrange/2:interval:xrange/2, -yrange/2:interval:yrange/2);
x_2 = x_2(:);
y_2 = y_2(:);
z_2_ref = - 0.5 * ones(size(x_2));
coords_2_ref = [x_2, y_2, z_2_ref];


% ##### Triangulation #####
tri_1 = delaunay(x_1, y_1);
bell_tri_1 = triangulation(tri_1, coords_1_bell);
ref_tri_1 = triangulation(tri_1, coords_1_ref);
tri_2 = delaunay(x_2, y_2);
ref_tri_2 = triangulation(tri_2, coords_2_ref);

% ##### Calculate Migration as Height #####
h1_list = zeros(length(tri_1), 1);
for i = 1:length(tri_1)
    tri_top = coords_1_bell(tri_1(i,:), :);
    tri_bottom = coords_1_ref(tri_1(i,:), :);
    DT = delaunayTriangulation([tri_top; tri_bottom]);
    [C, v] = convexHull(DT);
    area_top = calcArea(tri_top);
    area_bottom = calcArea(tri_bottom);
    h1_list(i) = heightUsingTruncatedConeModel(area_top, area_bottom, v);
end
migration_trunctedcone_1 = sum(h1_list)/length(h1_list)

tri_1from2 = y_to_x(tri_2);
h2_list = zeros(length(tri_2), 1);
for i = 1:length(tri_2)
    tri_top = coords_1_bell(tri_1from2(i,:), :);
    tri_bottom = coords_2_ref(tri_2(i,:), :);
    DT = delaunayTriangulation([tri_top; tri_bottom]);
    [C, v] = convexHull(DT);
    area_top = calcArea(tri_top);
    area_bottom = calcArea(tri_bottom);
    h2_list(i) = heightUsingTruncatedConeModel(area_top, area_bottom, v);
end
migration_trunctedcone_2 = sum(h2_list)/length(h2_list)

% ##### Direct Migration #####
migration_direct_1 = sum(z_1_bell + 0.5)/length(z_1_bell)
node_migration_2 = vecnorm(coords_1_bell(tri_1from2, :) - coords_2_ref(tri_2, :), 2, 2);
migration_direct_2 = sum(node_migration_2)/length(node_migration_2)



%% ##### Visualization #####
trisurf(tri_1, coords_1_bell(:, 1), coords_1_bell(:, 2), coords_1_bell(:, 3), 'Facecolor', color1, 'Facealpha', 0.3, 'Edgealpha', 0.3)
rotate3d on
hold on
trisurf(tri_1, coords_1_ref(:, 1), coords_1_ref(:, 2), coords_1_ref(:, 3), 'Facecolor', color2, 'Facealpha', 0.3, 'Edgealpha', 0.3)
tri_1_bell_coord = coords_1_bell(tri_1, :);
tri_1_ref_coord = coords_1_ref(tri_1, :);
scatter3(tri_1_bell_coord(:,1), tri_1_bell_coord(:,2), tri_1_bell_coord(:,3), 20, 'filled', 'MarkerFaceColor',color1, 'MarkerEdgeColor',color1);
scatter3(tri_1_ref_coord(:,1), tri_1_ref_coord(:,2), tri_1_ref_coord(:,3), 20, 'filled', 'MarkerFaceColor',color2, 'MarkerEdgeColor',color2);
for i = 1:length(tri_1_ref_coord)
    plot3([tri_1_bell_coord(i,1), tri_1_ref_coord(i,1)], [tri_1_bell_coord(i,2), tri_1_ref_coord(i,2)], [tri_1_bell_coord(i,3), tri_1_ref_coord(i,3)], 'k', 'LineWidth', 1);
end
print('simpleGeo1', '-dpng', '-r300')

trisurf(tri_1, coords_1_bell(:, 1), coords_1_bell(:, 2), coords_1_bell(:, 3), 'Facecolor', color1, 'Facealpha', 0.3, 'Edgealpha', 0.3)
rotate3d on
hold on
trisurf(tri_2, coords_2_ref(:, 1), coords_2_ref(:, 2), coords_2_ref(:, 3), 'Facecolor', color2, 'Facealpha', 0.3, 'Edgealpha', 0.3)
tri_1from2_coord = coords_1_bell(tri_1from2, :);
tri_2_coord = coords_2_ref(tri_2, :);
scatter3(tri_1from2_coord(:,1), tri_1from2_coord(:,2), tri_1from2_coord(:,3), 20, 'filled', 'MarkerFaceColor',color1, 'MarkerEdgeColor',color1);
scatter3(tri_2_coord(:,1), tri_2_coord(:,2), tri_2_coord(:,3), 20, 'filled', 'MarkerFaceColor',color2, 'MarkerEdgeColor',color2);
for i = 1:length(tri_2_coord)
    plot3([tri_1from2_coord(i,1), tri_2_coord(i,1)], [tri_1from2_coord(i,2), tri_2_coord(i,2)], [tri_1from2_coord(i,3), tri_2_coord(i,3)], 'k', 'LineWidth', 1);
end
print('simpleGeo2', '-dpng', '-r300')





%% ############################################################################
% Simple Example 3). Prism with Top/Bottom Surface Shape
% ############################################################################
figure()
points_1 = [0,0,0; 2,0,0; 0,1,0; 0,0,1; 4,0,1; 0,2,1];  
DT_4 = delaunayTriangulation(points_1);
[C_4, v_4] = convexHull(DT_4);
trisurf(C_4,DT_4.Points(:,1),DT_4.Points(:,2),DT_4.Points(:,3), 'Facealpha', 0.3);
h1 = heightUsingTruncatedConeModel(1, 4, v_4)
rotate3d on
% print('pillar1', '-dpng', '-r300')

figure()
points_2 = [0,0,0; 1,0,0; 0,1,0; 1,1,0; 0,0,1; 2,0,1; 0,2,1; 2,2,1];  
DT_2 = delaunayTriangulation(points_2);
[C_2, v_2] = convexHull(DT_2);
trisurf(C_2,DT_2.Points(:,1),DT_2.Points(:,2),DT_2.Points(:,3), 'Facealpha', 0.3);
h2 = heightUsingTruncatedConeModel(1, 4, v_2)
rotate3d on

figure()
points_3 = [0,0,0; 1,0,0; 0,1,0; 1,1,0; 0,0,1; 4,0,1; 0,1,1; 4,1,1];  
DT_3 = delaunayTriangulation(points_3);
[C_3, v_3] = convexHull(DT_3);
trisurf(C_3,DT_3.Points(:,1),DT_3.Points(:,2),DT_3.Points(:,3), 'Facealpha', 0.3);
h2 = heightUsingTruncatedConeModel(1, 4, v_3)
rotate3d on

figure()
points_1 = [-1,0,0; 1,0,0; 0,1,0; 0,0,1; 4,0,1; 0,2,1];  
DT_4 = delaunayTriangulation(points_1);
[C_4, v_4] = convexHull(DT_4);
trisurf(C_4,DT_4.Points(:,1),DT_4.Points(:,2),DT_4.Points(:,3), 'Facealpha', 0.3);
daspect([2, 2, 1])
h4 = heightUsingTruncatedConeModel(1, 4, v_4)
rotate3d on

figure()
points_1 =  [-4,0,0; 4,0,0; 0,0.25,0; 0,0,1; 4,0,1; 0,2,1]; 
DT_5 = delaunayTriangulation(points_1);
[C_5, v_5] = convexHull(DT_5);
trisurf(C_5,DT_5.Points(:,1),DT_5.Points(:,2),DT_5.Points(:,3), 'Facealpha', 0.3);
daspect([2, 2, 1])
h5 = heightUsingTruncatedConeModel(1, 4, v_5)
rotate3d on

figure()
points_1 = [-10,0,0; 10,0,0; 0,0.1,0; 0,0,1; 4,0,1; 0,2,1];  
DT_6 = delaunayTriangulation(points_1);
[C_6, v_6] = convexHull(DT_6);
trisurf(C_6,DT_6.Points(:,1),DT_6.Points(:,2),DT_6.Points(:,3), 'Facealpha', 0.3);
daspect([2, 2, 1])
h6 = heightUsingTruncatedConeModel(1, 4, v_6)
rotate3d on