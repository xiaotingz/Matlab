% """
% - In simpleware, every mesh triangle have only one label. 
% - Node the mesh triangles seem to be all unique, but there are duplicates
% in mesh nodes. 
% """

grain_id1 = 1;
mask_1 = Model1_surface_triangles_MatlabPartIndex == grain_id1;
tris_1 = Model1_surface_triangles(mask_1, :);
grain_id2 = 6;
mask_2 = Model1_surface_triangles_MatlabPartIndex == grain_id2;
tris_2 = Model1_surface_triangles(mask_2, :);
coords = Model1_surface_vertices;

trimesh(tris_1, coords(:,1), coords(:,2), coords(:,3), 'edgecolor', [0, 0.4470, 0.7410]);
hold on
trimesh(tris_2, coords(:,1), coords(:,2), coords(:,3), 'edgecolor', [0.8500    0.3250    0.0980]);
daspect([1,1,1])

rotate3d on
hold off

figure
trimesh(tris_1, coords(:,1), coords(:,2), coords(:,3), 'edgecolor', [0, 0.4470, 0.7410]);
rotate3d on
daspect([1,1,1])

figure
trimesh(tris_2, coords(:,1), coords(:,2), coords(:,3), 'edgecolor', [0.8500    0.3250    0.0980]);
rotate3d on
daspect([1,1,1])









