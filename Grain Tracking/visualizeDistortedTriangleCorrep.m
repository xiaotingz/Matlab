i = 1;
obj_facelabel_an4 = tracked_uniqueface_an4(i, :);
obj_facelabel_an5 = tracked_uniqueface_an5(i, :);
x_to_y = X_to_Y{i};

% ##### get objective triangles on the objective face #####
mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));

% ##### get id and coord of the objective triangles #####
face_node_id_an4 = tri_node_an4(mask_objface_an4, :);
face_node_uniqueid_an4 = unique(face_node_id_an4);
face_node_coord_an4 = node_coord_an4(face_node_uniqueid_an4,:);
face_node_id_an5 = tri_node_an5(mask_objface_an5, :);
face_node_uniqueid_an5 = unique(face_node_id_an5);
face_node_coord_an5 = node_coord_an5(face_node_uniqueid_an5,:);
% ----- construct triangle on face_an5 following the correspondence ------
map = containers.Map(face_node_uniqueid_an4, x_to_y);
face_tri_5from4 = zeros(size(face_node_id_an4));
for i = 1:size(face_tri_5from4, 1)
    for j = 1:3
        face_tri_5from4(i,j) = map(face_node_id_an4(i,j));
    end
end
% ----- calculate triangle edge lengths ------
face_triedges_5from4 = zeros(size(face_tri_5from4));
face_triedges_5from4(:,1) = vecnorm(face_node_coord_an5(face_tri_5from4(:,1), :) - face_node_coord_an5(face_tri_5from4(:,2), :), 2, 2);
face_triedges_5from4(:,2) = vecnorm(face_node_coord_an5(face_tri_5from4(:,2), :) - face_node_coord_an5(face_tri_5from4(:,3), :), 2, 2);
face_triedges_5from4(:,3) = vecnorm(face_node_coord_an5(face_tri_5from4(:,3), :) - face_node_coord_an5(face_tri_5from4(:,1), :), 2, 2);
% ----- get triangles with long edges ------
distort_tri = any(face_triedges_5from4 > 6, 2);
plot_node_an4 = node_coord_1(face_node_id_an4(distort_tri, :), :);
plot_node_an5 = face_node_coord_an5(face_tri_5from4(distort_tri, :), :);

% ----- visualize grain face ------
load('node_coord_180822.mat');
color1 = [0, 0.4470, 0.7410];
color2 = [0.9290, 0.6940, 0.1250];
%%
trisurf(face_node_id_an4, node_coord_1(:,1), node_coord_1(:,2), node_coord_1(:,3),'Facecolor',color1, 'Facealpha', 0.3, 'edgealpha', 0.3);
hold on
trisurf(face_node_id_an5, node_coord_2(:,1), node_coord_2(:,2), node_coord_2(:,3),'Facecolor',color2, 'Facealpha', 0.3, 'edgealpha', 0.3);
rotate3d on
% ----- visualize distorted triangles ------
scatter3(plot_node_an4(:,1), plot_node_an4(:,2), plot_node_an4(:,3), 20, 'filled', 'MarkerFaceColor',color1, 'MarkerEdgeColor',color1);
scatter3(plot_node_an5(:,1), plot_node_an5(:,2), plot_node_an5(:,3), 20, 'filled', 'MarkerFaceColor',color2, 'MarkerEdgeColor',color2);
plot3([plot_node_an4(1:3,1), plot_node_an5(1:3,1)], [plot_node_an4(1:3,2), plot_node_an5(1:3,2)], [plot_node_an4(1:3,3), plot_node_an5(1:3,3)], 'k', 'LineWidth', 1);

  
    
    
    
    
    
    
