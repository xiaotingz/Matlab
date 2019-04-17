function projected_nodes = projectKmeans(kmeans_center, nodes_on_surface)
% ############################################################################
% Input
%     - nodes_on_surface = [n, 3], the mesh nodes on surface
%     - kmeans_center = [m, 3], found by running k
% Output
%     - projected_nodes = [m, 3], the surface nodes that closest to the
%     centers found by kmeans
% ############################################################################
% ----------------------- load debug data -----------------------
% kmeans_center = kmeans_c_an4;
% nodes_on_surface = obj_nodecoordsmooth_an4;
% ---------------------------------------------------------------
m = size(kmeans_center, 1);
n = size(nodes_on_surface, 1);

diff_tmp1 = repmat(nodes_on_surface, 1, 1, m);
diff_tmp1 = permute(diff_tmp1, [1,3,2]);
diff_tmp2 = repmat(kmeans_center, 1, 1, n);
diff_tmp2 = permute(diff_tmp2, [3,1,2]);
diff = diff_tmp1 - diff_tmp2;
norm_diff = sum(diff .* diff, 3);


[~, projected_nodes] = min(norm_diff, [], 1);


% % ##### Visualization Check #####
% trisurf(reshape(obj_face{2,1}, 3, [])', node_coordsmooth_an4(:,1), node_coordsmooth_an4(:,2), node_coordsmooth_an4(:,3), 'Facealpha', 0.3, 'edgealpha', 0.3);
% hold on
% rotate3d on
% scatter3(kmeans_c_an4(:,1), kmeans_c_an4(:,2), kmeans_c_an4(:,3), 20, 'filled', 'MarkerFaceColor','k', 'MarkerEdgeColor','k');
% scatter3(obj_nodecoordsmooth_an4(projected_nodes,1), obj_nodecoordsmooth_an4(projected_nodes,2), obj_nodecoordsmooth_an4(projected_nodes,3), 20, 'filled', 'MarkerFaceColor','r', 'MarkerEdgeColor','r');

end