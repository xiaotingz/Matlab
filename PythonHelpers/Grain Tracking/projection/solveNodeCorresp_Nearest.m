function x_to_y = solveNodeCorresp_Nearest(X, Y)
% ##########################################################################
% * Input
%     - X = [m, 3], node coordinates on the before state face
%     - Y = [n, 3], node coordinates on the after state face
% * Output
%     - x_to_y = [m, 1], for each node in X, id of the corresponding node in Y
% * NOTE
%     - Make sure samples are aligned by changing ORIGIN.
%     - This script finds the nearest nodes for nodes on the smaller face.
% ##########################################################################
% % ----------------------- load debug data -----------------------
%     idx = 1500;
%     obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
%     obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
%     
%     % ##### Get Objective Triangles on the Objective Face #####
%     mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
%     mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));
% 
%     % ##### Get Id and Coord of the Objective Triangles #####
%     face_tri_nodeid_an4 = tri_node_an4(mask_objface_an4, :);
%     face_unique_nodeid_an4 = unique(face_tri_nodeid_an4);
%     face_node_coord_an4 = node_coord_an4(face_unique_nodeid_an4,:);
%     face_tri_nodeid_an5 = tri_node_an5(mask_objface_an5, :);
%     face_unique_nodeid_an5 = unique(face_tri_nodeid_an5);
%     face_node_coord_an5 = node_coord_an5(face_unique_nodeid_an5,:);
% 
% 
% X = face_node_coord_an4;
% Y = face_node_coord_an5;
% % ---------------------------------------------------------------

m = size(X,1);
n = size(Y,1);

% ####################################### Calc Node Pair Distance #######################################
diff_tmp1 = repmat(X, 1, 1, n);
diff_tmp1 = permute(diff_tmp1, [1,3,2]);
diff_tmp2 = repmat(Y, 1, 1, m);
diff_tmp2 = permute(diff_tmp2, [3,1,2]);
diff = diff_tmp1 - diff_tmp2;
norm_diff = sum(diff .* diff, 3);


% % ####################################### Keep Mininum Distance Pairs ####################################### 
% % ------------------------ Global View ------------------------ 
% % """
% % Work globaly. The problem is that at regions with large shape change,
% % good minimum distance pairs are rare. Most nodes are just systematically
% % off. 
% % """
% [x_to_y_tmpval, x_to_y_tmpidx] = min(norm_diff, [], 2);
% [y_to_x_tmpval, y_to_x_tmpidx] = min(norm_diff, [], 1);
% corresp_tmp1 = [(1:m)', x_to_y_tmpidx, x_to_y_tmpval];
% corresp_tmp2 = [y_to_x_tmpidx', (1:n)', y_to_x_tmpval'];
% corresp = intersect(corresp_tmp1(:, 1:2), corresp_tmp2(:, 1:2), 'rows');
% 
% % ##### Convert the pair to idx of initial face #####
% x_to_y = NaN([m, 1]);
% x_to_y(corresp(:,1)) = corresp(:,2);


% ------------------------ Local View ------------------------ 
% """
% Work locally, from the smaller face. 
% The procedure is indeed two-way: 1) For every node on the smaller face, find
% its minimum distance corresp node_larger_faceon. Then 2) gather the unique 
% node_larger_faceon. For each node_larger_faceon, if it has several min
% distance corresp nodes on an4 face, keep the one with minimum distance. 
% """
if m <= n
    [x_to_y_tmpval, x_to_y_tmpidx] = min(norm_diff, [], 2);
    corresp_tmp = [(1:m)', x_to_y_tmpidx, x_to_y_tmpval];
    corresp_unique = unique(corresp_tmp(:,2));
    corresp_unique = [zeros(length(corresp_unique), 1), corresp_unique];
    for i = 1:size(corresp_unique, 1)
        mask = (corresp_tmp(:,2) == corresp_unique(i, 2));
        candidates = corresp_tmp(mask, :);
        [~, idx] = min(candidates(:,3));
        corresp_unique(i,1) = candidates(idx, 1);
    end
    x_to_y = NaN([m, 1]);
    x_to_y(corresp_unique(:,1)) = corresp_unique(:,2);
else
    [y_to_x_tmpval, y_to_x_tmpidx] = min(norm_diff, [], 1);
    corresp_tmp = [y_to_x_tmpidx', (1:n)', y_to_x_tmpval'];
    corresp_unique = unique(corresp_tmp(:,1));
    corresp_unique = [corresp_unique, zeros(length(corresp_unique), 1)];
    for i = 1:size(corresp_unique, 1)
        mask = (corresp_tmp(:,1) == corresp_unique(i, 1));
        candidates = corresp_tmp(mask, :);
        [~, idx] = min(candidates(:,3));
        corresp_unique(i,2) = candidates(idx, 2);
    end
    x_to_y = NaN([n, 1]);
    x_to_y(corresp_unique(:,1)) = corresp_unique(:,2);
end




% % ####################################### Visual Checks ####################################### 
% face_node_info = getSingleFaceNodes(obj_facelabel_an4, obj_facelabel_an5);
% figure
% visualizeFace(face_node_info, x_to_y, corresp_unique);


% face_node_info_tmp = face_node_info;
% tmp = face_node_info_tmp{4,2};
% face_node_info_tmp{4,2} = face_node_info_tmp{4,1};
% face_node_info_tmp{4,1} = tmp;
% tmp = face_node_info_tmp{1,2};
% face_node_info_tmp{1,2} = face_node_info_tmp{1,1};
% face_node_info_tmp{1,1} = tmp;
% tmp = face_node_info_tmp{5,2};
% face_node_info_tmp{5,2} = face_node_info_tmp{5,1};
% face_node_info_tmp{5,1} = tmp;
% corresp = [1:length(y_to_x_tmpidx); y_to_x_tmpidx]';
% visualizeFace(face_node_info_tmp, y_to_x_tmpidx', corresp);

end



