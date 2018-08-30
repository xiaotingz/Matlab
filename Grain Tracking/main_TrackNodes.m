% % ############################################################################
% % Solve Matching
% % * Notes
% %     - In this file, the coordinates on face are not stored because the
% %     memory cost is too expensive. However, the face coordinates can be obtained by
% %     running getSingleFaceNodes.m
%       - Three strategies
%           (1) Solve corresp between all nodes directly
%           (2) First find the disconnected subgraphs, solve corresp between subgraphs then solve nodes corresp
%           (3) Sample a subset of points instead of using all points
%           The (2) and (3) strategy can be combined
% % ############################################################################
% file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
% file_an5 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d');
% load('look_up_table_an4_an5.mat')
% facelabel_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
% facelabel_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
%     % -- NOTE triNodes are indexes starting from zero 
% tri_node_an4 = 1 + double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
% tri_node_an5 = 1 + double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
% node_coord_an4 = double(h5read(file_an4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
% node_coord_an5 = double(h5read(file_an5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
% mask_inner_an4 = all(facelabel_an4 > 0, 2);
% mask_inner_an5 = all(facelabel_an5 > 0, 2);
% facelabel_an4 = facelabel_an4(mask_inner_an4, :);
% facelabel_an5 = facelabel_an5(mask_inner_an5, :);
% facelabel_an4 = sort(facelabel_an4, 2);
% facelabel_an5 = sort(facelabel_an5, 2);
% tri_node_an4 = tri_node_an4(mask_inner_an4, :);
% tri_node_an5 = tri_node_an5(mask_inner_an5, :);
% 
% % --- tracked_unique_face = trackUniqueFace(file_an4, file_an5, look_up_table, complete) ---
% [tracked_uniqueface_an4, tracked_uniqueface_an5] = trackUniqueFace(file_an4, file_an5, look_up_table, 1);
% n = length(tracked_uniquean4);
 
% %% ############################################################################
% (1) Dicret Corresp
% % ############################################################################
% %% ##### (1.1) Corresp of The Small Faces #####
% X_to_Y = cell(n,1);
% Y_to_X = cell(n,1);
% large_face_id = [];
% for i = 1:length(tracked_uniqueface_an4)
%     obj_facelabel_an4 = tracked_uniqueface_an4(i, :);
%     obj_facelabel_an5 = tracked_uniqueface_an5(i, :);
%     
%     % ##### Get Objective Triangles on the Objective Face #####
%     mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
%     mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));
% 
%     % ##### Get Id and Coord of the Objective Triangles #####
%     face_tri_nodeid_an4 = tri_node_an4(mask_objface_an4, :);
%     face_unique_nodeid_an4 = unique(face_tri_nodeid_an4);
%     face_node_coord_an4 = face_node_coord_an4(face_unique_nodeid_an4,:);
%     face_tri_nodeid_an5 = tri_node_an5(mask_objface_an5, :);
%     face_unique_nodeid_an5 = unique(face_tri_nodeid_an5);
%     face_node_coord_an5 = face_node_coord_an5(face_unique_nodeid_an5,:);
%     
%     if (length(face_unique_nodeid_an4) < 1000) && (length(face_unique_nodeid_an5) < 1000)
%         [x_to_y, y_to_x] = solveNodeCorresp(face_node_coord_an4, face_node_coord_an5);
%         X_to_Y{i} = int32(x_to_y);
%         Y_to_X{i} = int32(y_to_x)';
%     else 
%         large_face_id = [large_face_id, i];
%         X_to_Y{i} = NaN;
%         Y_to_X{i} = NaN;
%     end
%     disp(['small faces: ', num2str(i), '/7004']);
% end
% save('180820_smallFaces.mat', 'X_to_Y', 'Y_to_X', 'large_id');
% 
% 
% %% ##### (1.2) Corresp of The Large Faces #####
% for i = 1:length(large_face_id)
%     idx = large_face_id(i);
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
%     face_node_coord_an4 = face_node_coord_an4(face_unique_nodeid_an4,:);
%     face_tri_nodeid_an5 = tri_node_an5(mask_objface_an5, :);
%     face_unique_nodeid_an5 = unique(face_tri_nodeid_an5);
%     face_node_coord_an5 = face_node_coord_an5(face_unique_nodeid_an5,:);
%     
%     [x_to_y, y_to_x] = solveNodeCorresp(face_node_coord_an4, face_node_coord_an5);
%     X_to_Y_large{i} = int32(x_to_y);
%     Y_to_X_large{i} = int32(y_to_x);
%     
%     disp(['large faces: ', num2str(i), '/', num2str(large_sizee)]);
% end
% save('180820_largeFaces.mat', 'X_to_Y_large', 'Y_to_X_large');
% 

%% ############################################################################
% (2) First Find Disconnected Subgraphs 
%      - Not every subgraph will find its corresp. 
%           length(subgraph_pair) = min(length(subgraph_an4), length(subgraph_an5))
% % ############################################################################
load('180822_FaceCorresp');
load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822.mat');

face_piecewise = [];
subgraph_corresp_list = {};
idx_face_piecewise = 1;
for idx = 1:length(tracked_uniqueface_an4)
    obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
    obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);

    % ##### Get Objective Triangles on the Objective Face #####
    mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
    mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));

    % ##### Get Id and Coord of the Objective Triangles #####
    face_tri_nodeid_an4 = tri_node_an4(mask_objface_an4, :);
    face_unique_nodeid_an4 = unique(face_tri_nodeid_an4);
    face_node_coord_an4 = node_coord_an4(face_unique_nodeid_an4,:);
    face_tri_nodeid_an5 = tri_node_an5(mask_objface_an5, :);
    face_unique_nodeid_an5 = unique(face_tri_nodeid_an5);
    face_node_coord_an5 = node_coord_an5(face_unique_nodeid_an5,:);

    % ##### Find Disconnected Subgraphs and Solve Subgraph Correp #####
    [subgraph_an4, subgraph_an5] = findSubgraph(face_unique_nodeid_an4, face_unique_nodeid_an5, face_tri_nodeid_an4, face_tri_nodeid_an5);

    if length(unique(subgraph_an4)) > 1 || length(unique(subgraph_an5)) > 1
        face_piecewise = [face_piecewise, i];
    end
    
end

%% 
% % ############################################################################
% % Check: Visualize Subgraph
% % ############################################################################
color1 = [0, 0.4470, 0.7410];
color2 = [0.9290, 0.6940, 0.1250];
color3 = [0.8500, 0.3250, 0.0980];
color4 = [0.4660, 0.6740, 0.1880];
colors = [color1; color2; color3; color4];

idx = face_piecewise(randi(length(face_piecewise)));
obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);

% ##### Get Objective Triangles on the Objective Face #####
mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));

% ##### Get Id and Coord of the Objective Triangles #####
face_tri_nodeid_an4 = tri_node_an4(mask_objface_an4, :);
face_unique_nodeid_an4 = unique(face_tri_nodeid_an4);
face_node_coord_an4 = node_coord_an4(face_unique_nodeid_an4,:);
face_tri_nodeid_an5 = tri_node_an5(mask_objface_an5, :);
face_unique_nodeid_an5 = unique(face_tri_nodeid_an5);
face_node_coord_an5 = node_coord_an5(face_unique_nodeid_an5,:);

[subgraph_an4, subgraph_an5] = findSubgraph(face_unique_nodeid_an4, face_unique_nodeid_an5, face_tri_nodeid_an4, face_tri_nodeid_an5);
subgraph_corresp = solveSubgraphCorresp(subgraph_an4, subgraph_an5, face_node_coord_an4, face_node_coord_an5);


trisurf(face_tri_nodeid_an4, node_coord_an4(:,1), node_coord_an4(:,2), node_coord_an4(:,3),'Facecolor',color1, 'Facealpha', 0.3, 'edgealpha', 0.3);
hold on
rotate3d on
trisurf(face_tri_nodeid_an5, node_coord_an5(:,1), node_coord_an5(:,2), node_coord_an5(:,3),'Facecolor',color2, 'Facealpha', 0.3, 'edgealpha', 0.3);
if length(unique(subgraph_an4)) < length(unique(subgraph_an5))
    for i = 1:length(unique(subgraph_an4))
        scatter3(face_node_coord_an4((subgraph_an4==i), 1), face_node_coord_an4((subgraph_an4==i), 2), face_node_coord_an4((subgraph_an4==i), 3), ...
            10, 'filled', 'MarkerFaceColor',colors(i, :), 'MarkerEdgeColor',colors(i, :));
        scatter3(face_node_coord_an5((subgraph_an5==subgraph_corresp(i)), 1), face_node_coord_an5((subgraph_an5==subgraph_corresp(i)), 2), face_node_coord_an5((subgraph_an5==subgraph_corresp(i)), 3), ...
            10, 'filled', 'MarkerFaceColor',colors(i, :), 'MarkerEdgeColor',colors(i, :));
    end
else
    for i = 1:length(unique(subgraph_an5))
        scatter3(face_node_coord_an4((subgraph_an4==subgraph_corresp(i)), 1), face_node_coord_an4((subgraph_an4==subgraph_corresp(i)), 2), face_node_coord_an4((subgraph_an4==subgraph_corresp(i)), 3), ...
            10, 'filled', 'MarkerFaceColor',colors(i, :), 'MarkerEdgeColor',colors(i, :));
        scatter3(face_node_coord_an5((subgraph_an5==i), 1), face_node_coord_an5((subgraph_an5==i), 2), face_node_coord_an5((subgraph_an5==i), 3), ...
            10, 'filled', 'MarkerFaceColor',colors(i, :), 'MarkerEdgeColor',colors(i, :));
    end
end
hold off
print(['pair_', num2str(idx), '_subgraph'], '-dpng','-r300')





% % ############################################################################
% % Visualize Face Correspondences
% % ############################################################################
% idx = 6;
% % idx = 3102;
% 
% % ----- get the object face triangles and nodes -----
% obj_facelabel_an4 = tracked_uniquean4(idx, :);
% obj_facelabel_an5 = tracked_uniquean5(idx, :);
% node_info = getSingleFaceNodes(file_an4, obj_facelabel_an4, file_an5, obj_facelabel_an5);
% 
% % ----- plot with alpha shape-----
% figure(1)
% visualizeFace(node_info, X_to_Y{idx})
% 
% % ----- alpha shape -----
% % shp = alphaShape([node_info{4,1}; node_info{4,2}]);
% % shp.Alpha = 16.4938;
% % plot(shp)
% % alpha(.5)
% % rotate3d on
% volume = shp.volume;


