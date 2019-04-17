% ############################################################################
% Sampling Points: K-means
% ############################################################################
load('node_coord.mat');
load('4face_node_infos.mat')
obj_face = face_node_info_4;
obj_nodecoordsmooth_an4 = node_coordsmooth_an4(obj_face{2,1},:);
obj_nodecoordsmooth_an5 = node_coordsmooth_an5(obj_face{2,2},:);
[~,kmeans_c_an4] = kmeans(obj_nodecoordsmooth_an4, 10);
[~,kmeans_c_an5] = kmeans(obj_nodecoordsmooth_an5, 10);

% visualizeFace(obj_face, kmeans_c_an4, kmeans_c_an5)

nodes_kmeans_an4 = projectKmeans(kmeans_c_an4, obj_nodecoordsmooth_an4);
nodes_kmeans_an5 = projectKmeans(kmeans_c_an5, obj_nodecoordsmooth_an5);

% ##### Local Search to Optimize Node Positions #####
% --- search on face_an4 ---

    
    
% --- search on face_an5 ---





% ############################################################################
% Sampling Points: max(min(distance))
% ############################################################################
