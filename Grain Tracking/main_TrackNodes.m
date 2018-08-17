% ############################################################################
% Sampling Points: K-means
% ############################################################################
load('node_coord.mat');
load('4face_node_infos.mat')
obj_face = face_node_info_4;
obj_nodecoordsmooth_an4 = node_coordsmooth_an4(obj_face{2,1},:);
obj_nodecoordsmooth_an5 = node_coordsmooth_an5(obj_face{2,2},:);
[~,center_an4] = kmeans(obj_nodecoordsmooth_an4, 30);
[~,center_an5] = kmeans(obj_nodecoordsmooth_an5, 20);


% ############################################################################
% Sampling Points: max(min(distance))
% ############################################################################


%
% ############################################################################
% Solve Matching
% ############################################################################
X = center_an4;
Y = center_an5;
[x_to_y, y_to_x] = solveNodeCorresp(X, Y);

visualizeFace(obj_face, center_an4, center_an5, x_to_y)


