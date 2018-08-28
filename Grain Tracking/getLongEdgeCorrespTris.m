function [longedge_tri_node_an4, longedge_tri_node_an5, longestedge] = getLongEdgeCorrespTris(obj_facelabel_an4, obj_facelabel_an5, x_to_y, thres)
% ###########################################################################
% * Input 
%     - obj_facelabel_ = [2,1], facelabel of the objective face. 
%         The full facelabel_list is returned by TrackUniqueFace.m
%     - x_to_y = [n, 1], the correspondence
%         Returned by main_TrackNodes.m or solveNodeCorresp.m
%     - thres, a scalar value defining which triangle edges are long
%  * Output
%     - longedge_tri_node_ = [m, 3], node_id for the triangles which have long edges
%         _an5 are not mesh nodes that form mesh triangles in an5, but are the an5 nodes that
%         corresponds to the an4 mesh triangle nodes
%     - longestedge = [m, 1], the length of the longest triangle edge for the triangles in an5
%  * Notes
%     - This file is use together with visualizeFace.m obj='distort_tri' to
%     see if the corresponding triangles have distorted shape. If they do,
%     the volume/area model for migration distance are not valid. If the
%     triangles are not distorted, the volume/area model for migration
%     distance is valid.
% ###########################################################################
% ------------------ load data for debug --------------------
% load('180822_FaceCorresp');
% load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822.mat');
% idx = 1;
% obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
% obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
% x_to_y = X_to_Y{idx};
% thres = 6;
% -----------------------------------------------------------
load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822.mat');
color1 = [0, 0.4470, 0.7410];
color2 = [0.9290, 0.6940, 0.1250];

% ##### get objective triangles on the objective face #####
mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));

% ##### get id and coord of the objective triangles #####
face_tri_an4 = tri_node_an4(mask_objface_an4, :);
face_node_uniqueid_an4 = unique(face_tri_an4);
face_node_coord_an4 = node_coord_an4(face_node_uniqueid_an4,:);
face_tri_an5 = tri_node_an5(mask_objface_an5, :);
face_node_uniqueid_an5 = unique(face_tri_an5);
face_node_coord_an5 = node_coord_an5(face_node_uniqueid_an5,:);

% ##### construct triangle on face_an5 following the correspondence #####
map = containers.Map(face_node_uniqueid_an4, x_to_y);
face_tri_5from4 = zeros(size(face_tri_an4));
for i = 1:size(face_tri_5from4, 1)
    for j = 1:3
        face_tri_5from4(i,j) = map(face_tri_an4(i,j));
    end
end

% ##### calculate triangle edge lengths #####
face_triedges_5from4 = zeros(size(face_tri_5from4));
face_triedges_5from4(:,1) = vecnorm(face_node_coord_an5(face_tri_5from4(:,1), :) - face_node_coord_an5(face_tri_5from4(:,2), :), 2, 2);
face_triedges_5from4(:,2) = vecnorm(face_node_coord_an5(face_tri_5from4(:,2), :) - face_node_coord_an5(face_tri_5from4(:,3), :), 2, 2);
face_triedges_5from4(:,3) = vecnorm(face_node_coord_an5(face_tri_5from4(:,3), :) - face_node_coord_an5(face_tri_5from4(:,1), :), 2, 2);

% ##### get triangles with long edges #####
distort_tri = any(face_triedges_5from4 > thres, 2);
longedge_tri_node_an4 = face_tri_an4(distort_tri, :);
longedge_tri_node_an5 = face_node_uniqueid_an5(face_tri_5from4(distort_tri, :));
longestedge = max(face_triedges_5from4(any(face_triedges_5from4 > 6, 2), :),[], 2);

end