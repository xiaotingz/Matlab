function [distort_tri_node_an4, distort_tri_node_an5, values] = getDistortCorrespTris(obj_facelabel_an4, obj_facelabel_an5, x_to_y, obj, thres)
% ###########################################################################
% * Input 
%     - obj_facelabel_ = [2,1], facelabel of the objective face. 
%         The full facelabel_list is returned by TrackUniqueFace.m
%     - x_to_y = [n, 1], the correspondence
%         Returned by main_TrackNodes.m or solveNodeCorresp.m
%     - thres
%         a scalar value defining what is a long edge,
%         or what is a small angle
%     - obj = 'long_edge' or 'min_angle_diff'
%  * Output
%     - distort_tri_node_ = [m, 3], node_id for the triangles which have long
%     edges or small inner angle
%         _an5 are not mesh nodes that form mesh triangles in an5, but are the an5 nodes that
%         corresponds to the an4 mesh triangle nodes
%     - values = [m, 1], 
%         the length of the longest triangle edge for the triangles in an5,
%         or min_angle_diff, = min_angle_an5 - min_angle_an4
%  * Notes
%     - This file is use together with visualizeFace.m to see if the corresponding
%      triangles have distorted shape. 
%     - Run this file for all faces in TrackAnalysis.m
%     - If triangles are distorted, the volume/area model for migration distance have high error.
% ###########################################################################
% ------------------ load data for debug --------------------
% load('180822_FaceCorresp');
% load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822.mat');
% idx = 171;
% obj_facelabel_an4 = tracked_uniqueface_an4(idx, :);
% obj_facelabel_an5 = tracked_uniqueface_an5(idx, :);
% x_to_y = X_to_Y{idx};
% obj = 'min_angle_diff';
% thres = 20;
% -----------------------------------------------------------
load('/Users/xiaotingzhong/Dropbox/grainTracking_forCluster/180822.mat');

% ##### get objective triangles on the objective face #####
mask_objface_an4 = (facelabel_an4(:,1) == obj_facelabel_an4(1) & facelabel_an4(:,2) == obj_facelabel_an4(2));
mask_objface_an5 = (facelabel_an5(:,1) == obj_facelabel_an5(1) & facelabel_an5(:,2) == obj_facelabel_an5(2));

% ##### get id and coord of the objective triangles #####
face_tri_an4 = tri_node_an4(mask_objface_an4, :);
face_node_uniqueid_an4 = unique(face_tri_an4);
face_tri_an5 = tri_node_an5(mask_objface_an5, :);
face_node_uniqueid_an5 = unique(face_tri_an5);

% ##### construct triangle on face_an5 following the correspondence #####
map = containers.Map(face_node_uniqueid_an4, x_to_y);
face_tri_5from4 = zeros(size(face_tri_an4));
for i = 1:size(face_tri_5from4, 1)
    for j = 1:3
        face_tri_5from4(i,j) = map(face_tri_an4(i,j));
    end
end

if strcmp(obj, 'long_edge')
%     face_node_coord_an4 = node_coord_an4(face_node_uniqueid_an4,:);
    face_node_coord_an5 = node_coord_an5(face_node_uniqueid_an5,:);
    % ##### calculate triangle edge lengths #####
    face_triedges_5from4 = zeros(size(face_tri_5from4));
    face_triedges_5from4(:,1) = vecnorm(face_node_coord_an5(face_tri_5from4(:,1), :) - face_node_coord_an5(face_tri_5from4(:,2), :), 2, 2);
    face_triedges_5from4(:,2) = vecnorm(face_node_coord_an5(face_tri_5from4(:,2), :) - face_node_coord_an5(face_tri_5from4(:,3), :), 2, 2);
    face_triedges_5from4(:,3) = vecnorm(face_node_coord_an5(face_tri_5from4(:,3), :) - face_node_coord_an5(face_tri_5from4(:,1), :), 2, 2);
    % ##### get triangles with long edges #####
    mask_distort_tri = any(face_triedges_5from4 > thres, 2);
    distort_tri_node_an4 = face_tri_an4(mask_distort_tri, :);
    distort_tri_node_an5 = face_node_uniqueid_an5(face_tri_5from4(mask_distort_tri, :));
    values = max(face_triedges_5from4(any(face_triedges_5from4 > 6, 2), :),[], 2);
    
elseif strcmp(obj, 'min_angle_diff')
    min_angle_diff = - ones(length(face_tri_an4), 1);
    for i = 1:size(face_tri_an4, 1)
        tri_an4 = node_coord_an4(face_tri_an4(i,:),:);
        tri_5from4 = node_coord_an5(face_node_uniqueid_an5(face_tri_5from4(i,:)),:);
        minangle_tri_an4 = calcTriangleMinAngle(tri_an4);
        minangle_tri_an5 = calcTriangleMinAngle(tri_5from4);
        min_angle_diff(i) = abs(minangle_tri_an4 - minangle_tri_an5);
    end
    mask_distort_tri = (min_angle_diff > thres);
    distort_tri_node_an4 = face_tri_an4(mask_distort_tri,:);
    distort_tri_node_an5 = face_node_uniqueid_an5(face_tri_5from4(mask_distort_tri, :));
    values = min_angle_diff(mask_distort_tri);
end

end

% ##### Visual Check of Distort Triangles #####
% color1 = [0, 0.4470, 0.7410];
% color2 = [0.9290, 0.6940, 0.1250];
% % obj_face = getSingleFaceNodes(file_an4, obj_facelabel_an4, file_an5, obj_facelabel_an5);
% idx_tri = 9;
% visualizeFace(obj_face, x_to_y, distort_tri_node_an4(idx_tri, :), distort_tri_node_an5(idx_tri,:), 'distort_tri')
% 
% calcTriangleMinAngle(node_coord_an4(distort_tri_node_an4(idx_tri, :), :))
% calcTriangleMinAngle(node_coord_an5(distort_tri_node_an5(idx_tri, :), :))