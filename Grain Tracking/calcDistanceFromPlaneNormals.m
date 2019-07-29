% function [dists, std] = calcDistanceFromPlaneNormals(file, faces, target_normals, eps_curv, eps_area, eps_min_ang)
% ###########################################################################
% * Input
%     - faces = [n,2] = [label_a, label_b]
%         Facelabels of the faces of interest
%     - target_normals = [m, 3]
%         Plane normals of interest, e.g. [1,1,1], [1,1,0], [1,0,0]
%     - eps_
%         Filter bad triangles, this requires corresponding D3D filters to be run.
% * Output
%     - dists = [n,1] 
%         Weighted distance of triangle normals from target_normals
%     - std = [n, 1]
%         Standard deviation from the distances.
% * Note 
%     - This file is to be used in featurePrepOther.m
%     - This function aims to characterize the plane normal of a grain
%     face in a few scalar numbers. Scalars because we don't want deep NN.
%     A grain face include many different normal directions, so the scalars
%     are chosen as the avg(dist(triangle_normals, target_normals)) and 
%     std(dist(triangle_normals, target_normals)).
% ###########################################################################
% ------------------ load data for debug --------------------
file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin3_Hsmooth.dream3d';
eps_curv = 1;
eps_area = 7;
eps_min_ang = 10;
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190624_tracked_faces_full.mat', 'tracked_uniqueface_an4_full')
faces = tracked_uniqueface_an4_full;
clear tracked_uniqueface_an4_full
target_normals = [1,1,1;1,1,0;1,0,0;1,1,2];
% -----------------------------------------------------------




















