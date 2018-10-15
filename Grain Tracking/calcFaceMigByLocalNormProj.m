% function face_migration = MigrationByLocalNormProj(nodes_an4, nodes_an5, x_to_y)
% ############################################################################
% Input
%     - nodes_an4 = [n, 3], coordinates of the unique mesh nodes in an4
%     - nodes_an5 = [m, 3], coordinates of the unique mesh nodes in an5
%     - x_to_y = [n, 1], a mapping from nodes_an4 to nodes_an5
% Output
%     - face_migration, scalar, the average migration distance of the grain face
% Notes
%     - This script can be compared to the other two projection methods
%     also: MigrationBySVMPlaneProj.m and MigrationByPillarHeight.m
% ############################################################################
% ----------------------- load debug data -----------------------

% ---------------------------------------------------------------