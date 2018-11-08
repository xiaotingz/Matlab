function [mig_1_sign, mig_2_sign, mig_1_abs, mig_2_abs] = calcFaceMigByLocalNormProj(features)
% ############################################################################
% Input
%     - features = [n+m, 9], [face_id, node_id, coordinates, normals, cluster_id] 
%         see calcMigration.m for details
% Output
%     - mig_normal_proj = [mig_1_sign, mig_2_sign, mig_1_abs, mig_2_abs]
%         If interwining shouldn't be considered as migration, use
%         mig_sign. If interwining should be consider as migration, use
%         mig_abs
% Notes
%     - This script can be compared to the other two projection methods
%     also: MigrationBySVMPlaneProj.m and MigrationByPillarHeight.m
% ############################################################################

% ##### Prepare Data #####
coord_1 = features(features(:,1)==1, 3:5);
coord_1_corresp = features(features(:,1)==2, 3:5);
% coord_2 = features(features(:,1)==2, 3:5);
% coord_1_corresp = coord_2(x_to_y, :); 
normal_1 = features(features(:,1)==1, 6:8);
normal_1_corresp = features(features(:,1)==2, 6:8);
% normal_2 = features(features(:,1)==2, 6:8);
% normal_1_corresp = normal_2(x_to_y, :); 


% ##### The Direct Migration Vector #####
migration_direct = coord_1 - coord_1_corresp;

% ##### Project Migration to Triangle Normals #####
migration_1 = sum(migration_direct .* normal_1, 2);
migration_2 = sum(migration_direct .* normal_1_corresp, 2);

% ##### Return Migration Value #####
mig_1_sign = migration_1;
mig_2_sign = migration_2;
mig_1_abs = abs(migration_1);
mig_2_abs = abs(migration_2);
% mig_1_sign = sum(migration_1)/length(migration_1);
% mig_2_sign = sum(migration_2)/length(migration_2);
% mig_1_abs = sum(abs(migration_1))/length(migration_1);
% mig_2_abs = sum(abs(migration_2))/length(migration_2);
% mig_normal_proj = [mig_1_sign, mig_2_sign, mig_1_abs, mig_2_abs];

end





