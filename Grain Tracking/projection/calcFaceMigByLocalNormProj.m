function [mig_1_abs, mig_2_abs, mig_1_sign, mig_2_sign, dist_left] = calcFaceMigByLocalNormProj(features, x_to_y)
% ############################################################################
% Input
%     - features = [n+n', 10], [face_id, node_id, coordinates, normals, curves, cluster_id] 
%         see calcMigration.m for details
% Output
%     - mig_normal_proj = [mig_1_sign, mig_2_sign, mig_1_abs, mig_2_abs]
%         If interwining shouldn't be considered as migration, use
%         mig_sign. If interwining should be consider as migration, use
%         mig_abs
%     - move_left
%         Indicator variable, if the grain face has moved towards the left
%         grain. 
% Notes
%     - This script can be compared to the other two projection methods
%     also: MigrationBySVMPlaneProj.m, MigrationByLinearRegPlaneProj.m, and MigrationByPillarHeight.m
% ############################################################################

% ##### Prepare Data #####
% ----- If correspondence from full optimal transport  -----
if min(x_to_y) > 0
    coord_1 = features(features(:,1)==1, 3:5);
    coord_2 = features(features(:,1)==2, 3:5);
    coord_1_corresp = coord_2(x_to_y, :); 
    normal_1 = features(features(:,1)==1, 6:8);
    normal_2 = features(features(:,1)==2, 6:8);
    normal_1_corresp = normal_2(x_to_y, :); 
    curv_1 = features(features(:,1)==1, 9);
% ----- If correspondence from nearest tracking -----
else
    mask = x_to_y > 0;
    x_to_y = x_to_y(mask);
    coord_1 = features(features(:,1)==1, 3:5);
    coord_2 = features(features(:,1)==2, 3:5);
    coord_1 = coord_1(mask, :);
    coord_1_corresp = coord_2(x_to_y, :); 
    normal_1 = features(features(:,1)==1, 6:8);
    normal_2 = features(features(:,1)==2, 6:8);
    normal_1 = normal_1(mask, :);
    normal_1_corresp = normal_2(x_to_y, :); 
    curv_1 = features(features(:,1)==1, 9);
    curv_1 = curv_1(mask);
end

% ##### Migration Vector #####
mig_vec = coord_1_corresp - coord_1;


% ##### Migration Sign #####
is_convex = sign(curv_1);
angdiff_mig_normal = abs(atan2d(norm(cross(mig_vec,normal_1, 2)), dot(mig_vec,normal_1, 2)));
mig_along_normal = sign(angdiff_mig_normal - 90);
sign_proj = is_convex .* mig_along_normal;

% ##### Project Migration to Triangle Normals #####
% migration_1 = sum(mig_vec .* normal_1, 2);
% migration_2 = sum(mig_vec .* normal_1_corresp, 2);
migration_1 = dot(mig_vec, normal_1, 2);
migration_2 = dot(mig_vec, normal_1_corresp, 2);

% ##### Return Migration Value #####
% mig_1_sign = migration_1;
% mig_2_sign = migration_2;
% mig_1_abs = abs(migration_1);
% mig_2_abs = abs(migration_2);
mig_1_abs = sum(abs(migration_1))/length(migration_1);
mig_2_abs = sum(abs(migration_2))/length(migration_2);
mig_1_sign = sum(migration_1 .* sign_proj)/length(migration_1);
mig_2_sign = sum(migration_2 .* sign_proj)/length(migration_2);

% ##### Direction of migration #####
dist_left = sum(migration_1)/size(migration_1, 1);
% if dist_left > eps
%     move_left = 1;
% elseif dist_left < - eps
%     move_left = -1;
% else
%     move_left = 0;
% end

end





