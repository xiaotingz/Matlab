function [dist_proj_abs, dist_proj_sign, fit_ratio] = calcFaceMigByPillarHeight(features, eps_min_angle, x_to_y, tri_node_1, y_to_x, tri_node_2)
% ############################################################################
% Input
%     - features = [n+n', 10], [face_id, node_id, coordinates, normals, curves, cluster_id] 
%         including nodes of both faces
%     - x_to_y = [n, 1], a mapping from nodes_1 to nodes_2
%         this is to calculate direct node_pair_distance, for determining good median planes.
%     - tri_node_1 = [m+m', 3]
%         node_id for nodes on the grain face of previous states. 
%     - eps, scalar, 
%         the tolarence for pillar quality, here pillar quality is good if
%         top and bottom faces are of regular shape. esp = min(triangle_inner_angle)
% Output
%     - dist_proj_sign & dist_proj_abs
%         If interwining shouldn't be considered as migration, use
%         mig_sign. If interwining should be consider as migration, use
%         mig_abs
%     - fit_ratio
%        fit_ratio % of the triangles contributed to the pillar height model.
% Notes
%     - Usually, working from small face towards large face, or working from large toward small
%       , is not important since the distances would be similar so is just a scaling factor.
%       However, it's should work from large face towards small face in this pillar model, 
%       because big towards small means larger chance to find good pillars.
%     - This script can be compared to the other three projection methods, see also:
%       MigrationBySVMPlaneProj.m,  MigrationByLinearRegPlaneProj.m and MigrationByLocalNormProj.m
% ############################################################################
% ----------------------- load debug data -----------------------
% tri_node_1 = facetri_nodeid_an4;
% tri_node_2 = facetri_nodeid_an5;
% features, x_to_y, eps
% ---------------------------------------------------------------
% ##################################### Data Preparation #####################################
% """
% - If nodes correspondences are not ont-to-one, should work from the larger face towards to 
%   the smaller face, because there would be a higher chance to find good pillars. 
% - If nodes do have one-to-one corresp, then working either way is the same.
% """
% ----------------- convert the global_node_id in tri_node_ to local_node_id -----------------
mask_on_face1 = (features(:,1)==1);
mask_on_face2 = (features(:,1)==2);
node_id_global_1 = features(mask_on_face1, 2);
node_id_global_2 = features(mask_on_face2, 2);
map_1_local = containers.Map(node_id_global_1, (1:length(x_to_y)));
map_2_local = containers.Map(node_id_global_2, (1:length(y_to_x)));
map_1_corresp = containers.Map(node_id_global_1, x_to_y);
map_2_corresp = containers.Map(node_id_global_2, y_to_x);
tri_node_1_corresp = zeros(size(tri_node_1));
for i = 1:size(tri_node_1, 1)
    for j = 1:3
        tri_node_1_corresp(i,j) = map_1_corresp(tri_node_1(i,j));
        tri_node_1(i,j) = map_1_local(tri_node_1(i,j));
    end
end
tri_node_2_corresp = zeros(size(tri_node_2));
for i = 1:size(tri_node_2, 1)
    for j = 1:3
        tri_node_2_corresp(i,j) = map_2_corresp(tri_node_2(i,j));
        tri_node_2(i,j) = map_2_local(tri_node_2(i,j));
    end
end


%  -----------------  working from an4 to an5  -----------------
if nargin == 4 || ( nargin == 6 && length(x_to_y) >= length(y_to_x) )
    tri_node = tri_node_1;
    tri_node_corresp = tri_node_1_corresp;
    coords = features(mask_on_face1, 3:5);
    coords_other = features(mask_on_face2, 3:5);
    normals = features(mask_on_face1, 6:8);
    num_tris = size(tri_node, 1);
    
    % ----- if one-to-one corresp, clear triangles -----
    mask = all(tri_node_1_corresp>0, 2);
    tri_node = tri_node(mask, :);
    tri_node_corresp = tri_node_corresp(mask, :);

%  -----------------  working from an5 to an4  -----------------
else
    tri_node = tri_node_2;
    tri_node_corresp = tri_node_2_corresp;
    coords = features(mask_on_face2, 3:5);
    coords_other = features(mask_on_face1, 3:5);
    normals = features(mask_on_face2, 6:8);
    num_tris = size(tri_node, 1);
    
end



% ##################################### Find Good Pillars & Calculate Pillar Height #####################################
% """
% - Pillar quality: if the both the bottom and top face of the pillar are 
%   triangles of regular shape, mark as good quality.
% - Migration Sign: make a vector mig_avg = node_an5 - node_an4. If the direction
%   of mig_avg is consistent node_normal_an4, then sign = +.
% """
num_good_pillars = 0;
dist_proj_sign = 0;
dist_proj_abs = 0;
for i = 1:size(tri_node, 1)
    % ----------------- check pillar quality -----------------
    tri_coords = coords(tri_node(i,:),:);
    tri_coords_corresp = coords_other(tri_node_corresp(i, :),:);
    minangle_tri = calcTriangleMinAngle(tri_coords);
    minangle_tri_corresp = calcTriangleMinAngle(tri_coords_corresp);

    % ----------------- if good pillar, calc height -----------------
    if minangle_tri > eps_min_angle && minangle_tri_corresp > eps_min_angle
        num_good_pillars = num_good_pillars + 1;
        h = calcHeightTruncatedConeModel(tri_coords, tri_coords_corresp);
        
        dist_proj_abs = dist_proj_abs + h;
        mig_avgvec = sum(tri_coords - tri_coords_corresp)/3;
        normal_avg = sum(normals(tri_node(i,:),:))/3;
        ang_diff = abs( atan2d(norm(cross(mig_avgvec,normal_avg)), dot(mig_avgvec,normal_avg)) );
        if ang_diff < 90
            dist_proj_sign = dist_proj_sign + h;
        else
            dist_proj_sign = dist_proj_sign - h;
        end
    end
end

dist_proj_abs = dist_proj_abs / num_good_pillars;
dist_proj_sign = abs(dist_proj_sign) / num_good_pillars;
fit_ratio = num_good_pillars / num_tris;


end

