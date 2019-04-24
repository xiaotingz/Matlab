function [triple_line, tl_info] = findTripleLines(file, eps_area, eps_curv, eps_min_ang)
% ############################################################################
% * Output
%   - triple_line = [n, 3], n is the number of triple lines in the volume, 3
%       stands for the threeGrainIDs encapsulating the triple line.
%   - tl_info = [n, 4]
%       [:, 1:3] = weighted dihedral angle along the triple line.
%       [:, 4]   = length of the triple line.
% * Illustration of order
%                     |
%           g2, da_2  |  g3, da_3
%                    / \
%                   /   \
%                 g1, da_1
% ############################################################################
% ----------------------- load debug data -----------------------
% % clear
% % file = '/Users/xiaotingzhong/Desktop/Datas/synthetic/180502_CubicSingleEquiaxedOut.dream3d';
% file = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
% % '/Volumes/XIAOTING/Ni/An4new6_fixOrigin2_smooth.dream3d';
% % '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
% % '/Volumes/XIAOTING/Ni/An5new6_Hsmooth.dream3d'
% eps_area = 7;
% eps_curv = 1;
% eps_min_ang = 10;
% ---------------------------------------------------------------

% ##### Read in data ##### 
tri_node_raw = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coords = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
fl_raw = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
node_type = double(h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType'))';
tri_curv = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)';
tri_area = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
tri_min_ang = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';
tri_normal_raw = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'),-5)';
fl_idx = (1:size(fl_raw, 1))';

% ##### Find the triangles on tripleline #####
%     """
%     - treat QuadNodes as TLnodes
%     - the triple line triangles are those which have at least 2 TLnodes or QuadNodes 
%            TLnodes=3 happens for the triangles sitting on a corner 
%     - only mesh triangles of good quality are used.
%     """
node_type(node_type==4) = 3;
mask_tl_tri = (sum(node_type(tri_node_raw)==3, 2) >= 2);
mask_good_tri = (all(fl_raw>0, 2) & abs(tri_curv) < eps_curv & ...
                tri_area < eps_area & tri_min_ang > eps_min_ang);
mask_obj_tri = mask_tl_tri & mask_good_tri;
fl = fl_raw(mask_obj_tri, :);
tri_node = tri_node_raw(mask_obj_tri, :);
tri_normal = tri_normal_raw(mask_obj_tri, :);
fl_idx = fl_idx(mask_obj_tri);

% ##### Sort the triangles into groups of 3 #####
%     """
%     - creat a TLnode list to record the triple line nodes of each triangle,
%           group triangles according to the common triple line node 
%     - special case are the triangles having 3 quad nodes, then it's difficult to tell which are the two nodes in use
%           in such case, duplicate the triangle
%     """
mask_TLnodes = (node_type(tri_node)==3);
TLnodes_raw = tri_node.*mask_TLnodes;
TLnodes_raw = sort(TLnodes_raw,2);
% TLnodes(:,1) = (1:length(TLnodes))';
mask_3quadTris = ~any(TLnodes_raw == 0, 2);
TLnodes = zeros(length(TLnodes_raw) + sum(mask_3quadTris)*2, 3);
idx = 1;
for i = 1:length(TLnodes_raw)
    if TLnodes_raw(i, 1) == 0
        TLnodes(idx, 1) = i;
        TLnodes(idx, 2:3) = TLnodes_raw(i, 2:3);
        idx = idx + 1;
    else
        TLnodes(idx:idx+2, 1) = i;
        TLnodes(idx, 2:3) = TLnodes_raw(i, 1:2);
        TLnodes(idx+1, 2:3) = TLnodes_raw(i, 2:3);
        TLnodes(idx+2, 2) = TLnodes_raw(i, 1);
        TLnodes(idx+2, 3) = TLnodes_raw(i, 3);
        idx = idx + 3;
    end
end
% ---- TLnodes = [IDinTLnodes_raw, TLnode1, TLnode2] ----
TLnodes = sortrows(TLnodes,[2,3]);

% --- if a group contain more or less than 3 triangles, ignore ---
mask_groupOf3 = zeros(length(TLnodes),1);
cnt = 1;
for i = 1:(length(TLnodes)-1)
    if TLnodes(i,2)==TLnodes(i+1,2) && TLnodes(i,3)==TLnodes(i+1,3)
        cnt = cnt + 1;
    else
        if cnt == 3
            mask_groupOf3(i-2:i) = 1;
        end
        cnt = 1;
    end
end
% --- the last else wasn't excuted ---
if cnt == 3
    mask_groupOf3(i+1-2:i+1) = 1;
end
mask_groupOf3 = logical(mask_groupOf3);

% ##### get triple lines from the triangle groups #####
% """
% TLnodes = [tri_id, tlnode_id_1, tlnode_id_2]
% Note that the triangle winding must be kept correct for the angles to be correct. 
% """
tl_tri_fl = fl(TLnodes(mask_groupOf3, 1), :);
% tl_tri_nodeid = tri_node(TLnodes(mask_groupOf3, 1), :);
tl_node = TLnodes(mask_groupOf3, 2:3);
tl_tri_normal = tri_normal(TLnodes(mask_groupOf3, 1), :);
tl_tri_fl_idx = fl_idx(TLnodes(mask_groupOf3, 1));
if rem(length(tl_tri_fl), 3) ~= 0 
    warning('ERROR, length(TLtri_FLs)/3 is not integer')
end
tlgroup_id = zeros(length(tl_tri_fl)/3, 3);
tlgroup_da = zeros(length(tl_tri_fl)/3, 3);
tlgroup_length = zeros(length(tl_tri_fl)/3, 1);
tlgroup_fl_idx = zeros(length(tl_tri_fl)/3, 3);

for i = 1:length(tl_tri_fl)/3
    fl_group = tl_tri_fl((i-1)*3+1:(i-1)*3+3, :);
    
    if length(unique(fl_group)) == 3 && size(unique(sort(fl_group, 2), 'rows'),1) == 3
        tlgroup_fl_idx(i,:) = tl_tri_fl_idx((i-1)*3+1 : (i-1)*3+3);
        tlgroup_id(i,:) = unique(fl_group);
        tl_node_coord_1 = node_coords(tl_node((i-1)*3+1, 1), :);
        tl_node_coord_2 = node_coords(tl_node((i-1)*3+1, 2), :);
        tlgroup_length(i) = norm(tl_node_coord_1 - tl_node_coord_2);
% ----------------------------------------------------------------------------------
        % """ 
        % no really working 
        % """
        % --- calculate normal from coords, and keep winding consistent ---
%         node_group = tl_tri_nodeid((i-1)*3+1 : (i-1)*3+3, :);
%         common_nodes = intersect(node_group(1,:), node_group(2,:));
%         % ------------------------------------------------------------------------
%         % debug
%         common_nodes_23 = intersect(node_group(2,:), node_group(3,:));
%         common_nodes_13 = intersect(node_group(1,:), node_group(3,:));
%         if length(common_nodes) ~= 2 || ~all(common_nodes == common_nodes_23) || ~all(common_nodes_13 == common_nodes_23)
%             warning(['problem at group ', num2str(i)]);
%         end
%         % ------------------------------------------------------------------------
%         p1 = node_coords(common_nodes(1), :);
%         p2 = node_coords(common_nodes(2), :);
%         normal_group = zeros(3,3);
%         for j = 1:3
%             p3 = node_coords(setdiff(node_group(j,:), common_nodes), :);
%             normal_group(j, :) = cross(p1-p2, p1-p3);
%             normal_group(j, :) = normal_group(j, :) / norm(normal_group(j, :));
%         end
% ----------------------------------------------------------------------------------        

        % --- calcuate dihedral_angle, and keep the order consistent with the grain_id in triple_line ---
        normal_group = tl_tri_normal((i-1)*3+1 : (i-1)*3+3, :);
        for j = 1:3
            grain_id = tlgroup_id(i, j);
            mask = any(fl_group == grain_id, 2);
            fl_grain = fl_group(mask, :);
            normal_grain = normal_group(mask, :);
            % """ two faces of a grain, one normal point inwards, the other point outwards """
            if grain_id == fl_grain(1, 1)
                normal_grain(1, :) = - normal_grain(1, :);
            end
            if grain_id == fl_grain(2, 2)
                normal_grain(2, :) = - normal_grain(2, :);
            end
            tlgroup_da(i, j) = atan2d(norm(cross(normal_grain(1,:),normal_grain(2,:))),dot(normal_grain(1,:),normal_grain(2,:)));
%             % ------------------------------------------------------------------------
%             % debug
%             if size(normal_grain, 1) ~= 2
%                 warning(['problem 2 at group ', num2str(i)]);
%             end
%             % ------------------------------------------------------------------------
        end
    else
%         disp(['length(unique(FLgroup)) ~= 3 at FLgroup=', num2str(i)]);
    end
end

% ----- clean the data ----- 
% """
% Reasons for three-triangle groups may not be good: 1. artifitial group.
% 2. Actually a four--triangle group but got incorporated because one of
% the four triangles had bad quality and excluded before grouping.
% """
mask_valid_tl = (all(tlgroup_da>0, 2) & abs(sum(tlgroup_da, 2) - 360)<0.1);
tlgroup_id = tlgroup_id(mask_valid_tl, :);
tlgroup_da = tlgroup_da(mask_valid_tl, :);
tlgroup_length = tlgroup_length(mask_valid_tl, :);

% ##### get unique triple_line identifiers and the weighted turning angles #####
% ----- first sort triple line into groups ----- 
[tlgroup_id, sort_idx] = sortrows(tlgroup_id);
tlgroup_da = tlgroup_da(sort_idx, :);
tlgroup_length = tlgroup_length(sort_idx, :);
triple_line = [];
dihedral_angle = [];
len = [];

culm_da = [0,0,0];
culm_len = 0.0;
cnt = 0;
for i = 1:length(tlgroup_id)-1
    if all(tlgroup_id(i,:) == tlgroup_id(i+1,:))
        culm_da = culm_da + tlgroup_da(i, :);
        culm_len = culm_len + tlgroup_length(i);
        cnt = cnt + 1;
    else
        culm_da = culm_da + tlgroup_da(i, :);
        culm_len = culm_len + tlgroup_length(i);
        cnt = cnt + 1;
        triple_line = [triple_line; tlgroup_id(i,:)];
        dihedral_angle = [dihedral_angle; culm_da/cnt];
        len = [len; culm_len];
        
        culm_da = [0,0,0];
        culm_len = 0.0;
        cnt = 0;
    end
end
culm_da = culm_da + tlgroup_da(i+1, :);
culm_len = culm_len + tlgroup_length(i);
cnt = cnt + 1;
triple_line = [triple_line; tlgroup_id(i+1,:)];
dihedral_angle = [dihedral_angle; culm_da/cnt];
len = [len; culm_len];

tl_info = [dihedral_angle, len];

% triple_line_unique = unique(tlgroup_id, 'rows');
end


% %% ######################################### Check ######################################### 
% % """
% % Idea: plot a group of triangles to see if things are correct.
% % Helper variable: fl_idx_, with this data, TL triangle groups can be
% % visualized together with grains. 
% % Or, just use the following code to plot the nodes. 
% % """
% %% ##### Find bad groups #####
% notfull = abs(sum(dihedral_angle, 2) - 360) > 0.1;
% idx_notfull = (1:size(dihedral_angle,1))';
% idx_notfull = idx_notfull(notfull);
% print('an4_TLWeightedDihedralAnalge')
% tmp = length(idx_notfull)/length(dihedral_angle);
% disp([num2str(tmp), '   triple lines bad'])
% 
% %% ##### Display info of bad groups #####
% group = idx_notfull(randi(length(idx_notfull)));
% start = (group-1)*3 + 1;
% final = (group-1)*3 + 3;
% fl_group = tl_tri_fl(start:final, :)
% fl_idx_group = tl_tri_fl_idx(start:final)
% disp(['dihedral angles = ', num2str(dihedral_angle(group,:)), ...
%     ';  sum(angs) = ', num2str(sum(dihedral_angle(group,:)))]);
% 
% 
% %% ##### Plot nodes of the group #####
% tmp = unique(node_group(:));
% node_cnt = [tmp,histc(node_group(:),tmp)];
% mask_tl = node_cnt(:,2) > 1;
% coords_tlnode = node_coords(node_cnt(mask_tl, 1), :);
% coords_fnode = node_coords(node_cnt(~mask_tl, 1), :);
% 
% scatter3(coords_fnode(:,1), coords_fnode(:,2), coords_fnode(:,3), 'filled')
% hold on
% scatter3(coords_tlnode(:,1), coords_tlnode(:,2), coords_tlnode(:,3), 'filled')
% 
% daspect([1,1,1])
% rotate3d on






