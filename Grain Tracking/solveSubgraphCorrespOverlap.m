% function subgraph_corresp = solveSubgraphCorrespOverlap(file_1, file_2, X_to_Y, look_up_table, face_piecewise)
% ############################################################################
% * Input 
%     - normal_ = [3, 1], average triangle normal of the faces
%         The underlying logic is that a grain face's orientation change
%         can be quantified by the difference in average triangle normal,
%         which is not necessarily true actually.
%     - nodes_ = [n, 3], coordinates of mesh nodes
%     - 
%  * Output
%     - 
%  * Notes
%     - Prototype in subgraphCorrespAnalysis: Project Nodes Onto Basis Plane, 
% ############################################################################
% % ------------------ load data for debug --------------------
load('180927_tmp.mat', 'X_to_Y', 'file_an4', 'file_an5', 'onestate_multipiece', 'look_up_table')
file_1 = file_an4;
file_2 = file_an5;
clear file_an4 file_an5
% %-----------------------------------------------------------
% ##### Data Preparation #####
facelabel_1 = double(h5read(file_1,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
facelabel_2 = double(h5read(file_2,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
tri_normal_1 = double(h5read(file_1,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
tri_normal_2 = double(h5read(file_2,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals'))';
    % -- NOTE triNodes are indexes starting from zero --
tri_node_1 = 1 + double(h5read(file_1,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
tri_node_2 = 1 + double(h5read(file_2,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coord_1 = double(h5read(file_1,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
node_coord_2 = double(h5read(file_2,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
% ---- filter bad data -----
mask_1 = all(facelabel_1 > 0, 2);
facelabel_1 = facelabel_1(mask_1, :);
tri_node_1 = tri_node_1(mask_1, :);
tri_normal_1 = tri_normal_1(mask_1, :);
mask_2 = all(facelabel_2 > 0, 2);
facelabel_2 = facelabel_2(mask_2, :);
tri_node_2 = tri_node_2(mask_2, :);
tri_normal_2 = tri_normal_2(mask_2, :);
color1 = [0, 0.4470, 0.7410];
color2 = [0.9290, 0.6940, 0.1250];

[tracked_uniqueface_1, tracked_uniqueface_2] = trackUniqueFace(file_1, file_2, look_up_table, 'use_complete_faces');
% ----- find piecewise faces -----
% face_piecewise = [];
% num_pieces = [];
% for i = 1:length(tracked_uniqueface_1)
%     obj_facelabel_1 = tracked_uniqueface_1(i, :);
%     obj_facelabel_2 = tracked_uniqueface_2(i, :);
%     mask_1_sub1 = (facelabel_1(:,1) == obj_facelabel_1(1) & facelabel_1(:,2) == obj_facelabel_1(2));
%     mask_1_sub2 = (facelabel_1(:,1) == obj_facelabel_1(2) & facelabel_1(:,2) == obj_facelabel_1(1));
%     mask_1 = (mask_1_sub1 | mask_1_sub2);
%     mask_2_sub1 = (facelabel_2(:,1) == obj_facelabel_2(1) & facelabel_2(:,2) == obj_facelabel_2(2));
%     mask_2_sub2 = (facelabel_2(:,1) == obj_facelabel_2(2) & facelabel_2(:,2) == obj_facelabel_2(1));
%     mask_2 = (mask_2_sub1 | mask_2_sub2);
%     facenode_id_1 = unique(tri_node_1(mask_1, :));
%     facenode_id_2 = unique(tri_node_2(mask_2, :));
%     facetri_nodeid_1 = tri_node_1(mask_1, :);
%     facetri_nodeid_2 = tri_node_2(mask_2, :);
%     
%     [subgraph_1, subgraph_2] = findSubgraph(facenode_id_1, facenode_id_2, facetri_nodeid_1, facetri_nodeid_2);
%     m = length(unique(subgraph_1));
%     n = length(unique(subgraph_2));
%     if m > 1 || n > 1
%         face_piecewise = [face_piecewise; i];
%         num_pieces = [num_pieces; max(subgraph_an4), max(subgraph_an5)];
%     end
% end


% ##### Corresp of Subgraphs #####
% for i = 1:length(onestate_multipiece)
i = 1;
    % ----- find face node info -----
    idx = onestate_multipiece(i);
    obj_facelabel_1 = tracked_uniqueface_1(idx, :);
    obj_facelabel_2 = tracked_uniqueface_2(idx, :);
    mask_1_sub1 = (facelabel_1(:,1) == obj_facelabel_1(1) & facelabel_1(:,2) == obj_facelabel_1(2));
    mask_1_sub2 = (facelabel_1(:,1) == obj_facelabel_1(2) & facelabel_1(:,2) == obj_facelabel_1(1));
    mask_1 = (mask_1_sub1 | mask_1_sub2);
    mask_2_sub1 = (facelabel_2(:,1) == obj_facelabel_2(1) & facelabel_2(:,2) == obj_facelabel_2(2));
    mask_2_sub2 = (facelabel_2(:,1) == obj_facelabel_2(2) & facelabel_2(:,2) == obj_facelabel_2(1));
    mask_2 = (mask_2_sub1 | mask_2_sub2);
    facenode_id_1 = unique(tri_node_1(mask_1, :));
    facenode_id_2 = unique(tri_node_2(mask_2, :));
    facetri_nodeid_1 = tri_node_1(mask_1, :);
    facetri_nodeid_2 = tri_node_2(mask_2, :);
    facetri_normal_1 = [tri_normal_1(mask_1_sub1, :); -tri_normal_1(mask_1_sub2, :)];
    facetri_normal_2 = [tri_normal_2(mask_2_sub1, :); -tri_normal_2(mask_2_sub2, :)];
    
    [subgraph_1, subgraph_2] = findSubgraph(facenode_id_1, facenode_id_2, facetri_nodeid_1, facetri_nodeid_2);

    % ----- determine the basis plane and the candidate planes -----
    if max(subgraph_1) == 1 && max(subgraph_2) > 1
        % ----- make a dictionary to convert node subgraph info to triangles -----
        map = containers.Map(facenode_id_2, subgraph_2);
        mask_tri_subgraph = zeros(length(facetri_nodeid_2), 1);
        for i = 1:length(facetri_nodeid_2)
            % -- because all nodes of a triangle should belong to the same subgraph, no need to loop them all --
            mask_tri_subgraph(i) = map(facetri_nodeid_2(i));
        end
        
        % ----- write infomation information to basis plane and candidate plane -----        
        node_id_bp = facenode_id_1;
        node_id_cp1 = facenode_id_2(subgraph_2 == 1);
        node_id_cp2 = facenode_id_2(subgraph_2 == 2);
        tri_normal_bp = facetri_normal_1;
        tri_normal_cp1 = facetri_normal_2(mask_tri_subgraph==1, :);
        tri_normal_cp2 = facetri_normal_2(mask_tri_subgraph==2, :);
        
    elseif max(subgraph_1) > 1 && max(subgraph_2) == 1
        % ----- make a dictionary to convert node subgraph info to triangles -----
        map = containers.Map(facenode_id_1, subgraph_1);
        mask_tri_subgraph = zeros(length(facetri_nodeid_1), 1);
        for i = 1:length(facetri_nodeid_1)
            % -- because all nodes of a triangle should belong to the same subgraph, no need to loop them all --
            mask_tri_subgraph(i) = map(facetri_nodeid_1(i));
        end
        
        % ----- write infomation information to basis plane and candidate plane -----        
        node_id_bp = facenode_id_2;
        node_id_cp1 = facenode_id_1(subgraph_1 == 1);
        node_id_cp2 = facenode_id_1(subgraph_1 == 2);
        tri_normal_bp = facetri_normal_2;
        tri_normal_cp1 = facetri_normal_1(mask_tri_subgraph==1, :);
        tri_normal_cp2 = facetri_normal_1(mask_tri_subgraph==2, :);
    end
    clearvars -except node_id_bp node_id_cp1 node_id_cp2 tri_normal_bp tri_normal_cp1 tri_normal_cp2
    
%%    % ----- calc basis plane normal -----
    normal_avg_1 = sum(tri_normal_1(mask_1_sub1, :)) - sum(tri_normal_1(mask_1_sub2, :));
    normal_avg_1 = normal_avg_1/norm(normal_avg_1);
    normal_avg_2 = sum(tri_normal_2(mask_2_sub1, :)) - sum(tri_normal_2(mask_2_sub2, :));
    normal_avg_2 = normal_avg_2/norm(normal_avg_2);
    normal_bp = normal_avg_1;

    % ----- calc basis plane origion from centroid of basis plane -----
    centroid_basis = sum(nodes_basis)/length(nodes_basis);

    [nodes_proj2d_1, native2d_x] = project3DPointsTo2DPlane(nodes_1, centroid_bp, normal_bp, 'bp');

    % ----- polygon and intersection -----
    pgon_1 = polyshape(nodes_proj2d_1(bound_1, 1), nodes_proj2d_1(bound_1, 2));
    pgon_2 = polyshape(nodes_proj2d_2(bound_2, 1), nodes_proj2d_2(bound_2, 2));
    pgon_inter = intersect(pgon_1,pgon_2);

    % ----- normal direction difference -----
    % """
    % Notice switching symmetry, max(normal_diff)=90°
    % """
    norm_diff = abs(atand(norm(cross(normal_avg_1, normal_avg_2))/dot(normal_avg_1, normal_avg_2)));

    % ----- overlap score -----
    % """
    % The logic is that overlapping area and consistent normal are equally
    % important. This logic maybe questionable.
    % """
    score = 0.5 * (area(pgon_inter)/area(pgon_2)) + 0.5 * (1 - norm_diff/90);



% end
