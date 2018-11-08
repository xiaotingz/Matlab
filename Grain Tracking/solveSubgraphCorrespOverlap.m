function subgraph_corresp = solveSubgraphCorrespOverlap(subgraph_id_an4,subgraph_id_an5, ...
facetri_normal_an4, facetri_normal_an5, facetri_nodeid_an4, facetri_nodeid_an5, facenode_id_an4, facenode_id_an5, ...
facenode_coord_an4, facenode_coord_an5, facetri_area_an4, facetri_area_an5)
% ############################################################################
% * Input 
%     - see main_TrackNodes.m
%  * Output
%     - subgraph_corresp = [n, 3]
%         The correspondence between subgraph_id, corresp_score
%  * Notes
%     - Used in main_TrackNodes.m
%     - Dependency: project3DPointsTo2DPlane.m
% ############################################################################
% % ------------------ load data for debug --------------------
% %-----------------------------------------------------------
%%
% ----- threshold for determining good corresp_pair and same_state_overlap -----
thres_dup = 0.2;
thres_keepN = 0.5;

% ######################################## Section I, One Face Has Only One Piece ########################################
if max(subgraph_id_an4) == 1 || max(subgraph_id_an5) == 1
    % ##### determine the basis plane and the candidate planes #####
    if max(subgraph_id_an4) == 1 && max(subgraph_id_an5) > 1
        % ----- convert node subgraph info to triangles -----
        map_an4 = containers.Map(facenode_id_an5, subgraph_id_an5);
        mask_tri_subgraph_an5 = zeros(length(facetri_nodeid_an5), 1);
        for i = 1:length(facetri_nodeid_an5)
            % -- because all nodes of a triangle should belong to the same subgraph, no need to loop them all --
            mask_tri_subgraph_an5(i) = map_an4(facetri_nodeid_an5(i));
        end

        % ----- write infomation information to basis plane and candidate plane -----        
        node_coord_bp = facenode_coord_an4;
        tri_normal_bp = facetri_normal_an4;
        normal_avg_bp = sum(tri_normal_bp);
        normal_avg_bp = normal_avg_bp/norm(normal_avg_bp);
        area_bp = sum(facetri_area_an4);
        n_sg_cp = max(subgraph_id_an5);
        node_coord_cp = cell(n_sg_cp, 1);
        normal_avg_cp = zeros(n_sg_cp, 3);
        area_cp = zeros(n_sg_cp, 1);
        for i = 1:n_sg_cp
            node_coord_cp{i} = facenode_coord_an5(subgraph_id_an5 == i, :);
            tri_normal_cp_i = facetri_normal_an5(mask_tri_subgraph_an5==i, :);
            normal_avg_cp(i,:) = sum(tri_normal_cp_i)/length(tri_normal_cp_i);
            normal_avg_cp(i,:) = normal_avg_cp(i,:)/norm(normal_avg_cp(i,:));
            area_cp(i) = sum(facetri_area_an5(mask_tri_subgraph_an5==i, :));
        end

        % ----- record the basis plane is in which state ----- 
        basis_state = 1;

    elseif max(subgraph_id_an4) > 1 && max(subgraph_id_an5) == 1
        % ----- convert node subgraph info to triangles -----
        map_an4 = containers.Map(facenode_id_an4, subgraph_id_an4);
        mask_tri_subgraph_an4 = zeros(length(facetri_nodeid_an4), 1);
        for i = 1:length(facetri_nodeid_an4)
            % -- because all nodes of a triangle should belong to the same subgraph, no need to loop them all --
            mask_tri_subgraph_an4(i) = map_an4(facetri_nodeid_an4(i));
        end

        % ----- write infomation information to basis plane and candidate plane -----        
        node_coord_bp = facenode_coord_an5;
        tri_normal_bp = facetri_normal_an5;
        normal_avg_bp = sum(tri_normal_bp);
        normal_avg_bp = normal_avg_bp/norm(normal_avg_bp);
        area_bp = sum(facetri_area_an5);
        n_sg_cp = max(subgraph_id_an4);
        node_coord_cp = cell(n_sg_cp, 1);
        normal_avg_cp = zeros(n_sg_cp, 3);
        area_cp = zeros(n_sg_cp, 1);
        for i = 1:n_sg_cp
            node_coord_cp{i} = facenode_coord_an4(subgraph_id_an4 == i, :);
            tri_normal_cp_i = facetri_normal_an4(mask_tri_subgraph_an4==i, :);
            normal_avg_cp(i,:) = sum(tri_normal_cp_i)/length(tri_normal_cp_i);
            normal_avg_cp(i,:) =  normal_avg_cp(i,:)/norm(normal_avg_cp(i,:));
            area_cp(i) = sum(facetri_area_an4(mask_tri_subgraph_an4==i, :));
        end

        % ----- record the basis plane is in which state ----- 
        basis_state = 2;

    end
    % clearvars -except node_id_bp node_id_cp1 node_id_cp2 tri_normal_bp tri_normal_cp1 tri_normal_cp2

    % ----- basis plane normal -----
    normal_bp = normal_avg_bp;

    % ----- basis plane origion from centroid of basis plane -----
    centroid_bp = sum(node_coord_bp)/size(node_coord_bp, 1);


    % ##### Projection Score #####
    % ----- project nodes onto the basis plane -----
    [node_proj2d_bp, native2d_x] = project3DPointsTo2DPlane(node_coord_bp, centroid_bp, normal_bp, 'bp');
    node_proj2d_cp = cell(n_sg_cp, 1);
    for i = 1:n_sg_cp
        [node_proj2d_cp{i}, ~] = project3DPointsTo2DPlane(node_coord_cp{i}, centroid_bp, normal_bp, native2d_x);
    end

    % ----- polygon from boundary of the projected nodes and the intersections -----
    bound_bp = boundary(node_proj2d_bp(:,1), node_proj2d_bp(:,2));
    pgon_bp = polyshape(node_proj2d_bp(bound_bp, 1), node_proj2d_bp(bound_bp, 2));
    pgon_cp = cell(n_sg_cp, 1);
    score_overlap = zeros(n_sg_cp, 1);
    for i = 1:n_sg_cp
        bound_cp_i = boundary(node_proj2d_cp{i}(:,1), node_proj2d_cp{i}(:,2));
        pgon_cp{i} = polyshape(node_proj2d_cp{i}(bound_cp_i, 1), node_proj2d_cp{i}(bound_cp_i, 2));
        % """
        % Either (cp covered by bp) | (bp covered by cp) should be considered as perfect overlap
        % """
        denorminator = min(area(pgon_bp), area(pgon_cp{i}));
        score_overlap(i) = area(intersect(pgon_bp,pgon_cp{i}))/denorminator;
    end

    % ----- normal direction difference -----
    % """
    % - Notice switching symmetry, max(normal_diff)=90°
    % - One can also do a Baysian estimation for score_norm can calculate it from
    % KL divergence. However, probably not necessary. 
    % """
    inter_angle = zeros(n_sg_cp, 1);
    for i = 1:n_sg_cp
        inter_angle(i) = abs(atand(norm(cross(normal_avg_bp, normal_avg_cp(i,:)))/dot(normal_avg_bp, normal_avg_cp(i,:))));
    end
    score_norm = (1 - inter_angle/45);

    % ----- total score -----
    % """
    % - The logic is that 1. overlapping area should be large 2. normal should be consistant. 
    % """
    score = (0.5 * score_overlap + 0.5 * score_norm);

    
    % ##### Subgraph Overlap Within cp #####
    % ----- check if subgraphs of cp are overlapping themselves -----
    % """
    % If yes, chose the one as the larger and closer one. 
    % """
    % dup_score = area(intersect(pgon_cp{1},pgon_cp{2}))/min(area(pgon_cp{1}), area(pgon_cp{2}));
    cp_dup_checklist = combnk(1:n_sg_cp, 2);
    cp_dup_score = zeros(size(cp_dup_checklist, 1), 1);
    for i = 1:size(cp_dup_score, 1)
        idx_an4 = cp_dup_checklist(i, 1);
        idx_an5 = cp_dup_checklist(i, 2);
        cp_dup_score(i) = area(intersect(pgon_cp{idx_an4},pgon_cp{idx_an5}))/min(area(pgon_cp{idx_an4}), area(pgon_cp{idx_an5}));
    end

    mask_keep = score > thres_keepN;
    array_keep = (1:n_sg_cp)';
    % ----- no overlap -----
    if all(cp_dup_score < thres_dup)
        idx_keep = array_keep(mask_keep);
        N = sum(mask_keep);
        tmp = ones(N, 1);
        % ----- if basis plane is in the first state -----
        if basis_state == 1
            subgraph_corresp = [tmp, idx_keep];
        % ----- if basis plane is in the second state -----
        else
            subgraph_corresp = [idx_keep, tmp];
        end

    % ----- some overlap -----
    else
        % ----- check each overlapping pair, remove the face with smaller score -----
        % """
        % - The score is now weighted by relative area
        % - Designed as neither large expansion nor serious shrinkage is likely. 
        % Expand to twice as large have the same possibility as shrink to 1/2 the original area. 
        % """
        area_relative = min([area_cp./area_bp, area_bp./area_cp], [], 2);
        weight_areachange = gaussmf(area_relative,[0.5, 1]);
        score_relative = score.*weight_areachange;
        mask_overlap_pair = (cp_dup_score > thres_dup);
        idx_overlap_pair = cp_dup_checklist(mask_overlap_pair, :);
%         score_overlap_pair = score_relative(idx_overlap_pair);
%         for i = 1:size(idx_overlap_pair)
%             [~, idx_remove] = min(score_overlap_pair);
%             mask_keep(idx_overlap_pair(idx_remove)) = false;
%         end
        [~, idx_remove] = min(score_relative(idx_overlap_pair), [], 2);
        idx_remove = idx_overlap_pair(:, idx_remove);
        mask_keep(idx_remove) = false;
        
        idx_keep = array_keep(mask_keep);
        N = sum(mask_keep);
        tmp = ones(N, 1);
        % ----- if basis plane is in the first state -----
        if basis_state == 1
            subgraph_corresp = [tmp, idx_keep];
        % ----- if basis plane is in the second state -----
        else
            subgraph_corresp = [idx_keep, tmp];
        end
    end

    % ##### Final Check, Make Sure There Is At Least One Corresp #####
    % ----- no good fitting -----
    if sum(mask_keep) == 0
        [~, idx_keep] = max(score);
        % ----- if basis plane is in the first state -----
        if basis_state == 1
            subgraph_corresp = [1, idx_keep];
        % ----- if basis plane is in the second state -----
        else
            subgraph_corresp = [idx_keep, 1];
        end
    end


% ######################################## Section II, Both Faces Have Multiple Pieces ########################################
else
    % ##### Determine A Basis State by The Largest Single Piece #####
    % ----- convert node subgraph info to triangles -----
    map_an4 = containers.Map(facenode_id_an4, subgraph_id_an4);
    mask_tri_subgraph_an4 = zeros(length(facetri_nodeid_an4), 1);
    for i = 1:length(facetri_nodeid_an4)
        mask_tri_subgraph_an4(i) = map_an4(facetri_nodeid_an4(i));
    end
    map_an5 = containers.Map(facenode_id_an5, subgraph_id_an5);
    mask_tri_subgraph_an5 = zeros(length(facetri_nodeid_an5), 1);
    for i = 1:length(facetri_nodeid_an5)
        mask_tri_subgraph_an5(i) = map_an5(facetri_nodeid_an5(i));
    end
    % ----- calculate area of each piece (subgraph) -----
    n_sg_an4 = max(subgraph_id_an4);
    n_sg_an5 = max(subgraph_id_an5);
    area_an4 = zeros(n_sg_an4, 1);
    area_an5 = zeros(n_sg_an5, 1);
    for i = 1:n_sg_an4
        area_an4(i) = sum(facetri_area_an4(mask_tri_subgraph_an4==i));
    end
    for i = 1:n_sg_an5
        area_an5(i) = sum(facetri_area_an5(mask_tri_subgraph_an5==i));
    end
    % ----- determine a basis state -----
    if max(area_an4) > max(area_an5)
        basis_state = 1;
    else
        basis_state = 2;
    end
    
    % ##### Repeat Section I for Every Piece of The Basis State #####
    if basis_state == 1
        n_sg_bp = n_sg_an4;
        n_sg_cp = n_sg_an5;
        area_bp = area_an4;
        area_cp = area_an5;
        % ----- collect data for basis planes -----
        node_coord_bp = cell(n_sg_bp, 1);
        normal_avg_bp = zeros(n_sg_bp, 3);
        centroid_bp = zeros(n_sg_bp, 3);
        for i = 1:n_sg_bp
            node_coord_bp{i} = facenode_coord_an4(subgraph_id_an4 == i, :);
            tri_normal_bp_i = facetri_normal_an4(mask_tri_subgraph_an4==i, :);
            normal_avg_bp(i, :) = sum(tri_normal_bp_i)/length(tri_normal_bp_i);
            centroid_bp(i, :) = sum(node_coord_bp{i})/size(node_coord_bp{i}, 1);
        end
        % ----- collect data for candidate planes -----
        node_coord_cp = cell(n_sg_cp, 1);
        normal_avg_cp = zeros(n_sg_cp, 3);
        for i = 1:n_sg_cp
            node_coord_cp{i} = facenode_coord_an5(subgraph_id_an5 == i, :);
            tri_normal_cp_i = facetri_normal_an5(mask_tri_subgraph_an5==i, :);
            normal_avg_cp(i, :) = sum(tri_normal_cp_i)/length(tri_normal_cp_i);
        end 
        
    else
        n_sg_bp = n_sg_an5;
        n_sg_cp = n_sg_an4;
        area_bp = area_an5;
        area_cp = area_an4;
        % ----- collect data for basis planes -----
        node_coord_bp = cell(n_sg_bp, 1);
        normal_avg_bp = zeros(n_sg_bp, 3);
        centroid_bp = zeros(n_sg_bp, 3);
        for i = 1:n_sg_bp
            node_coord_bp{i} = facenode_coord_an5(subgraph_id_an5 == i, :);
            tri_normal_bp_i = facetri_normal_an5(mask_tri_subgraph_an5==i, :);
            normal_avg_bp(i, :) = sum(tri_normal_bp_i)/length(tri_normal_bp_i);
            centroid_bp(i, :) = sum(node_coord_bp{i})/size(node_coord_bp{i}, 1);
        end
        % ----- collect data for candidate planes -----
        node_coord_cp = cell(n_sg_cp, 1);
        normal_avg_cp = zeros(n_sg_cp, 3);
        for i = 1:n_sg_cp
            node_coord_cp{i} = facenode_coord_an4(subgraph_id_an4 == i, :);
            tri_normal_cp_i = facetri_normal_an4(mask_tri_subgraph_an4==i, :);
            normal_avg_cp(i, :) = sum(tri_normal_cp_i)/length(tri_normal_cp_i);
        end
    end
    
    % ##### Projection Score #####
    % ----- initialization -----
    score_overlap = zeros(n_sg_cp, n_sg_bp);
    score_norm = zeros(n_sg_cp, n_sg_bp);
    score = zeros(n_sg_cp, n_sg_bp);
    subgraph_corresp = [];
    area_cp_mat = repmat(area_cp, 1, n_sg_bp);
    area_bp_mat = repmat(area_bp', n_sg_cp, 1);
    area_relative = reshape([area_cp_mat./area_bp_mat, area_bp_mat./area_cp_mat], n_sg_cp, n_sg_bp, 2);
    area_relative = min(area_relative, [], 3);
    for idx_bp = 1:n_sg_bp
        % ----- project nodes onto the basis plane -----
        [node_proj2d_bp, native2d_x] = project3DPointsTo2DPlane(node_coord_bp{idx_bp}, ...
            centroid_bp(idx_bp, :), normal_avg_bp(idx_bp, :), 'bp');
        node_proj2d_cp = cell(n_sg_cp, 1);
        for i = 1:n_sg_cp
            [node_proj2d_cp{i}, ~] = project3DPointsTo2DPlane(node_coord_cp{i}, centroid_bp(idx_bp, :), normal_avg_bp(idx_bp, :), native2d_x);
        end
        
        % ----- polygon from boundary of the projected nodes and the intersections -----
        bound_bp = boundary(node_proj2d_bp(:,1), node_proj2d_bp(:,2));
        pgon_bp = polyshape(node_proj2d_bp(bound_bp, 1), node_proj2d_bp(bound_bp, 2));
        pgon_cp = cell(n_sg_cp, 1);
        for i = 1:n_sg_cp
            bound_cp_i = boundary(node_proj2d_cp{i}(:,1), node_proj2d_cp{i}(:,2));
            pgon_cp{i} = polyshape(node_proj2d_cp{i}(bound_cp_i, 1), node_proj2d_cp{i}(bound_cp_i, 2));
            denorminator = min(area(pgon_bp), area(pgon_cp{i}));
            score_overlap(i, idx_bp) = area(intersect(pgon_bp,pgon_cp{i}))/denorminator;
        end
        
        % ----- normal direction difference -----
        inter_angle = zeros(n_sg_cp, 1);
        for i = 1:n_sg_cp
            inter_angle(i) = abs(atand(norm(cross(normal_avg_bp(idx_bp, :), normal_avg_cp(i,:)))/dot(normal_avg_bp(idx_bp, :), ...
                normal_avg_cp(i,:))));
        end
        score_norm(:, idx_bp) = (1 - inter_angle/45);

        % ----- total score -----
        score(:, idx_bp) = (0.5 * score_overlap(:, idx_bp) + 0.5 * score_norm(:, idx_bp));
        
        % ----- check if subgraphs of cp are overlapping themselves -----
        cp_dup_checklist = combnk(1:n_sg_cp, 2);
        cp_dup_score = zeros(size(cp_dup_checklist, 1), 1);
        for i = 1:size(cp_dup_score, 1)
            idx_an4 = cp_dup_checklist(i, 1);
            idx_an5 = cp_dup_checklist(i, 2);
            cp_dup_score(i) = area(intersect(pgon_cp{idx_an4},pgon_cp{idx_an5}))/min(area(pgon_cp{idx_an4}), area(pgon_cp{idx_an5}));
        end
        
        % ----- if overlap, adjust the scores to assure keep only one corresp -----
        weight_areachange = gaussmf(area_relative(:, idx_bp),[0.5, 1]);
        score_relative = score(:, idx_bp).*weight_areachange;
        mask_overlap_pair = (cp_dup_score > thres_dup);
        idx_overlap_pair = cp_dup_checklist(mask_overlap_pair, :);
        [~, idx_remove] = min(score_relative(idx_overlap_pair), [], 2);
        idx_remove = idx_overlap_pair(:, idx_remove);
        score(idx_remove, idx_bp) = 0;
    end
    
    % ##### Write Corresp From score #####
    mask_corresp = (score > thres_keepN);
    % ----- if no good corresp, make sure at least the largest piece have a corresp -----
    if sum(mask_corresp) == 0
        [~, corresp_bp] = max(area_bp);
        [~, corresp_cp] = max(score(:, corresp_bp));
    % ----- if correps have been found, just use -----
    else
        [corresp_cp, corresp_bp] = find(mask_corresp==1);
    end
    
    if basis_state == 1
        subgraph_corresp = [corresp_bp, corresp_cp];
    else
        subgraph_corresp = [corresp_cp, corresp_bp];
    end
end     

end

% ################################### Visual Check: 3D and projected 2D grain face ###################################
% disp(['idx = ', num2str(idx)]);
% info_score = 'score = ';
% info_score_norm = 'score_norm = ';
% for i = 1:n_sg_cp
%     info_score = [info_score, num2str(score(i)), '  '];
%     info_score_norm = [info_score_norm, num2str(score_norm(i)), '  '];
% end
% info_dup = num2str([cp_dup_checklist, cp_dup_score]);
% disp(info_score);
% disp(info_score_norm);
% disp('dupplication sores:')
% disp(info_dup);
% disp('subgraph corresp:')
% disp(num2str(subgraph_corresp));
% disp(' ');
% 
% colors = get(gca,'colororder');
% x_to_y = X_to_Y{idx};
% face_node_info = getSingleFaceNodes(tracked_uniqueface_an4(idx,:), tracked_uniqueface_an5(idx,:));
% visualizeFace(face_node_info, x_to_y)
% hold on
% coord = facenode_coord_an5;
% id = subgraph_id_an5;
% scatter3(coord(id==1, 1), coord(id==1, 2), coord(id==1, 3), 'MarkerFaceColor', 'r')


% figure
% plot(pgon_bp, 'FaceColor', colors(1, :))
% hold on
% for i = 1:n_sg_cp
%     plot(pgon_cp{i}, 'FaceColor', colors(2, :))
% end
% daspect([1 1 1])

% %% ----- visual check: the 2d projection local coordinates consistent -----
% figure
% [test1, ~] = project3DPointsTo2DPlane(facenode_coord_an5, centroid_bp, normal_bp, native2d_x);
% scatter(test1(:,1), test1(:,2), 'MarkerFaceColor', colors(1, :))
% hold on
% test2 = project3DPointsTo2DPlane(node_coord_cp{1}, centroid_bp, normal_bp, native2d_x);
% test3 = project3DPointsTo2DPlane(node_coord_cp{2}, centroid_bp, normal_bp, native2d_x);
% scatter(test2(:,1), test2(:,2), 'MarkerEdgeColor', colors(2, :))
% scatter(test3(:,1), test3(:,2), 'MarkerEdgeColor', colors(2, :))
% plot(pgon_bp, 'FaceColor', colors(1, :))
% plot(pgon_cp{1}, 'FaceColor', colors(2, :))
% plot(pgon_cp{2}, 'FaceColor', colors(2, :))

