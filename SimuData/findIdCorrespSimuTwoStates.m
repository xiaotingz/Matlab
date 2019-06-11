function [corresp_simu_new_45, id_new_an4, id_new_an5] = findIdCorrespSimuTwoStates...
            (id_new_an4, id_new_an5, id_corresp_an4, id_corresp_an5, unique_id_given_an5)
% ##################################################################
% * Notes
%     - After re-id grains in the simulation data, this script finds the
%     correspondence between id_news. 
%     - NOTE this algorithm is not good. It finds corresp for whichever
%     state that have fewer grains but this is not good. the magnitude of
%     score should be considered and pairs with low scores should be
%     ignored.
%     - checked with simu_An4_new_id.dream3d and simu_An5_new_id.dream3d,
%     most corresps should be good.
% * Input
%     - id_new_ = [n1, n2, n3]
%           returned by findNewId.m
%     - id_corresp = {m, 1}
%           returned findIdCorrespSameStateOldNew.m
%     - unique_id_given_an5 = unique(id_given_an5);
% * Output
%     - corresp_simu_new_45 = [grain_id_an4, grain_id_an5]
%         Valid corresp between the two simulations, values from id_new_
%     - id_new_ = [n1, n2, n3]
%         cleaned by size threshold, smaller than 27 are assigned as -1.
%         see findSizesAndCentroids.m
% ##################################################################
% % ----------------------- load debug data -----------------------
% 
% % ---------------------------------------------------------------

% % ############################# Find size and centroid for id_new #############################
[sizes_an4, centroids_an4, id_new_an4] = findSizesAndCentroids(id_new_an4);
[sizes_an5, centroids_an5, id_new_an5] = findSizesAndCentroids(id_new_an5);

% ############################# Track id_new between two simulation states #############################
good_id_new_an4 = unique(id_new_an4);
good_id_new_an5 = unique(id_new_an5);
good_id_new_an4(1:2) = [];
good_id_new_an5(1:2) = [];
% """ id_given = 1 is the out large periphery """
unique_id_given_an5(1) = [];
corresp_simu_new_45 = containers.Map(good_id_new_an4, - ones(size(good_id_new_an4)));


for i = 2:length(unique_id_given_an5)
    % """ start from an5 because some of an4 grains may have disappeared """
    id_given = unique_id_given_an5(i);
    grains_new_an4 = id_corresp_an4(id_given);
    grains_new_an5 = id_corresp_an5(id_given);
    grains_new_an4 = grains_new_an4(ismember(grains_new_an4, good_id_new_an4));
    grains_new_an5 = grains_new_an5(ismember(grains_new_an5, good_id_new_an5));
    
    % ----- multiple islands in either state -----
    if ~ isempty(grains_new_an4) && ~ isempty(grains_new_an5)
        if length(grains_new_an4) > 1 || length(grains_new_an5) > 1
            % """ work from an4 to an5, start from the biggest grains """
            if length(grains_new_an5) > length(grains_new_an4)
                data_an4 = [grains_new_an4', sizes_an4(grains_new_an4), centroids_an4(grains_new_an4, :)];
                data_an4 = sortrows(data_an4, 2);
                data_an5 = [grains_new_an5', sizes_an5(grains_new_an5), centroids_an5(grains_new_an5, :)];
                while size(data_an4, 1) > 0
                    % ----- compute score -----
                    size_diffs = data_an5(:,2) - data_an4(1, 2);
                    size_diffs = size_diffs ./ max(size_diffs);
                    centroid_diffs = data_an5(:,3:5) - repmat(data_an4(1, 3:5), size(data_an5, 1), 1);
                    centroid_diffs = vecnorm(centroid_diffs, 2, 2);
                    centroid_diffs = centroid_diffs ./ max(centroid_diffs);
                    score = size_diffs + centroid_diffs;
                    % ----- record corresp -----
                    [~, idx] = min(score);
                    corresp_simu_new_45(data_an4(1,1)) = data_an5(idx, 1);
                    data_an4(1, :) = [];
                    data_an5(idx, :) = [];
                end
            % """ work from an5 to an4, start from the biggest grains """
            else
                data_an4 = [grains_new_an4', sizes_an4(grains_new_an4), centroids_an4(grains_new_an4, :)];
                data_an5 = [grains_new_an5', sizes_an5(grains_new_an5), centroids_an5(grains_new_an5, :)];
                data_an5 = sortrows(data_an5, 2);
                while size(data_an5, 1) > 0
                    % ----- compute score -----
                    size_diffs = data_an4(:,2) - data_an5(1, 2);
                    size_diffs = size_diffs ./ max(size_diffs);
                    centroid_diffs = data_an4(:,3:5) - repmat(data_an5(1, 3:5), size(data_an4, 1), 1);
                    centroid_diffs = vecnorm(centroid_diffs, 2, 2);
                    centroid_diffs = centroid_diffs ./ max(centroid_diffs);
                    score = size_diffs + centroid_diffs;
                    % ----- record corresp -----
                    [~, idx] = min(score);
                    corresp_simu_new_45(data_an4(idx, 1)) = data_an5(1,1);
                    data_an5(1, :) = [];
                    data_an4(idx, :) = [];
                end
            end

        % ----- one island in both states -----
        else
            corresp_simu_new_45(grains_new_an4) = grains_new_an5;
        end
    end

end
corresp_simu_new_45 = [cell2mat(corresp_simu_new_45.keys)', cell2mat(corresp_simu_new_45.values)'];
corresp_simu_new_45 = corresp_simu_new_45(all(corresp_simu_new_45>0, 2), :);

end