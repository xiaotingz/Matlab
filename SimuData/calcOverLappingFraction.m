function ret = calcOverLappingFraction(file_1, file_2, look_up_table_tmp)
% ##########################################################################
% * Notes
%     - This script is to prepare the geometric and topological features. 
%     - Preparation of the topology features: see main.m in /Topologies
% ##########################################################################
% ----------------------- load debug data -----------------------
% file_1 = '/Users/xiaotingzhong/Desktop/Datas/Ni_simu/simu_An4__seg.dream3d';
% file_2 = '/Users/xiaotingzhong/Desktop/Datas/Ni_simu/simu_An5_clean_seg.dream3d';

% file_1 = '/Users/xiaotingzhong/Desktop/Datas/Ni_simu/simu_An4_clean_seg.dream3d';
% file_2 = '/Users/xiaotingzhong/Desktop/Datas/Ni_simu/simu_An5_clean_seg.dream3d';
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_simu/d3d_simu_corresp.mat', 'corresp_d3d_45')
% look_up_table_tmp = corresp_d3d_45;
% ---------------------------------------------------------------

grain_id_1 = squeeze(h5read(file_1, '/DataContainers/ImageDataContainer/CellData/FeatureIds'));
grain_id_2 = squeeze(h5read(file_2, '/DataContainers/ImageDataContainer/CellData/FeatureIds'));
size_1 = size(grain_id_1, 1);
size_2 = size(grain_id_1, 2);
size_3 = size(grain_id_1, 3);

if narg == 3
    % ----- sort look_up_table -----
    num_unique_grains_1 = length(unique(grain_id_1));
    % """ note for grain 0 in d3d files"""
    if ismember(0, grain_id_1)
        look_up_table = [(0:num_unique_grains_1 - 1)', - ones(num_unique_grains_1, 1)];
        look_up_table(look_up_table_tmp(:,1) + 1, 2) = look_up_table_tmp(:,2);
    else
        look_up_table = [(1:num_unique_grains_1)', - ones(num_unique_grains_1, 1)];
        look_up_table(look_up_table_tmp(:,1), 2) = look_up_table_tmp(:,2);
    end
    look_up_table = containers.Map(look_up_table(:,1), look_up_table(:,2));
    
    for i = 1:size_1
        for j = 1:size_2
            for k = 1:size_3
                grain_id_1(i, j, k) = look_up_table(grain_id_1(i, j, k));
            end
        end
    end
end

% """ note grain 1 is the bad one in Jade's simulation """
tmp = (grain_id_1 == grain_id_2 & grain_id_1 > 0);
ret = sum(tmp(:)) / (size_1 * size_2 * size_3) * 100;

end





