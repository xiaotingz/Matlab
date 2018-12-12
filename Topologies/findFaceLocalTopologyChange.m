function [avg_nndiff, max_nndiff] = findFaceLocalTopologyChange(file_1, file_2, faces, look_up_table)
% ############################################################################
% * Inputs 
%     - faces = [n, 2], returned by trackUniqueFace.m in /Grain Tracking
%     - look_up_table, returned by Aditi's code
% * Outputs
%     - avg_nndiff: average(number of neighbors difference of connected grains)
%     - max_nndiff: max(number of neighbors difference of connected grains)
% * Note
%     - Local topology change is defined as the number of faces change of grains. 
%     - !!! Topology change due to field of view can yield large noise !!!
%     - Dependency: findTripleLines.m
% ############################################################################
% ----------------------- load debug data -----------------------
% file_1 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
% file_2 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d';
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/181107.mat', 'tracked_uniqueface_an4')
% faces = tracked_uniqueface_an4;
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/look_up_table_an4_an5.mat')
% 
% clear tracked_uniqueface_an4 
% ---------------------------------------------------------------

num_neigh_1 = h5read(file_1, '/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')';
num_neigh_2 = h5read(file_2, '/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')';
surf_grain_1 = h5read(file_1, '/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures')';
surf_grain_2 = h5read(file_2, '/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures')';
num_neigh_1(1) = [];
num_neigh_2(1) = [];
surf_grain_1(1) = [];
surf_grain_2(1) = [];

% ##### Find Triple Lines #####
triple_line = findTripleLines(file_1);

% ##### Prepare num neigh difference (nn_diff) of grains #####
% """
% - It's assumed that if a grain can't be found in the next state, its
% nn_diff = - num_neigh_initialstate.
% - Sign of nn_diff: + = grow, - = shrink
% """
look_up_table = sortrows(look_up_table, 1);
nn_diff = - num_neigh_1;
nn_diff(look_up_table(:,1)) = (num_neigh_2(look_up_table(:,2)) - num_neigh_1(look_up_table(:,1)));

avg_nndiff = zeros(size(faces, 1), 1);
max_nndiff = zeros(size(faces, 1), 1);

for i = 1:size(faces, 1)
    % ##### Get Connected Grains Of This Face #####
    mask_TL = (sum(ismember(triple_line, faces(i, :)), 2) == 2);
    connected_grains = unique(triple_line(mask_TL, :));

    % ##### ?(Num Neighbors) of the Connected Grains #####
    nn_diff_connected_grains = nn_diff(connected_grains);
    avg_nndiff(i) = sum(nn_diff_connected_grains)/length(nn_diff_connected_grains);
    [~, idx] = max(abs(nn_diff_connected_grains));
    max_nndiff(i) = nn_diff_connected_grains(idx);
end    

% """
% next step: filter out the grains due to field of view change 
% """

end



