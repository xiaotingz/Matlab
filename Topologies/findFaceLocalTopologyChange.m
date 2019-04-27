function [avg_nng_diff, max_nng_diff] = findFaceLocalTopologyChange(file_1, file_2, faces, look_up_table, triple_line)
% ############################################################################
% * Inputs 
%     - faces = [n, 2]
%         Should be of the PREVIOUS STATE. Returned by trackUniqueFace.m in /Grain Tracking
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
% file_1 = file_an4;
% file_2 = file_an5;
% faces = faces_an4;
% triple_line = triple_line_an4;
% ---------------------------------------------------------------

num_neigh_1 = h5read(file_1, '/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')';
num_neigh_2 = h5read(file_2, '/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')';
surf_grain_1 = h5read(file_1, '/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures')';
surf_grain_2 = h5read(file_2, '/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures')';
num_neigh_1(1) = [];
num_neigh_2(1) = [];
surf_grain_1(1) = [];
surf_grain_2(1) = [];

% ##### Prepare num neigh difference (nn_diff) of grains #####
% """
% - It's assumed that if a grain can't be found in the next state, its
% nn_diff = - num_neigh_initialstate.
% - Sign of nn_diff: + = grow, - = shrink
% """
look_up_table = sortrows(look_up_table, 1);
nn_diff = - num_neigh_1;
nn_diff(look_up_table(:,1)) = (num_neigh_2(look_up_table(:,2)) - num_neigh_1(look_up_table(:,1)));

avg_nng_diff = zeros(size(faces, 1), 1);
max_nng_diff = zeros(size(faces, 1), 1);

for i = 1:size(faces, 1)
    % ##### Get Connected Grains Of This Face #####
    mask_TL = (sum(ismember(triple_line, faces(i, :)), 2) == 2);
    connected_grains = unique(triple_line(mask_TL, :));
    % ----- If No Good TL On Face, Just Use The Two Label Grains -----
    if sum(connected_grains) == 0
        connected_grains = faces(i, :);
    end
    % ##### (Num Neighbors) of the Connected Grains #####
    nn_diff_connected_grains = nn_diff(connected_grains);
    avg_nng_diff(i) = sum(nn_diff_connected_grains)/length(nn_diff_connected_grains);
    [~, idx] = max(abs(nn_diff_connected_grains));
    max_nng_diff(i) = nn_diff_connected_grains(idx);

end    

% """
% next step: filter out the grains due to field of view change 
% """

end



