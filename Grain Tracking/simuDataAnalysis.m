% file_1 = '/Users/xiaotingzhong/Desktop/Datas/Ni_simu/simu_An4__seg.dream3d';
% % file_2 = '/Users/xiaotingzhong/Desktop/Datas/Ni_simu/simu_An5_clean_seg.dream3d';
% 
% grain_id_1 = squeeze(h5read(file_1, '/DataContainers/ImageDataContainer/CellData/FeatureIds'));
% % grain_id_2 = squeeze(h5read(file_2, '/DataContainers/ImageDataContainer/CellData/FeatureIds'));


% tmp1 = P_an4(:,:,1);
% tmp2 = P_an5(:,:,1);
% 
% tmp = (P_an4 == P_an5);
% sum(tmp(:)) / (size(P_an4, 1) * size(P_an4, 2) * size(P_an4, 3)) * 100



% Outline:
%   1. From simulation data, assign new grainIDs
%           Ignore new grains whose size < min_size_thres.
%   2. Tracked the new_ID between an4 and an5. 
%           Try use centroid position for differentiating the pieces.
%   3. Track between new_ID & D3D_ID


test = containers.Map('KeyType', 'double', 'ValueType', 'any');
test(1.0) = [2.0, 3.0];
test(2.0) = 1;






