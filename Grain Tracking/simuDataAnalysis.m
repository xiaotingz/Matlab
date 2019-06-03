file_1 = '/Users/xiaotingzhong/Desktop/Datas/Ni_simu/simu_An4__seg.dream3d';
% file_2 = '/Users/xiaotingzhong/Desktop/Datas/Ni_simu/simu_An5_clean_seg.dream3d';

grain_id_1 = squeeze(h5read(file_1, '/DataContainers/ImageDataContainer/CellData/FeatureIds'));
% grain_id_2 = squeeze(h5read(file_2, '/DataContainers/ImageDataContainer/CellData/FeatureIds'));


% tmp1 = P_an4(:,:,1);
% tmp2 = P_an5(:,:,1);
% 
% tmp = (P_an4 == P_an5);
% sum(tmp(:)) / (size(P_an4, 1) * size(P_an4, 2) * size(P_an4, 3)) * 100