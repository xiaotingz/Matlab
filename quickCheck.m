file = '/Users/xiaotingzhong/Desktop/Datas/STO_1470/070617_V4_misA3ReconsAgrain/070617_sub2_reconsAgain10.dream3d';
num_of_neigh = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
max(num_of_neigh)
misA = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/MisorientationList'));
min(misA)

misA = sort(misA);
misA(1:5)