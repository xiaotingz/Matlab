file = '/Users/xiaotingzhong/Desktop/Datas/Mg/05_GB_crystallography.dream3d';

num_neigh = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
size = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters')).';

num_neigh(1) = [];
size(1) = [];

id = (1:length(size))';

data = [id, size, num_neigh];
data = sortrows(data, [2, 3], 'descend');


csvwrite('Mg_5_GB_grains.txt', data)

