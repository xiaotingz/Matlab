file_1 = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub1_reconsHsmooth.dream3d');
file_2 = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub2_reconsHsmooth.dream3d');

area_1 = double(h5read(file_1, '/DataContainers/TriangleDataContainer/FaceData/FaceAreas'))';
area_2 = double(h5read(file_2, '/DataContainers/TriangleDataContainer/FaceData/FaceAreas'))';
fl_1 = double(h5read(file_1, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
fl_2 = double(h5read(file_2, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';

triAreas_1 = [fl_1, area_1];
triAreas_2 = [fl_2, area_2];
triAreas_1 = sortrows(triAreas_1, -3);
triAreas_2 = sortrows(triAreas_2, -3);

d = h5read(file_1,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters');
[~, ind] = sort(d, 'descend');






