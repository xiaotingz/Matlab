file_1 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6.dream3d');
file_2 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_mesh.dream3d');

avgEA_1 = round(double(h5read(file_1,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters')),2);
avgEA_2 = round(double(h5read(file_2,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters')),2);

comp = (avgEA_1 == avgEA_2);
sum(comp(:))