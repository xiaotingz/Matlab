file_1 = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub1_reconsHsmooth.dream3d');
file_2 = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub2_reconsHsmooth.dream3d');

h5write(file_1, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList', xsmooth_1);
h5write(file_2, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList', xsmooth_2);


