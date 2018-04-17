% --------------------------- write the smoothed coordinates back --------------------------- 
% coorFile1 = load('/Users/xiaotingzhong/Desktop/Datas/STO_1350C_4volumes_forXiaoting/20180112_D3Dfiles_andTriangles/0729_xsmooth.mat');
newCoor = xsmooth_2;
d3dfile1 = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub2_recons.dream3d');
% coordinates1 = h5read(d3dfile1,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList');
h5write(d3dfile1,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList',single(newCoor));

% --------------------------- get the biggest grain ID --------------------------- 
d = h5read(d3dfile1,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters');
[~, ind] = sort(d, 'descend');

'sub2, largest grain ID:'







