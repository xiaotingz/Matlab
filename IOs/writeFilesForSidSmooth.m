clear
InFile = ('/Users/xiaotingzhong/Desktop/Datas/STO_1350C_4volumes_forXiaoting/20180112_D3Dfiles_andTriangles/STO1350C_D3DFiles/meancurvature_STO_1350C_400nm_0802.dream3d');
OutFile_node = ('/Users/xiaotingzhong/Desktop/Datas/STO_1350C_4volumes_forXiaoting/20180112_D3Dfiles_andTriangles/STO1350C_D3DFiles/meancurvature_STO_1350C_400nm_0802_nodeType.txt');
OutFile_facelabel = ('/Users/xiaotingzhong/Desktop/Datas/STO_1350C_4volumes_forXiaoting/20180112_D3Dfiles_andTriangles/STO1350C_D3DFiles/meancurvature_STO_1350C_400nm_0802_facelabel.txt');

nodeType = h5read(InFile, '/DataContainers/TriangleDataContainer/VertexData/NodeType').';
% Notice the size of facelabel=[2, n] in order to be written correctly by fprintf
facelabel = h5read(InFile, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels');

fileID = fopen(OutFile_node,'w');
fprintf(fileID, '%d\n', nodeType);
fclose(fileID);

fileID = fopen(OutFile_facelabel,'w');
fprintf(fileID, '%d  %d\n', facelabel);
fclose(fileID);