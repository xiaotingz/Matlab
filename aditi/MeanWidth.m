% Mean width of a grain calculation
clear
file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_crop.dream3d');

%  -------------------------- load data -------------------------- 
SharedEdgeList = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedEdgeList')).';%edge formed between vertices
SharedTriList = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')).';
SharedVertexList = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList')).';
nodetype =double(h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType')).';
Facelabel_triangles = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
FaceNormals = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceNormals')).';
FeaturefaceId= double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FeatureFaceId')).';
FaceLabel_grain= double(h5read(file,'/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')).';
Numtriangles = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceFeatureData/NumTriangles')).';
tri_curv = h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures').';
tri_area = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-8).';

 vertex_data = [SharedVertexList, nodetype];
 %colnames1 = {'X','Y','Z','Node'};
 %vertex_data=array2table(vertex_data,'VariableNames',colnames1);%containes the x,y,z coordinates of each vertex and the nodetype.
 triangle_data = [Facelabel_triangles,SharedTriList, FaceNormals,FeaturefaceId,tri_curv, tri_area];
 %colnames2={'Flabel1','Flabel2','V1','V2','V3','N1','N2','N3','FaceId','triangle curvature','triangle area'};
% triangle_data = array2table(triangle_data,'VariableNames',colnames2)
 
% %finding the common triangles to every edge
% C(:,1:2)= intersect(SharedTriList(:,1:3),SharedEdgeList(:,1:2))

%edge lengths of each triangle
Edge_data=calcEdgeLength(SharedEdgeList,SharedVertexList);

% %calculating triangle IDs
% tri_id = triangleID(triangle_data);
% 
% 









