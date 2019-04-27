% #####################################################
% data shape for the c++ script
%     fl, ntype, tri = [n, 2/1/3]
%     xdat = [3, n]
% #####################################################
file_1 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/D3Ds/An4new6_fixOrigin3_Hsmooth.dream3d');
xdat_1 = double(h5read(file_1, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'));
tri_1 = 1 + double(h5read(file_1, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')).';
fl_1 = double(h5read(file_1, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
ntype_1 = double(h5read(file_1, '/DataContainers/TriangleDataContainer/VertexData/NodeType'))';


f_1 = find( any( fl_1 <= 0, 2 ) );
fl_1( f_1, : ) = [];
tri_1( f_1, : ) = [];

% xsmooth_1 = HierarchicalSmooth( xdat_1, tri_1, fl_1, ntype_1 );

% file_2 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/D3Ds/An5new6_Hsmooth.dream3d');
% xdat_2 = double(h5read(file_2, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'));
% tri_2 = 1 + double(h5read(file_2, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
% fl_2 = double(h5read(file_2, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
% ntype_2 = double(h5read(file_2, '/DataContainers/TriangleDataContainer/VertexData/NodeType'))';
% 
% f_2 = find( any( fl_2 <= 0, 2 ) );
% fl_2( f_2, : ) = [];
% tri_2( f_2, : ) = [];

% xsmooth_2 = HierarchicalSmooth( xdat_2, tri_2, fl_2, ntype_2 );

% h5write(file_1, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList', xsmooth_an4);
% h5write(file_2, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList', xsmooth_an5);


%% ##### extract data of a single grain #####
file_1 = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub1_recons.dream3d');
xdat_1 = double(h5read(file_1, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'));
tri_1 = 1 + double(h5read(file_1, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')).';
fl_1 = double(h5read(file_1, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
ntype_1 = double(h5read(file_1, '/DataContainers/TriangleDataContainer/VertexData/NodeType'))';

ID = 407;
mask = find( any( fl_1 ==ID, 2) );
fl = fl_1( mask, : );
tri = tri_1( mask, : );
xdat = xdat_1;
ntype = ntype_1;

% trisurf( tri, xdat(1,:), xdat(2,:), xdat(3,:) );
% rotate3d on