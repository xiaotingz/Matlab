% file_1 = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub1_recons.dream3d');
% xdat_1 = double(h5read(file_1, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'));
% tri_1 = 1 + double(h5read(file_1, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')).';
% fl_1 = double(h5read(file_1, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
% ntype_1 = double(h5read(file_1, '/DataContainers/TriangleDataContainer/VertexData/NodeType'));
% 
% 
% f_1 = find( any( fl_1==-1 | fl_1==0, 2 ) );
% fl_1( f_1, : ) = [];
% tri_1( f_1, : ) = [];

% xsmooth_1 = HierarchicalSmooth( xdat_1, tri_1, fl_1, ntype_1 );

% file_2 = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub2_recons.dream3d');
% xdat_2 = double(h5read(file_2, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'));
% tri_2 = 1 + double(h5read(file_2, '/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList')).';
% fl_2 = double(h5read(file_2, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
% ntype_2 = double(h5read(file_2, '/DataContainers/TriangleDataContainer/VertexData/NodeType'));

% f_2 = find( any( fl_2==-1 | fl_2==0, 2 ) );
% fl_2( f_2, : ) = [];
% tri_2( f_2, : ) = [];

% xsmooth_2 = HierarchicalSmooth( xdat_2, tri_2, fl_2, ntype_2 );


xdat = load( '/Users/xiaotingzhong/Desktop/programs/HierarchicalSmooth_matlab/examples/ex0/SharedVertexList.txt' );
xdat = xdat';
tri = int32(1 + load( '/Users/xiaotingzhong/Desktop/programs/HierarchicalSmooth_matlab/examples/ex0/SharedTriList.txt' ));
fl = int32(load( '/Users/xiaotingzhong/Desktop/programs/HierarchicalSmooth_matlab/examples/ex0/FaceLabels.txt' ));
f = find( any( fl==-1 | fl==0, 2 ) );
fl( f, : ) = [];
tri( f, : ) = [];
ntype = int32(load( '/Users/xiaotingzhong/Desktop/programs/HierarchicalSmooth_matlab/examples/ex0/NodeType.txt' ));
% HierarchicalSmoothFrontEnd(tri, xdat, fl, ntype) 


% 
% save('180311_STO1470sub1_xsmooth', 'xsmooth_1');
% save('180311_STO1470sub2_xsmooth', 'xsmooth_2');