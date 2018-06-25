
% --- why there is no TL=[5, 71,129] in the synthetic volume? ---

%  ##### Read in information #####
% --- dims_N = dims_sample + 1 ---
dims_N = [129, 129, 129];
% --- get the objQN ---
objQN = [5, 71, 129, 235];
% --- D3D info ---
file = '/Users/xiaotingzhong/Desktop/Datas/synthetic/180502_CubicSingleEquiaxedOut.dream3d';
NType_D3D = double(h5read(file,'/DataContainers/TriangleDataContainer/VertexData/NodeType'))';
NC_D3D = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
NC_D3D = [(1:length(NC_D3D))', NC_D3D];
triNodes = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
mask_QN_D3D = (NType_D3D == 4);
QNC_D3D = NC_D3D(mask_QN_D3D, :);
mask_TN_D3D = (NType_D3D == 3);
TNC_D3D = NC_D3D(mask_TN_D3D, :);

FLs = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';


%  ##### objQN's ID in D3D #####
QNList(2:5) = sort(QNList(2:5), 2);
mask_objQN = ismember(QNList(:,2:5), objQN, 'rows');
objQN = QNList(mask_objQN, :);
[X, Y, Z] = ind2sub(dims_N, objQN(:,1));
objN_a = ([X, Y, Z] - 1)*0.5;
mask_objNa_D3DID = ismember(QNC_D3D(:, 2:4), objN_a, 'rows');
objN_a_D3DID = QNC_D3D(mask_objNa_D3DID);
% --- Get the other node from PARAVIEW, here it'a TL, may need CHANGE ---
objN_b = ([X, Y-1, Z] - 1)*0.5;
mask_objNb_D3DID = ismember(TNC_D3D(:, 2:4), objN_b, 'rows');
objN_b_D3DID = TNC_D3D(mask_objNb_D3DID);
% --- Get face labels of triangles sharing the two objNodes ---
mask_objTris = ((any(triNodes==objN_a_D3DID(1), 2) + any(triNodes==objN_b_D3DID(1), 2)) == 2);
FLs_objTris = FLs(mask_objTris, :)


% clearvars -except FLs_objTris
