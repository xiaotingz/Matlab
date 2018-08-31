num_of_neigh = h5read('E:\Feb.20\uncropped\Feb.20 sf4.dream3d','/VoxelDataContainer/FIELD_DATA/NumNeighbors');
facelabel = h5read('E:\Feb.20\uncropped\Feb.20 sf4.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels').';
% neighborList = h5read('E:\Feb.20\uncropped\Feb.20 sf4.dream3d','/VoxelDataContainer/FIELD_DATA/NeighborList');

neighof_badgrain = max(num_of_neigh);
[r c] = find(num_of_neigh == neighof_badgrain);  % find bad data

data_raw = facelabel.';
tmp1 = data_raw(1,:);
tmp2 = data_raw(2,:);

boolean1 = (data_raw(1,:) == 0 | data_raw(2,:) == 0);
data_with0(1,:) = tmp1(boolean1);
data_with0(2,:) = tmp2(boolean1);
percent_0 = length(data_with0)/length(facelabel);

boolean2 = (data_raw(1,:) < 0 | data_raw(2,:) < 0);
data_withneg3(1,:) = tmp1(boolean2);
data_withneg3(2,:) = tmp2(boolean2);
percent_neg3 = length(data_withneg3)/length(facelabel);

%for bad datas, num_of_neighbors started from grainid=0 ==> r-1
boolean3 = (tmp1 == r-1 | tmp2 == r-1);
data_with_badgrain(1,:) = tmp1(boolean3);
data_with_badgrain(2,:) = tmp2(boolean3);
percent_badgrain = length(data_with_badgrain)/length(facelabel);

clear c tmp1 tmp2 tmp3 boolean1  boolean2 boolean3


% Dec.3:  E:\Dream3D\Datas\12_5_2014\Austenite - 3D TWIP_Jan2013_SQ ANG Files_Output\before Dec.16\Dec_3_curvature_laplacian.dream3d
% austenite-curvaturetry:   E:\Dream3D\Datas\12_5_2014\Austenite - 3D TWIP_Jan2013_SQ ANG Files_Output\before Dec.16\austenite-curvaturetry.dream3d
% Dec. 16:  E:\Dream3D\Datas\12_5_2014\Austenite - 3D TWIP_Jan2013_SQ ANG Files_Output\Dec.16--crashes when computing curvature\Dec.16_surfacemesh04.dream3d
% Dec.9:   E:\Dream3D\Datas\12_5_2014\Austenite - 3D TWIP_Jan2013_SQ ANG Files_Output\E:\Dream3D\Datas\12_5_2014\Austenite - 3D TWIP_Jan2013_SQ ANG Files_Output\Dec.16--crashes when computing curvature\Dec.16 statistics06.dream3d