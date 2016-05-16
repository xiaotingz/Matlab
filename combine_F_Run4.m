%% average --- CurvDistri file 
% load data
Ferrite_1 = ('/Users/xiaotingzhong/Desktop/Datas/setTo0/Jan31_Fca0_curDistri10.dream3d');
Run4_2 = ('/Users/xiaotingzhong/Desktop/Datas/Ferrite/Ferrite_Run4/setTo0/Jan31_F_Run4_CurvDistri10.dream3d');

CurvDistri1 = h5read(Ferrite_1,'/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
Counter1 = h5read(Ferrite_1,'/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCDCounters');
CurvDistri2 = h5read(Run4_2,'/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
Counter2 = h5read(Run4_2,'/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCDCounters');

% average data from the two volumes
aveCurvDistri = (CurvDistri1(:,1).*Counter1(:,1) + CurvDistri2(:,1).*Counter2(:,1)) ./ (Counter1(:,1) + Counter2(:,1));

% modify index to use gbcd_graph
load('TwoIndex_res10.mat')
    % link the two system. Sort as gbcd_index in its natural order, then graph_index can be used to place D3D gbcd data
d3d_to_graph = [gbcd_index,graph_index];
d3d_to_graph = sortrows(d3d_to_graph);

% % Prepare to converget: read GBCD data in D3D.
data_read = aveCurvDistri;

data_converted = zeros(length(data_read),1);
for i = 1:length(data_read)
    data_converted(d3d_to_graph(i,2)) = data_read(i);
end
% 
% fileID = fopen('F_CurvDistri10Ave.txt','w');
% fprintf(fileID,'%12.8f\n',data_converted);
% fclose(fileID);


calculated = textread('F_CurvDistri10Ave.txt');

%% average --- GBCD file
Ferrite_1 = ('/Users/xiaotingzhong/Desktop/Datas/setTo0/Jan31_Fca0_gbcd10.dream3d');
Run4_2 = ('/Users/xiaotingzhong/Desktop/Datas/Ferrite/Ferrite_Run4/setTo0/Jan31_F_Run4_GBCD10.dream3d');

GBCD1 = h5read(Ferrite_1,'/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
GBCD2 = h5read(Run4_2,'/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
triArea1 = roundn(h5read(Ferrite_1,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-8);
triArea2 = roundn(h5read(Run4_2,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-8);

MRD1 = length(GBCD1)/sum(triArea1);
MRD2 = length(GBCD2)/sum(triArea2);
MRDave = length(GBCD1)/(sum(triArea1)+sum(triArea2));

aveGBCD = (GBCD1/MRD1 + GBCD2/MRD2) * MRDave;

% convert to gbcd_graph order
load('TwoIndex_res10.mat')
    % link the two system. Sort as gbcd_index in its natural order, then graph_index can be used to place D3D gbcd data
d3d_to_graph = [gbcd_index,graph_index];
d3d_to_graph = sortrows(d3d_to_graph);

% % Prepare to converget: read GBCD data in D3D.
data_read = aveGBCD;

data_converted = zeros(length(data_read),1);
for i = 1:length(data_read)
    data_converted(d3d_to_graph(i,2)) = data_read(i);
end

% fileID = fopen('F_gbcd10Ave.txt','w');
% fprintf(fileID,'%12.8f\n',data_converted);
% fclose(fileID);

%% Stick data in DREAM3D file 
% start with (Jan.31 Ferrite) data, stick Run4's data at the end of each list then let D3D compute
% Needed data: label, normal, eulerAngle, curvature, eulerAngle
Ferrite_1 = ('/Users/xiaotingzhong/Desktop/Datas/Ferrite/Ferrite_stick.dream3d');
Run4_2 = ('/Users/xiaotingzhong/Desktop/Datas/Ferrite/Ferrite_Run4/setTo0/Jan31_F_Run4_setTo0.dream3d');

facelabel1 = double(h5read(Ferrite_1,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
triNormal1 = double(h5read(Ferrite_1,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceNormals'));
euler1 = h5read(Ferrite_1,'/VoxelDataContainer/FIELD_DATA/EulerAngles');
triCurv1 = h5read(Ferrite_1,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures');
triArea1 = roundn(h5read(Ferrite_1,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-8);

facelabel2 = double(h5read(Run4_2,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
triNormal2 = double(h5read(Run4_2,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceNormals'));
euler2 = h5read(Run4_2,'/VoxelDataContainer/FIELD_DATA/EulerAngles');
triCurv2 = h5read(Run4_2,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures');
triArea2 = roundn(h5read(Run4_2,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-8);

% need to change Run4's faceLabel to make it correctly ID for finding eulerAngle
facelabel2 = facelabel2.';
for i = 1:length(facelabel2)
    if facelabel2(i,1) > 0 && facelabel2(i,2) > 0
        facelabel2(i,1) = facelabel2(i,1) + 1114;
        facelabel2(i,2) = facelabel2(i,2) + 1114;
    end
end
facelabel2 = facelabel2.';

facelabel = [facelabel1,facelabel2];
triNormal = [triNormal1,triNormal2];
triCurv = [triCurv1;triCurv2];
triArea = [triArea1;triArea2];

% eulerAngle is special becase grainID is changed
euler2(:,1) = [];
euler = [euler1,euler2];

% Rewrite HDF5 file
h5write(Ferrite_1,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels',facelabel);
h5write(Ferrite_1,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceNormals',triNormal);
h5write(Ferrite_1,'/VoxelDataContainer/FIELD_DATA/EulerAngles',euler);
h5write(Ferrite_1,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures',triCurv);
h5write(Ferrite_1,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas',triArea);

% extent data size in HDF5 file
fileattrib(Ferrite_1,'+w')
fileName = H5F.open(Ferrite_1,'H5F_ACC_RDWR','H5P_DEFAULT');
dset_id = H5D.open(fileName,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures');
H5D.set_extent(dset_id,5000000);
H5D.close(dset_id);
H5F.close(fileName);
