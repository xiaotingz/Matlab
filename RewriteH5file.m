%% data structure of V4 & v6
% % V4
% facelabel = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
% curvature_of_triangle = abs(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'));
% GBCD_c1 = h5read(file_c1,'/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
% triangle_area = roundn(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-5);

% % V6
% facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
% curvature_of_triangle = abs(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'));
% triangle_area = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5);


%% Rewrite .dream3d file --- set triArea=0 for triangles on free surface

clear
file = '/Users/xiaotingzhong/Desktop/Datas/STO_1470/070617_V4_misA3ReconsAgrain/070617_sub2_reconsAgain_statsA0.dream3d';

facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
% facelabel = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels')).';
% curvature_of_triangle = abs(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'));
% triangle_area = roundn(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-5);
triangle_area = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5);

for i = 1:length(facelabel)
    if facelabel(i,1) <= 0 || facelabel(i,2) <=0
        triangle_area(i) = 0;
    end
end

% h5write(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures',curvature_of_triangle);
h5write(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas',triangle_area);

%% Rewrite .dream3d file --- seperate BB/BS/SS
% STO_1470, first calculate the average grain size for combining subset1
% and subset2

Size_sub1 = h5read('/Users/xiaotingzhong/Desktop/STO_1470/Modified/BB_BS_SS_GBCurvD_area0/subset1_BB_A0.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters');
Size_sub2 = h5read('/Users/xiaotingzhong/Desktop/STO_1470/Modified/BB_BS_SS_GBCurvD_area0/subset1_BB_A0.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters');

aveSize = (sum(Size_sub1) + sum(Size_sub2))/(length(Size_sub1) + length(Size_sub2));
% aveSize = 774.6
%% only for Big-Big grains
clear
readFile = '/Users/xiaotingzhong/Desktop/STO_1470/Modified/BB_BS_SS_GBCurvD_area0/subset2_BB_A0.dream3d';

facelabel = double(h5read(readFile,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels')).';
% curvature_of_triangle = abs(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'));
areas = roundn(h5read(readFile,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-5);
Size = h5read(readFile,'/VoxelDataContainer/FIELD_DATA/NumCells');

Size(1,:) = [];
aveSize = 775;

areas_filtered = zeros(size(areas));
for i = 1:length(facelabel)
    if facelabel(i,1) > 0 && facelabel(i,2) > 0
        if Size(facelabel(i,1)) > aveSize && Size(facelabel(i,2)) > aveSize
            areas_filtered(i) = areas(i);
        end
    end
end

h5write(readFile,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas',areas_filtered);

%% Rewrite .dream3d file --- set triArea!=0 only for Big-Small grains
clear
readFile = '/Users/xiaotingzhong/Desktop/STO_1470/Modified/BB_BS_SS_GBCurvD_area0/subset2_BS_A0.dream3d';

facelabel = double(h5read(readFile,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels')).';
% curvature_of_triangle = abs(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'));
areas = roundn(h5read(readFile,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-5);
Size = h5read(readFile,'/VoxelDataContainer/FIELD_DATA/NumCells');

Size(1,:) = [];
aveSize = 775;

areas_filtered = zeros(size(areas));
for i = 1:length(facelabel)
    if facelabel(i,1) > 0 && facelabel(i,2) > 0
        if Size(facelabel(i,1)) > aveSize && Size(facelabel(i,2)) <= aveSize
            areas_filtered(i) = areas(i);
        elseif Size(facelabel(i,1)) <= aveSize && Size(facelabel(i,2)) > aveSize
            areas_filtered(i) = areas(i);
        end
    end
end

h5write(readFile,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas',areas_filtered);





%% Rewrite .dream3d file --- set triArea!=0 only for Small-Small grains

clear
readFile = '/Users/xiaotingzhong/Desktop/STO_1470/Modified/BB_BS_SS_GBCurvD_area0/subset2_SS_A0.dream3d';

facelabel = double(h5read(readFile,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels')).';
% curvature_of_triangle = abs(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'));
areas = roundn(h5read(readFile,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-5);
Size = h5read(readFile,'/VoxelDataContainer/FIELD_DATA/NumCells');

Size(1,:) = [];
aveSize = 775;

areas_filtered = zeros(size(areas));
for i = 1:length(facelabel)
    if facelabel(i,1) > 0 && facelabel(i,2) > 0
        if Size(facelabel(i,1)) <= aveSize && Size(facelabel(i,2)) <= aveSize
            areas_filtered(i) = areas(i);
        end
    end
end

h5write(readFile,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas',areas_filtered);

%% END of BB/BS/SS



%% Rewrite .dream3d file --- set triArea=0 for triangles on TripleLines
clear
% read in data
file = '/Volumes/RESEARCH/Dec.4 TripleLineSpheres/no_TriLineTri/Jul20_bad_noTLT.dream3d';
nodeInfo = textread('/Volumes/RESEARCH/Dec.4 TripleLineSpheres/no_TriLineTri/Jul20_A_node.txt');
triInfo = textread('/Volumes/RESEARCH/Dec.4 TripleLineSpheres/no_TriLineTri/Jul20_A_triangle.txt');

triangle_area = roundn(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-5);
% first lines are #nodes and #triangles, not useful and interrupt loop, delete
nodeInfo(1,:) = [];
triInfo(1,:) = [];
% modify triangleID and nodeID to start from 1, for MATLAB convenience
nodeInfo(:,1) = nodeInfo(:,1) + 1;
triInfo(:,1) = triInfo(:,1) + 1;
triInfo(:,2) = triInfo(:,2) + 1;
triInfo(:,3) = triInfo(:,3) + 1;
triInfo(:,4) = triInfo(:,4) + 1;

% nodeInfo: 2nd column is nodeType
    % if only 2 is used, triangles on inside boundaries and not on triple lines nor quad points
    % if 2, 3, 4 are used, triangles on all inside boundaries
% triInfo: 2nd-4th column is the nodeID for its nodes

% prepare a boolean list to record if keep this triangle or not
% start with 1s! (if start with 0s: one of the node satisfy, triangle keep)
boolean = ones(length(triInfo),1);
% for every triangle, get its three nodeID, check the nodeType, if not 2, delete
for i = 1 : length(triInfo)
    for j = 2 : 4
        tmp = triInfo(i,j);
        if nodeInfo(tmp,2) ~= 2
            triangle_area(i) = 0;
        end
    end
end

h5write(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas',triangle_area);


% areaZero2 = roundn(h5read('/Volumes/RESEARCH/Grain Curvature/setTo_a0/Jul20_A_a0.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-5);









%% Analyze characterization BinValue
% file_c1a0 = '/Volumes/RESEARCH/Simple Geometry/Oct.7 No1-10/setTo1/CurvTo1_AreaTo0/No10_c1a0_CurvDistri.dream3d';
file_c1 = '/Volumes/RESEARCH/Simple Geometry/Oct.7 No1-10/setTo1/No10_qm_CurvatureTo1/No10_qm_TriCurv=1_CurvDistri.dream3d';

% GBCD_c1a0 = h5read(file_c1a0,'/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
% Counter_c1a0 = h5read(file_c1a0,'/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCDCounters');
GBCD_c1 = h5read(file_c1,'/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
Counter_c1 = h5read(file_c1,'/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCDCounters');

% holderB_c1a0 = GBCD_c1a0(:,1);
holderB_c1 = GBCD_c1(:,1);
bool = zeros(size(holderB_c1));
for i = 1:length(holderB_c1)
    if holderB_c1(i) ~= 0
        bool(i) = 1;
    end
end
bool = logical(bool);

% Bins_c1a0 = holderB_c1a0(bool);
Bins_c1 = holderB_c1(bool);
% holderC_c1a0 = Counter_c1a0(:,1);
holderC_c1 = Counter_c1(:,1);
% Counter_c1a0 = holderC_c1a0(bool);
Counter_c1 = holderC_c1(bool);

%% double check triangle set
file_c1a0 = '/Volumes/RESEARCH/Simple Geometry/Oct.7 No1-10/setTo1/CurvTo1_AreaTo0/No10_c1a0.dream3d';

facelabel = double(h5read(file_c1a0,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels')).';
curvature_of_triangle = abs(h5read(file_c1a0,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'));

bool_label = zeros(length(facelabel),1);
bool_curvature = zeros(length(facelabel),1);

for i = 1 : length(facelabel)
    if (facelabel(i,1) == 1 && facelabel(i,2) == 2) || (facelabel(i,1) == 2 && facelabel(i,2) == 1)
        bool_label(i) = 1;
    elseif curvature_of_triangle(i) ~= 0
        bool_curvature(i) = 1;
    end
end




