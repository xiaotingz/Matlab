curvFile = ('/Users/xiaotingzhong/Desktop/Datas/STO_1350C_4volumes_forXiaoting/mesh_STO_1350C_400nm_0729_rohrerpipe_curv.dream3d');
OutFile = strcat('/Users/xiaotingzhong/Desktop/Datas/STO_1350C_4volumes_forXiaoting/mesh_STO_1350C_400nm_0729_rohrerpipe_curvTriangles.ph');


FaceLabel_raw = h5read(curvFile, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels').';
curv_raw = h5read(curvFile, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures').';
tri_bool = ones(length(curv_raw),1);
for i = 1:length(FaceLabel_raw)
    if FaceLabel_raw(i,1) < 0 || FaceLabel_raw(i,2) < 0
        tri_bool(i) = 0;
    end
end
tri_bool = logical(tri_bool);
curv = curv_raw(tri_bool);
curv = abs(curv);

fileID = fopen(OutFile,'w');
fprintf(fileID,'# Triangles Produced from DREAM3D version 6.4\n');
fprintf(fileID,'# Column 1-3:    right hand average orientation (phi1, PHI, phi2 in RADIANS)\n');
fprintf(fileID,'# Column 4-6:    left hand average orientation (phi1, PHI, phi2 in RADIANS)\n');
fprintf(fileID,'# Column 7-9:    triangle normal\n');
fprintf(fileID,'# Column 10:      surface area\n');
fprintf(fileID,'# Column 11:      triangle curvature\n');

format = '%7.4f    %7.4f    %7.4f    %7.4f    %7.4f  %7.4f  %7.4f  %7.4f  %7.4f   %7.4f   %7.4f\n';
for i = 1:length(VarName1)
    fprintf(fileID,format,VarName1(i),VarName2(i),VarName3(i),VarName4(i),VarName5(i),VarName6(i), ...
        VarName7(i),VarName8(i),VarName9(i),VarName10(i), curv(i));
end
fclose('all');
%% This is for the STO_1350 data, where there are facelabel=0
curvFile = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub2_GBCD.dream3d');
OutFile = strcat('/Users/xiaotingzhong/Desktop/Datas/180806/STOsub2_full.ph');


FaceLabel_raw = h5read(curvFile, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels').';
curv_raw = h5read(curvFile, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures').';
tri_bool = ones(length(curv_raw),1);
for i = 1:length(FaceLabel_raw)
    if FaceLabel_raw(i,1) < 0 || FaceLabel_raw(i,2) < 0
        tri_bool(i) = 0;
    end
end
tri_bool = logical(tri_bool);
curv = curv_raw(tri_bool);
curv = abs(curv);

FaceLabel = FaceLabel_raw(tri_bool,:);
tri_bool2 = ones(length(FaceLabel),1);
for i = 1: length(FaceLabel)
    if FaceLabel(i,1) == 0 || FaceLabel(i,2) == 0
        tri_bool2(i) = 0;
    end
end

fileID = fopen(OutFile,'w');
fprintf(fileID,'# Triangles Produced from DREAM3D version 6.4\n');
fprintf(fileID,'# Column 1-3:    right hand average orientation (phi1, PHI, phi2 in RADIANS)\n');
fprintf(fileID,'# Column 4-6:    left hand average orientation (phi1, PHI, phi2 in RADIANS)\n');
fprintf(fileID,'# Column 7-9:    triangle normal\n');
fprintf(fileID,'# Column 10:      surface area\n');
fprintf(fileID,'# Column 11:      triangle curvature\n');

format = '%7.4f    %7.4f    %7.4f    %7.4f    %7.4f  %7.4f  %7.4f  %7.4f  %7.4f   %7.4f   %7.4f\n';
for i = 1:length(VarName1)
    if tri_bool2(i) > 0
        fprintf(fileID,format,VarName1(i),VarName2(i),VarName3(i),VarName4(i),VarName5(i),VarName6(i), ...
            VarName7(i),VarName8(i),VarName9(i),VarName10(i), curv(i));
    end
end
fclose('all');




%% ################################ Combining triangles in two volumes ################################ 
file_1 = '/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311_STO1470sub1_GBCD.dream3d';
file_2 = '/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311_STO1470sub2_GBCD.dream3d';

EA_1 = h5read(file_1, '/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles').';
EA_2 = h5read(file_2, '/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles').';
EA_1(1, :) = [];
EA_2(1, :) = [];

fl_1 = h5read(file_1, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels').';
normal_1 = h5read(file_1, '/DataContainers/TriangleDataContainer/FaceData/FaceNormals').';
area_1 = h5read(file_1, '/DataContainers/TriangleDataContainer/FaceData/FaceAreas').';
curv_1 = abs(h5read(file_1, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures').');
fl_2 = h5read(file_2, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels').';
normal_2 = h5read(file_2, '/DataContainers/TriangleDataContainer/FaceData/FaceNormals').';
area_2 = h5read(file_2, '/DataContainers/TriangleDataContainer/FaceData/FaceAreas').';
curv_2 = abs(h5read(file_2, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures').');

mask = all(fl_1 > 0, 2) & curv_1 < 100;
fl_1 = fl_1(mask, :);
normal_1 = normal_1(mask, :);
area_1 = area_1(mask, :);
curv_1 = curv_1(mask, :);
mask = all(fl_2 > 0, 2) & curv_2 < 100;
fl_2 = fl_2(mask, :);
normal_2 = normal_2(mask, :);
area_2 = area_2(mask, :);
curv_2 = curv_2(mask, :);

format = '%7.4f    %7.4f    %7.4f    %7.4f    %7.4f  %7.4f  %7.4f  %7.4f  %7.4f   %7.4f   %7.4f\n';
fileID = fopen('STO_GBHDtris_combined.txt','w');
for i = 1:length(fl_1)
    g1 = fl_1(i, 1);
    g2 = fl_1(i, 2);
    fprintf(fileID, format, EA_1(g1, 1),EA_1(g1, 2),EA_1(g1, 3),EA_1(g2, 1),EA_1(g2, 2),EA_1(g2, 3), ...
            normal_1(i, 1), normal_1(i, 2), normal_1(i, 3), area_1(i), curv_1(i));
end

for i = 1:length(fl_2)
    g1 = fl_2(i, 1);
    g2 = fl_2(i, 2);
    fprintf(fileID, format, EA_2(g1, 1),EA_2(g1, 2),EA_2(g1, 3),EA_2(g2, 1),EA_2(g2, 2),EA_2(g2, 3), ...
            normal_2(i, 1), normal_2(i, 2), normal_2(i, 3), area_2(i), curv_2(i));
end
fclose(fileID);









