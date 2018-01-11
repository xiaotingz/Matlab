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
curvFile = ('/Users/xiaotingzhong/Desktop/Datas/STO_1350C_4volumes_forXiaoting/Mesh_STO_1350C_0802_v62_curv.dream3d');
OutFile = strcat('/Users/xiaotingzhong/Desktop/Datas/STO_1350C_4volumes_forXiaoting/Mesh_STO_1350C_0802_v62_curvTriangles.ph');


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



