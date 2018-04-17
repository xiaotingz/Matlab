% Instruction
% - This script is to combine seperate EBSD volumes into a large for for computing grain boundary curvature distribution(GBCurvD).
% - Before running this script, the user should have two files ready. 1. the GBCD triangles file from DREAM3D. 
%      2. the corresponding DREAM3D file that cotain the triangle curvatures. 
% - Note the GBCD triangles file's header need to be changed before running this script. 
%     Replace the header part of each txt file with the following (each line for a volume):
%         EA1_1  EA1_2  EA1_3  EA2_1  EA2_2  EA2_3  n1     n2     n3      area1
%         EA1_4  EA1_5  EA1_6  EA2_4  EA2_5  EA2_6  n4     n5     n6      area2
%         EA1_7  EA1_8  EA1_9  EA2_7  EA2_8  EA2_9  n7     n8     n9      area3
%         EA1_10 EA1_11 EA1_12 EA2_10 EA2_11 EA2_12 n10    n11    n12    area4
% - This script is designed to combine 4 seperate volumes. If the number of volumes to combine is different, the script needs to be changed.
% - BEFORE RUNNING THE SCRIPT:
%     1. change the directory of the .dream3d file (include the file name itself)
%             Note this scipt is for DREAM3D version6. If different version is used, may need to change the variable tree in 'h5read' command.
%     2. change the header of the GBCD triangle files and input each volume by MATLAB's 'Import Data' wizard. 
%             if it's not reading the columns correctly, change the 'Column deliminters', should be space.
%     3. change the fileID (at the end of the script), which would be the output file.
%% inputs that need to be specified
curvFile_1 = '/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub1_GBCD.dream3d';
curvFile_2 = '/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub2_tillSlice29_GBCD.dream3d';
% curvFile_3 = '/Users/xiaotingzhong/Desktop/Datas/STO_1350C_4volumes_forXiaoting/111217/STO1350C_D3DFiles/meancurvature_STO_1350C_400nm_0801.dream3d';
% curvFile_4 = '/Users/xiaotingzhong/Desktop/Datas/STO_1350C_4volumes_forXiaoting/111217/STO1350C_D3DFiles/meancurvature_STO_1350C_400nm_0802.dream3d';


% -------------------------------------------------------------------
% reading in the curvature data and make it the same length as the other variables.
FaceLabel_raw_1 = h5read(curvFile_1, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels').';
curv_raw_1 = h5read(curvFile_1, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures').';
tri_bool_1 = ones(length(curv_raw_1),1);
for i = 1:length(FaceLabel_raw_1)
    if FaceLabel_raw_1(i,1) < 0 || FaceLabel_raw_1(i,2) < 0
        tri_bool_1(i) = 0;
    end
end
tri_bool_1 = logical(tri_bool_1);
curv_1 = curv_raw_1(tri_bool_1);
curv_1 = abs(curv_1);

FaceLabel_raw_2 = h5read(curvFile_2, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels').';
curv_raw_2 = h5read(curvFile_2, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures').';
tri_bool_2 = ones(length(curv_raw_2),1);
for i = 1:length(FaceLabel_raw_2)
    if FaceLabel_raw_2(i,1) < 0 || FaceLabel_raw_2(i,2) < 0
        tri_bool_2(i) = 0;
    end
end
tri_bool_2 = logical(tri_bool_2);
curv_2 = curv_raw_2(tri_bool_2);
curv_2 = abs(curv_2);

% FaceLabel_raw_3 = h5read(curvFile_3, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels').';
% curv_raw_3 = h5read(curvFile_3, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures').';
% tri_bool_3 = ones(length(curv_raw_3),1);
% for i = 1:length(FaceLabel_raw_3)
%     if FaceLabel_raw_3(i,1) < 0 || FaceLabel_raw_3(i,2) < 0
%         tri_bool_3(i) = 0;
%     end
% end
% tri_bool_3 = logical(tri_bool_3);
% curv_3 = curv_raw_3(tri_bool_3);
% curv_3 = abs(curv_3);
% 
% FaceLabel_raw_4 = h5read(curvFile_4, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels').';
% curv_raw_4 = h5read(curvFile_4, '/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures').';
% tri_bool_4 = ones(length(curv_raw_4),1);
% for i = 1:length(FaceLabel_raw_4)
%     if FaceLabel_raw_4(i,1) < 0 || FaceLabel_raw_4(i,2) < 0
%         tri_bool_4(i) = 0;
%     end
% end
% tri_bool_4 = logical(tri_bool_4);
% curv_4 = curv_raw_4(tri_bool_4);
% curv_4 = abs(curv_4);

% curv = [curv_1; curv_2; curv_3; curv_4];
curv = [curv_1; curv_2];


% delete the unnecessary data
% clear curv_1 curv_2 curv_3 curv_4 curv_raw_1 curv_raw_2 curv_raw_3 curv_raw_4 ...
%     FaceLabel_1 FaceLabel_2 FaceLabel_3 FaceLabel_4 curvFile_1 curvFile_2 curvFile_3 curvFile_4;

% EA_A1 = [EA1_1; EA1_4; EA1_7; EA1_10];
% EA_A2 = [EA1_2; EA1_5; EA1_8; EA1_11];
% EA_A3 = [EA1_3; EA1_6; EA1_9; EA1_12];
% EA_B1 = [EA2_1; EA2_4; EA2_7; EA2_10];
% EA_B2 = [EA2_2; EA2_5; EA2_8; EA2_11];
% EA_B3 = [EA2_3; EA2_6; EA2_9; EA2_12];
% n_1 = [n1; n4; n7; n10];
% n_2 = [n2; n5; n8; n11];
% n_3 = [n3; n6; n9; n12];
% area = [area1; area2; area3; area4];
EA_A1 = [EA1_1; EA1_4];
EA_A2 = [EA1_2; EA1_5];
EA_A3 = [EA1_3; EA1_6];
EA_B1 = [EA2_1; EA2_4];
EA_B2 = [EA2_2; EA2_5];
EA_B3 = [EA2_3; EA2_6];
n_1 = [n1; n4];
n_2 = [n2; n5];
n_3 = [n3; n6];
area = [area1; area2];


toPrint_raw = [EA_A1, EA_A2, EA_A3, EA_B1, EA_B2, EA_B3, n_1, n_2, n_3, area, curv];

% -------------------------------------------------------------------
% get the bool array to exclude all the EA=0 data
FaceLabel_1 = FaceLabel_raw_1(tri_bool_1,:);
tri_bool2_1 = ones(length(FaceLabel_1),1);
for i = 1: length(FaceLabel_1)
    if FaceLabel_1(i,1) == 0 || FaceLabel_1(i,2) == 0
        tri_bool2_1(i) = 0;
    end
end

FaceLabel_2 = FaceLabel_raw_2(tri_bool_2,:);
tri_bool2_2 = ones(length(FaceLabel_2),1);
for i = 1: length(FaceLabel_2)
    if FaceLabel_2(i,1) == 0 || FaceLabel_2(i,2) == 0
        tri_bool2_2(i) = 0;
    end
end

% FaceLabel_3 = FaceLabel_raw_3(tri_bool_3,:);
% tri_bool2_3 = ones(length(FaceLabel_3),1);
% for i = 1: length(FaceLabel_3)
%     if FaceLabel_3(i,1) == 0 || FaceLabel_3(i,2) == 0
%         tri_bool2_3(i) = 0;
%     end
% end
% 
% FaceLabel_4 = FaceLabel_raw_4(tri_bool_4,:);
% tri_bool2_4 = ones(length(FaceLabel_4),1);
% for i = 1: length(FaceLabel_4)
%     if FaceLabel_4(i,1) == 0 || FaceLabel_4(i,2) == 0
%         tri_bool2_4(i) = 0;
%     end
% end

% tri_bool2 = logical([tri_bool2_1; tri_bool2_2; tri_bool2_3; tri_bool2_4]);
tri_bool2 = logical([tri_bool2_1; tri_bool2_2]);

% clear FaceLabel_raw_1 FaceLabel_raw_2 FaceLabel_raw_3 FaceLabel_raw_4 tri_bool_1 tri_bool_2 ...
%     FaceLabel_1 FaceLabel_2 FaceLabel_3 FaceLabel_4 tri_bool_3 tri_bool_4 tri_bool2_1 tri_bool2_2 tri_bool2_3 tri_bool2_4;

toPrint = toPrint_raw(tri_bool2,:);

%% -------------------------------------------------------------------
% Write the files
fileID = fopen('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470full_tillSlice29_GBCD.txt','w');

fprintf(fileID,'# Triangles Produced from DREAM3D version 6.4\n');
fprintf(fileID,'# Column 1-3:    right hand average orientation (phi1, PHI, phi2 in RADIANS)\n');
fprintf(fileID,'# Column 4-6:    left hand average orientation (phi1, PHI, phi2 in RADIANS)\n');
fprintf(fileID,'# Column 7-9:    triangle normal\n');
fprintf(fileID,'# Column 10:      surface area\n');
fprintf(fileID,'# Column 11:      triangle curvature\n');

format = '%7.4f    %7.4f    %7.4f    %7.4f    %7.4f  %7.4f  %7.4f  %7.4f  %7.4f   %7.4f   %7.4f\n';
for i = 1:length(toPrint)
%         fprintf(fileID,format,EA_A1(i),EA_A2(i),EA_A3(i),EA_B1(i),EA_B2(i),EA_B3(i), ...
%             n_1(i),n_2(i),n_3(i),area(i), curv(i));
    fprintf(fileID, format, toPrint(i,1), toPrint(i,2), toPrint(i,3), toPrint(i,4), toPrint(i,5), toPrint(i,6), ...
        toPrint(i,7), toPrint(i,8), toPrint(i,9), toPrint(i,10), toPrint(i,11));
end
fclose('all');
