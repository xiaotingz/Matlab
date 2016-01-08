%% distribution data in gmt file: mean value and SD

% d3d = textread('AxisD100_gmt_1.dat');
gbcd_graph = textread('No10_gmt_1.dat');
% Aug_s4 = textread('Aug24_s4_GBCD_gmt_1.dat');

% d3d_mean = sum(d3d(:,3))/length(d3d)
graph_gbcd_mean = sum(gbcd_graph(:,3))/length(gbcd_graph)
% DistriMean_s4 = sum(Aug_s4(:,3))/length(Aug_s4);
% 
% total_d3d = 0;
% for i = 1 : length(d3d)
%     total_d3d = total_d3d + (d3d_mean - d3d(i,3))^2;
% end
% d3d_SD = sqrt(total_d3d/length(d3d))

total_gbcd_graph = 0;
for i = 1 : length(gbcd_graph)
    total_gbcd_graph = total_gbcd_graph + (graph_gbcd_mean - gbcd_graph(i,3))^2;
end
graph_gbcd_SD = sqrt(total_gbcd_graph/length(gbcd_graph))
% 
% total_s4 = 0;
% for i = 1 : length(Aug_s4)
%     total_s4 = total_s4 + (DistriMean_s4 - Aug_s4(i,3))^2;
% end
% DistriSD_s4 = sqrt(total_s4/length(Aug_s4));


%% Character distribution data in DREAM3D file: mean value and SD
clear

% % version4
CurvDistri = h5read('/Volumes/RESEARCH/Oct.26 MinNeigh=4/Oct.26 Ferrite/Oct26_F_CurvDistri.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
% % version6
% CurvDistri = h5read('/Volumes/RESEARCH/Oct.7 No1-10/No1/No1_CurvDistri.dream3d','/DataContainers/TriangleDataContainer/FaceEnsembleData/GBCD');

j = 1;
for i = 1 : length(CurvDistri)
    if CurvDistri(i,1) ~= 0 || CurvDistri(i,2) ~= 0;
        data(j,1) = i;
        data(j,2) = CurvDistri(i,1);
        data(j,3) = CurvDistri(i,2);
        j = j + 1;
    end
end

mean = sum(data(:,2))/length(data)
total = 0;
for i = 1 : length(data)
    total = total + (mean - data(i,2))^2;
end
SD = sqrt(total/length(data))

% MaxValue = max(data(:,2));
% MinValue = min(data(:,2));



%% curvature data in DREAM3D file
%     mean1 is using #Triangles
%     mean2 is using TriangleAreas
clear

file = ('/Volumes/RESEARCH/Simple Geometry/Oct.7 t1-t10/curvature distributions/t9_CurvDistri.dream3d');
facelabel = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
curvature_of_triangle = abs(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'));
triangle_area = roundn(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-5);

data_raw = [facelabel; curvature_of_triangle.'; triangle_area.'];
tmp1 = data_raw(1,:);
tmp2 = data_raw(2,:);
tmp3 = data_raw(3,:);
tmp4 = data_raw(4,:);

% shouldn't exclude extreme curvature here: dream3d didn't when binning.
boolean = data_raw(1,:) > 0 & data_raw(2,:) > 0; %& data_raw(3,:) < 10 & data_raw(3,:) > -10;
data_cleared(1,:) = tmp1(boolean);
data_cleared(2,:) = tmp2(boolean);
data_cleared(3,:) = tmp3(boolean);
data_cleared(4,:) = tmp4(boolean);

data_cleared = data_cleared.';

% % average according to #Triangles
% mean1 = sum(data_cleared(:,3))/length(data_cleared(:,3));

% average according to TriArea
total = 0;
for i = 1:length(data_cleared)
    total = total + data_cleared(i,3) * data_cleared(i,4);
end
mean = total/sum(data_cleared(:,4))

% SD according to TriArea
total = 0;
for i = 1:length(data_cleared)
    total = total + (mean - data_cleared(i,3))^2 ;
end
SD = sqrt(total/length(data_cleared))

length(data_cleared)


%% compare the .dat results
No1 = textread('s4_GBCDRes9_gmt_1.dat');
No10 = textread('s4_GBCDRes9_converted_gmt_1.dat');
No1(1,:) = [];
No10(1,:) = [];

ave_D3D = sum(No1(:,3))/length(No1(:,3));
ave_graph_gbcd = sum(No10(:,3))/length(No10(:,3));
max_D3D = max(No1(:,3));
max_graph_gbcd = max(No10(:,3));

difference = roundn((No1(:,3) - No10(:,3)),-3);
% ave_diff = sum(abs(difference))/length(difference);
max_diff = max(difference);