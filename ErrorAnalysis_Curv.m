%% distribution data in gmt file: mean value and SD

Aug_r1 = textread('Aug_v1_CurvDistri_gmt_1.dat');
Aug_r2 = textread('Aug27_v2_CurvDistri_gmt_1.dat');
% Aug_s4 = textread('Aug24_s4_GBCD_gmt_1.dat');

DistriMean_r1 = sum(Aug_r1(:,3))/length(Aug_r1);
DistriMean_r2 = sum(Aug_r2(:,3))/length(Aug_r2);
% DistriMean_s4 = sum(Aug_s4(:,3))/length(Aug_s4);

total_r1 = 0;
for i = 1 : length(Aug_r1)
    total_r1 = total_r1 + (DistriMean_r1 - Aug_r1(i,3))^2;
end
DistriSD_r1 = sqrt(total_r1/length(Aug_r1));

total_r2 = 0;
for i = 1 : length(Aug_r2)
    total_r2 = total_r2 + (DistriMean_r2 - Aug_r2(i,3))^2;
end
DistriSD_r2 = sqrt(total_r2/length(Aug_r2));
% 
% total_s4 = 0;
% for i = 1 : length(Aug_s4)
%     total_s4 = total_s4 + (DistriMean_s4 - Aug_s4(i,3))^2;
% end
% DistriSD_s4 = sqrt(total_s4/length(Aug_s4));


%% GBCD/CurvDistri data in DREAM3D file: mean value and SD

CurvDistri = h5read('/Volumes/RESEARCH/Simple Geometry Jun.21_Aug/Aug.27 v1_v4abs/Aug_v1_GBCD.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');

j = 1;
for i = 1 : length(CurvDistri)
    if CurvDistri(i,1) ~= 0 || CurvDistri(i,2) ~= 0;
        data(j,1) = i;
        data(j,2) = CurvDistri(i,1);
        data(j,3) = CurvDistri(i,2);
        j = j + 1;
    end
end

mean = sum(data(:,2))/length(data);
total = 0;
for i = 1 : length(data)
    total = total + (mean - data(i,2))^2;
end
SD = sqrt(total/length(data));

MaxValue = max(data(:,2));
MinValue = min(data(:,2));



%% curvature data in DREAM3D file
%     mean1 is using #Triangles
%     mean2 is using TriangleAreas

facelabel = double(h5read('/Volumes/RESEARCH/Simple Geometry Jun.21_Aug/Aug.24 r1_v4abs/Aug24_r1_ForCurvDistri.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
curvature_of_triangle = abs(h5read('/Volumes/RESEARCH/Simple Geometry Jun.21_Aug/Aug.24 r1_v4abs/Aug24_r1_ForCurvDistri.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'));
triangle_area = roundn(h5read('/Volumes/RESEARCH/Simple Geometry Jun.21_Aug/Aug.24 r1_v4abs/Aug24_r1_ForCurvDistri.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-8);

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

% average according to #Triangles
mean1 = sum(data_cleared(:,3))/length(data_cleared(:,3));

% average according to TriArea
total = 0;
for i = 1:length(data_cleared)
    total = total + data_cleared(i,3) * data_cleared(i,4);
end
mean2 = total/sum(data_cleared(:,4));

% SD according to TriArea
total = 0;
for i = 1:length(data_cleared)
    total = total + (mean2 - data_cleared(i,3))^2 ;
end
SD = sqrt(total/length(data_cleared));