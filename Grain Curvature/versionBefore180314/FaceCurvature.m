% important data:
    % data_final [LabelA, LabelB, FaceArea, FaceCurvature]
clear
% subset1
% file = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/070617_V4_misA3ReconsAgrain/070617_sub1_stats.dream3d');
% centro_file = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/070617_V4_misA3ReconsAgrain/070617_sub1_stats.dream3d');
% X=232*0.3; Y=129*0.3; Z=36*0.3;
% subset2
file = ('/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/Jan31_Ac_setTo0.dream3d');
centro_file = ('/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/Jan31_Ac_setTo0.dream3d');
X=434*0.15; Y=267*0.15; Z=100*0.2;

% centroids
% Austenite --- X=434*0.15; Y=267*0.15; Z=100*0.2;
% Ferrite   --- X=234*0.15; Y=267*0.15; Z=68*0.2;
% STO_1470_sub1 --- X=232*0.3; Y=129*0.3; Z=36*0.3;
% STO_1470_sub2 --- X=213*0.3; Y=297*0.3; Z=40*0.3;
% Mg_undeform --- X=167*0.15; Y=133*0.15; Z=73*0.15;
% Mg_small --- X=167*0.15; Y=133*0.15; Z=17*0.15;
% MNK_Ti --- X=431*0.5; Y=109*0.5; Z=104*0.3;

% load data_v4
centroids = roundn(h5read(centro_file,'/VoxelDataContainer/FIELD_DATA/Centroids'),-5).';
grain_diameter_raw = roundn(h5read(centro_file,'/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
facelabel = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
curvature_of_triangle = h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures');
num_of_neigh = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
neighborList = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NeighborList'));
triangle_area_raw = roundn(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-8);
data_raw = [facelabel; curvature_of_triangle.'; triangle_area_raw.'];

% load data_v6
% centroids = abs(roundn(h5read(centro_file,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).');
% grain_diameter_raw = roundn(h5read(centro_file,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters'),-5);
% num_of_neigh = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors'));
% neighborList = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList'));
% facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'));
% curvature_of_triangle = h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures');
% triangle_area_raw = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-8);
% data_raw = [facelabel; curvature_of_triangle; triangle_area_raw];



tmp1 = data_raw(1,:);
tmp2 = data_raw(2,:);
tmp3 = data_raw(3,:);
tmp4 = data_raw(4,:);

% get rid of bad datas: 
    % 1.facelabel <= 0; 2.extreme curvature value; 3.NaN(sometimes)
boolean1 = data_raw(1,:) > 0 & data_raw(2,:) > 0 & ~isnan(data_raw(3,:)) & data_raw(3,:) < 100 & data_raw(3,:) > -100;
data_cleared(1,:) = tmp1(boolean1);
data_cleared(2,:) = tmp2(boolean1);
data_cleared(3,:) = tmp3(boolean1);
data_cleared(4,:) = tmp4(boolean1);

% prepare for face curvature calculation: sort the data by facelabel
data_sorted = sortrows(data_cleared.');
clear tmp1 tmp2 tmp3 tmp4 boolean1 data_raw data_cleared

% calculate face curvature
tmp5 = unique(data_sorted(:,1:2),'rows');
data_final = zeros(length(tmp5),4);
j = 0;
k = 1;
counter = 2;
s_area = data_sorted(j+counter-1,4);
s_cur = data_sorted(j+counter-1,3)*data_sorted(j+counter-1,4);

for i = 2:length(data_sorted)
    
%   next trianlge still in the same face
    if data_sorted(j+counter-1,1) == data_sorted(j+counter,1) && data_sorted(j+counter-1,2) == data_sorted(j+counter,2)
       s_cur = s_cur + data_sorted(j+counter,3)*data_sorted(j+counter,4);
       s_area = s_area + data_sorted(j+counter,4);
       counter = counter + 1;
       
%   next trianlge in another face: record former one and restart counting
    else
        data_final(k,1) = data_sorted(j+counter-1,1);
        data_final(k,2) = data_sorted(j+counter-1,2);
        data_final(k,3) = s_area;
        data_final(k,4) = s_cur/s_area;
        k = k+1;   % check k = length(x)
        j = j + counter -1;
        counter = 2;
        s_area = data_sorted(j+counter-1,4);
        s_cur = data_sorted(j+counter-1,3)*data_sorted(j+counter-1,4);
    end
end
% looped all but no comparison to execute else for the last face ==> write manually
data_final(k,1) = data_sorted(j+counter-1,1);
data_final(k,2) = data_sorted(j+counter-1,2);
data_final(k,3) = s_area;
data_final(k,4) = s_cur/s_area;

% % figure('name','closer look') %change axis value to look closer
% % scatter(data_final(:,3), data_final(:,4),2,'filled');
% % xlabel('face area','FontSize',13);
% % ylabel('face curvature','FontSize',13);
% % axis([0 10 -5 5]);
