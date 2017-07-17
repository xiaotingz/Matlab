clear
file = '/Users/xiaotingzhong/Desktop/STO_1470/Modified/subset2/subset2_curv_A0.dream3d';

facelabel = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels')).';
% curvature_of_triangle = abs(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'));
areas = roundn(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-5);
Size = h5read(file,'/VoxelDataContainer/FIELD_DATA/NumCells');

Size(1,:) = [];
aveSize = 775;

cnt_total = 0;
for i = 1:length(facelabel)
    if facelabel(i,1) <= 0 || facelabel(i,2) <=0
        areas(i) = 0;
        cnt_total = cnt_total + 1;
    end
end
total1 = sum(areas);
cnt_total = length(facelabel) - cnt_total;

cnt_BB = 0;
areas_filtered_BB = zeros(size(areas));
for i = 1:length(facelabel)
    if facelabel(i,1) > 0 && facelabel(i,2) > 0
        if Size(facelabel(i,1)) > aveSize && Size(facelabel(i,2)) > aveSize
            areas_filtered_BB(i) = areas(i);
            cnt_BB = cnt_BB + 1;
        end
    end
end
BB1 = sum(areas_filtered_BB);

cnt_BS = 0;
areas_filtered_BS = zeros(size(areas));
for i = 1:length(facelabel)
    if facelabel(i,1) > 0 && facelabel(i,2) > 0
        if Size(facelabel(i,1)) > aveSize && Size(facelabel(i,2)) <= aveSize
            areas_filtered_BS(i) = areas(i);
            cnt_BS = cnt_BS + 1;
        elseif Size(facelabel(i,1)) <= aveSize && Size(facelabel(i,2)) > aveSize
            areas_filtered_BS(i) = areas(i);
            cnt_BS = cnt_BS + 1;
        end
    end
end
BS1 = sum(areas_filtered_BS);

cnt_SS = 0;
areas_filtered_SS = zeros(size(areas));
for i = 1:length(facelabel)
    if facelabel(i,1) > 0 && facelabel(i,2) > 0
        if Size(facelabel(i,1)) <= aveSize && Size(facelabel(i,2)) <= aveSize
            areas_filtered_SS(i) = areas(i);
            cnt_SS = cnt_SS + 1;
        end
    end
end
SS1 = sum(areas_filtered_SS);

cnt_total   
cnt_BB
cnt_BS
cnt_SS
% 
total1_1 = sum(areas_filtered_BB) + sum(areas_filtered_BS) + sum(areas_filtered_SS)

data_BB2 = h5read('/Users/xiaotingzhong/Desktop/STO_1470/Modified/BB_BS_SS_GBCurvD_area0/subset2_BB_A0.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas');
data_BS2 = h5read('/Users/xiaotingzhong/Desktop/STO_1470/Modified/BB_BS_SS_GBCurvD_area0/subset2_BS_A0.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas');
data_SS2 = h5read('/Users/xiaotingzhong/Desktop/STO_1470/Modified/BB_BS_SS_GBCurvD_area0/subset2_SS_A0.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas');
BB2 = sum(data_BB2)
BS2 = sum(data_BS2)
SS2 = sum(data_SS2)

total = sum(BB2) + sum(BS2) + sum(SS2)
% 
% %% Individual check
% clear
% file = '/Users/xiaotingzhong/Desktop/STO_1470/Modified/BB_BS_SS_GBCurvD/subset2_curv_BB.dream3d';
% label =  h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels').';
% tri_curv = roundn(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'),-5);
% Size = h5read(file,'/VoxelDataContainer/FIELD_DATA/NumCells');
% Size(1,:) = [];

