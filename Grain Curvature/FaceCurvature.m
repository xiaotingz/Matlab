% important data:
    % data_final [LabelA, LabelB, FaceArea, FaceCurvature]

% load data
facelabel = double(h5read('/Volumes/RESEARCH/Grain Curvature/Jun.24 Ferrite/Jun.24 ForCurv_centroid.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
curvature_of_triangle = h5read('/Volumes/RESEARCH/Grain Curvature/Jun.24 Ferrite/Jun.24 ForCurv_centroid.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures');
num_of_neigh = double(h5read('/Volumes/RESEARCH/Grain Curvature/Jun.24 Ferrite/Jun.24 ForCurv_centroid.dream3d','/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
neighborList = double(h5read('/Volumes/RESEARCH/Grain Curvature/Jun.24 Ferrite/Jun.24 ForCurv_centroid.dream3d','/VoxelDataContainer/FIELD_DATA/NeighborList'));
triangle_area_raw = roundn(h5read('/Volumes/RESEARCH/Grain Curvature/Jun.24 Ferrite/Jun.24 ForCurv_centroid.dream3d','/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-8);

data_raw = [facelabel; curvature_of_triangle.'; triangle_area_raw.'];

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
