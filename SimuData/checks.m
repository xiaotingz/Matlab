file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_simu/simu_An4_clean_seg.dream3d';
file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_simu/simu_An5_clean_seg.dream3d';

% file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin3_Hsmooth.dream3d';
% file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_cropToAn4_Hsmooth.dream3d';

%%
% fl_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
% fl_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
% 
% fl_an4 = fl_an4(all(fl_an4>0, 2), :);
% fl_an5 = fl_an5(all(fl_an5>0, 2), :);

% rng('shuffle')
% idx = randi(size(tracked_uniqueface_an4, 1));
idx = 977;
disp('------------------------------- ')
dispFacePairInfo(file_an4, file_an5, tracked_uniqueface_an4, tracked_uniqueface_an5, idx)
disp(['area_an4 = ', num2str(face_area_an4(idx)), ';    area_diff = ', num2str(face_area_diff(idx))])
disp('  ')


data = [(1:length(face_area_an4))', face_area_an4, face_area_diff];
data = sortrows(data, -3);


%%
fl_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
fl_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
area_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceData/FaceAreas')';
area_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceData/FaceAreas')';

mask = all(fl_an4>0, 2);
fl_an4 = fl_an4(mask, :);
fl_an4 = sort(fl_an4, 2);
area_an4 = area_an4(mask);
mask = all(fl_an5>0, 2);
fl_an5 = fl_an5(mask, :);
fl_an5 = sort(fl_an5, 2);
area_an5 = area_an5(mask);

% load('/Users/xiaotingzhong/Desktop/Datas/Ni_simu/simu_features_geo_topo.mat', ...
%     'tracked_uniqueface_an4', 'tracked_uniqueface_an5', 'face_area_an4');
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190425_Hsmooth_geo_topo_an5crop2.mat')
total_area_an4 = sum(area_an4);
tracked_uniqueface_an4 = sort(tracked_uniqueface_an4, 2);
mask_tracked = ismember(fl_an4, tracked_uniqueface_an4, 'rows');

sum(area_an4(mask_tracked)) / total_area_an4 * 100
sum(face_area_an4) / total_area_an4 * 100

%%
num_neigh_an4 = h5read(file_an4,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors').';
num_neigh_an5 = h5read(file_an5,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors').';
size_an4 = roundn(h5read(file_an4,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters'),-5).';
size_an5 = roundn(h5read(file_an5,'/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters'),-5).';

num_neigh_an4(1) = [];
num_neigh_an5(1) = [];
size_an4(1) = [];
size_an5(1) = [];

num_neighs = [num_neigh_an4(corresp_d3d_45(:, 1)), num_neigh_an5(corresp_d3d_45(:, 2))];
sizes = [size_an4(corresp_d3d_45(:, 1)), size_an5(corresp_d3d_45(:, 2))];
F_mF_diff_tracked = F_mF_diff(corresp_d3d_45(:, 1));

fileID = fopen('simu_grains_data.txt','w');
fprintf(fileID,'%s, %s, %s, %s\n', ...
    'size_diff', 'size_an4', 'num_faces_an4', 'F_mF');
for i = 1:length(sizes)
    fprintf(fileID, '%6.3f, %6.3f, %6.3f, %6.3f\n', ...
        sizes(i, 2) - sizes(i, 1), sizes(i, 1), num_neighs(i, 1), F_mF_diff_tracked(i));
end
fclose(fileID);





%%
set(0,'defaultAxesFontSize',18)
histogram(da_len_w_an4(:,3), 'FaceColor', [0.5, 0.5, 0.5])
ylabel('simu\_dihedral\_angle\_opposite\_an4')
print('opp_dihedral_ang','-dtiff','-r300')

figure()
histogram(da_len_w_an5(:,3) - da_len_w_an4(:,3), 'FaceColor', [0.5, 0.5, 0.5])
ylabel('simu\_dihedral\_angle\_opposite\_diff')
print('opp_dihedral_ang_diff','-dtiff','-r300')






