z100_mimic = '/Users/xiaotingzhong/Desktop/Datas/STO_1425/041917_V4_mimic/z100_mimic_stats.dream3d';
z200_mimic = '/Users/xiaotingzhong/Desktop/Datas/STO_1425/041917_V4_mimic/z200_mimic_stats.dream3d';
z100_Mandy = '/Users/xiaotingzhong/Desktop/Datas/STO_1425/Medeleine/Final_MNKelly_selected/STO_1425C_50nmslices_100nmEBSDseparation_featuresaftermesh.dream3d';
z200_Mandy = '/Users/xiaotingzhong/Desktop/Datas/STO_1425/Medeleine/Final_MNKelly_selected/STO_1425_100nmslice200nmEBSD_meshcheck_featureanalysis.dream3d';
z100_default = '/Users/xiaotingzhong/Desktop/Datas/STO_1425/041917_V6_default/041917_v6default_GBCD.dream3d';

%% Triangle stats

z100_mimic_label = double(h5read(z100_mimic, '/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels')).';
z200_mimic_label = double(h5read(z200_mimic, '/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels')).';
% z100_Mandy_label = double(h5read(z100_Mandy, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
% z200_Mandy_label = double(h5read(z200_Mandy, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
z100_mimic_tri = h5read(z100_mimic, '/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas');
z200_mimic_tri = h5read(z200_mimic, '/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas');
% z100_Mandy_tri = h5read(z100_Mandy, '/DataContainers/TriangleDataContainer/FaceData/FaceAreas');
% z200_Mandy_tri = h5read(z200_Mandy, '/DataContainers/TriangleDataContainer/FaceData/FaceAreas');


z100_mimic_surfTriCnt = 0;
z100_mimic_surfTriSum = 0;
for i = 1:length(z100_mimic_label)
    if z100_mimic_label(i,1) <= 0 || z100_mimic_label(i,2) <=0
        z100_mimic_surfTriSum = z100_mimic_surfTriSum + z100_mimic_tri(i);
        z100_mimic_surfTriCnt = z100_mimic_surfTriCnt + 1;
    end
end

z200_mimic_surfTriCnt = 0;
z200_mimic_surfTriSum = 0;
for i = 1:length(z200_mimic_label)
    if z200_mimic_label(i,1) <= 0 || z200_mimic_label(i,2) <=0
        z200_mimic_surfTriSum = z200_mimic_surfTriSum + z200_mimic_tri(i);
        z200_mimic_surfTriCnt = z200_mimic_surfTriCnt + 1;
    end
end
% 
% z100_Mandy_surfTriCnt = 0;
% z100_Mandy_surfTriSum = 0;
% for i = 1:length(z100_Mandy_label)
%     if z100_Mandy_label(i,1) <= 0 || z100_Mandy_label(i,2) <=0
%         z100_Mandy_surfTriSum = z100_Mandy_surfTriSum + z100_Mandy_tri(i);
%         z100_Mandy_surfTriCnt = z100_Mandy_surfTriCnt + 1;
%     end
% end
% 
% z200_Mandy_surfTriCnt = 0;
% z200_Mandy_surfTriSum = 0;
% for i = 1:length(z200_Mandy_label)
%     if z200_Mandy_label(i,1) <= 0 || z200_Mandy_label(i,2) <=0
%         z200_Mandy_surfTriSum = z200_Mandy_surfTriSum + z200_Mandy_tri(i);
%         z200_Mandy_surfTriCnt = z200_Mandy_surfTriCnt + 1;
%     end
% end


z100_mimic_numTris = length(z100_mimic_tri)
% z100_Mandy_numTris = length(z100_Mandy_tri)
z100_mimic_surfTriCnt
% z100_Mandy_surfTriCnt
z100_mimic_TriSum = sum(z100_mimic_tri)
% z100_Mandy_TriSum = sum(z100_Mandy_tri)
z100_mimic_surfTriSum 
% z100_Mandy_surfTriSum 

z200_mimic_numTris = length(z200_mimic_tri)
% z200_Mandy_numTris = length(z200_Mandy_tri)
z200_mimic_surfTriCnt
% z200_Mandy_surfTriCnt
z200_mimic_TriSum = sum(z200_mimic_tri)
% z200_Mandy_TriSum = sum(z200_Mandy_tri)
z200_mimic_surfTriSum
% z200_Mandy_surfTriSum





%% Cell stats
% z100_mimic_cells = double(h5read(z100_mimic,'/VoxelDataContainer/FIELD_DATA/NumCells'));
% z200_mimic_cells = double(h5read(z200_mimic,'/VoxelDataContainer/FIELD_DATA/NumCells'));
z100_Mandy_cells = double(h5read(z100_Mandy,'/DataContainers/ImageDataContainer/CellFeatureData/NumCells'));
z200_Mandy_cells = double(h5read(z200_Mandy,'/DataContainers/ImageDataContainer/CellFeatureData/NumCells'));
% z100_default_cells = double(h5read(z100_default,'/DataContainers/ImageDataContainer/CellFeatureData/NumCells'));


% z100_mimic_sumCell = sum(z100_mimic_cells)
z100_Mandy_sumCell = sum(z100_Mandy_cells)

% z200_mimic_sumCell = sum(z200_mimic_cells)
z200_Mandy_sumCell = sum(z200_Mandy_cells)

% z100_default_sumCell = sum(z100_default_cells)

% z100_mimic_numGrains = length(z100_mimic_cells) - 1
z100_Mandy_numGrains = length(z100_Mandy_cells) - 1
% z200_mimic_numGrains = length(z200_mimic_cells) - 1
z200_Mandy_numGrains = length(z200_Mandy_cells) - 1


% z100_mimic_sumCell/z100_Mandy_sumCell 
% z200_mimic_sumCell/z200_Mandy_sumCell



%% #Faces plot
z100_mimic = '/Users/xiaotingzhong/Desktop/Datas/STO_1425/041917_V4_mimic/z100_mimic_stats_a0_D.dream3d';
z200_mimic = '/Users/xiaotingzhong/Desktop/Datas/STO_1425/041917_V4_mimic/z200_mimic_stats_a0_D.dream3d';
z100_mimic_Faces = double(h5read(z100_mimic,'/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'));
z200_mimic_Faces = double(h5read(z200_mimic,'/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'));
z100_mimic_Faces(1) = [];
z200_mimic_Faces(1) = [];

figure(1)
histogram(z100_mimic_Faces)
title('z100 #Faces','FontSize',21)
xlabel('#Faces','FontSize',21);
ylabel('#Grains','FontSize',21);
ax = gca;
set(ax,'fontsize',19)

figure(2)
histogram(z200_mimic_Faces)
title('z200 #Faces','FontSize',21)
xlabel('#Faces','FontSize',21);
ylabel('#Grains','FontSize',21);
ax = gca;
set(ax,'fontsize',19)