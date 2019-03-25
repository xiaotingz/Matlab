% %  -------------------------- bounding boxes -------------------------- 
% % Austenite --- X=434*0.15; Y=267*0.15; Z=100*0.2;
% % Ferrite   --- X=234*0.15; Y=267*0.15; Z=68*0.2;
% % Run4 --- X=201*0.15; Y=335*0.15; Z=108*0.2;
% % STO_1470_sub1 --- X1=232*0.3; Y1=129*0.3; Z1=36*0.3;
% % STO_1470_sub2 --- X2=213*0.3; Y2=297*0.3; Z2=40*0.3;
% % Mg_undeform --- X=167*0.15; Y=133*0.15; Z=73*0.15;
% % Mg_small --- X=167*0.15; Y=133*0.15; Z=17*0.15;
% % MNK_Ti --- X=431*0.5; Y=109*0.5; Z=104*0.3;
% % STO_1350_0729 --- X=375*0.4; Y=250*0.4; Z=28*0.4;
% % STO_1350_0730 --- X=375*0.4; Y=250*0.4; Z=57*0.4;
% % STO_1350_0801 --- X=375*0.4; Y=250*0.4; Z=21*0.4;
% % STO_1350_0802 --- X=250*0.4; Y=250*0.4; Z=51*0.4;
%  --------------------------------------------------------------------
clear
clc 
% ################### load data, one volumes ###################
% file = ('/Users/xiaotingzhong/Desktop/Datas/180704_steelDataResegment/180704_A2_noCentroidAlign_noResegment.dream3d');
% file = '/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/Jan31_Aa0_gbcd10.dream3d';
% file = ('/Users/xiaotingzhong/Desktop/cdf_V/NiAn4_An4new6_fixedOrigin.dream3d');
% % ---- criterion = 'centroidPos' | 'touchingFS' | 'numFaces' | 'NN_centoridPos' | 'NN_touchingFS' ----
% criterion = 'none';
% X=1000; Y=2670; Z=10000;
% % zeroG_F = 17;
% keyWord = 'Ni An4 - ';
% 
% % ----- V6 data -----
% num_of_neigh = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
% volumes = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumElements')).';
% num_of_neigh(1) = []; volumes(1) = [];
% % ----- V4 data -----
% num_of_neigh = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
% volumes = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NumCells'));
% num_of_neigh(1) = []; volumes(1) = [];
% 
% grain_ForCal = filterGrains(criterion, file, X,Y,Z);
% num_of_neigh = num_of_neigh(grain_ForCal);
% volumes = volumes(grain_ForCal);
% %%

% ################### load data, two volumes ###################
%  ----- V4 data -----
% file_1 = ('/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/Jan31_Fca0_gbcd10.dream3d');
% file_2 = ('/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/Jan31_F_Run4_GBCD10.dream3d');
% % ---- criterion = 'centroidPos' | 'touchingFS' | 'numFaces' | 'NN_centoridPos' | 'NN_touchingFS' ----
% criterion = 'none';
% % X1=234*0.15; Y1=267*0.15; Z1=68*0.2;
% % X2=201*0.15; Y2=335*0.15; Z2=108*0.2;
% zeroG_F = 17;
% figureName = ['Jan31_Fca0_gbcd10_comb_Jan31_F_Run4_GBCD10', criterion];
% 
% num_of_neigh_1 = double(h5read(file_1,'/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
% volumes_1 = double(h5read(file_1,'/VoxelDataContainer/FIELD_DATA/NumCells'));
% num_of_neigh_2 = double(h5read(file_2,'/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
% volumes_2 = double(h5read(file_2,'/VoxelDataContainer/FIELD_DATA/NumCells'));
% num_of_neigh_1(1) = []; volumes_1(1) = [];
% num_of_neigh_2(1) = []; volumes_2(1) = [];
% ----- V6 data -----
% file_1 = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub1_GBCD_originCorrect.dream3d');
% file_2 = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub2_GBCD_originCorrect.dream3d');
file_1 = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311_STO1470sub1_GBCD_originCorrect.dream3d');
file_2 = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311_STO1470sub2_GBCD_originCorrect.dream3d');



X1=232*0.3; Y1=129*0.3; Z1=36*0.3;
X2=213*0.3; Y2=297*0.3; Z2=40*0.3;
% % ---- criterion = 'centroidPos' | 'touchingFS' | 'numFaces' | 'NN_centoridPos' | 'NN_touchingFS' ----
criterion = 'no';
keyWord = 'STO1470 combined - ';
% zeroG_F = 16;

num_of_neigh_1 = double(h5read(file_1,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
volumes_1 = double(h5read(file_1,'/DataContainers/ImageDataContainer/CellFeatureData/NumElements')).';
num_of_neigh_2 = double(h5read(file_2,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
volumes_2 = double(h5read(file_2,'/DataContainers/ImageDataContainer/CellFeatureData/NumElements')).';
num_of_neigh_1(1) = []; volumes_1(1) = [];
num_of_neigh_2(1) = []; volumes_2(1) = [];
% ------------------

grain_ForCal_1 = filterGrains(criterion, file_1, X1,Y1,Z1);
grain_ForCal_2 = filterGrains(criterion, file_2, X2,Y2,Z2);
num_of_neigh = [num_of_neigh_1(grain_ForCal_1); num_of_neigh_2(grain_ForCal_2)];
volumes = [volumes_1(grain_ForCal_1); volumes_2(grain_ForCal_2)];


% ##### calculation #####
toPlot = [num_of_neigh, volumes];


toPlot = sortrows(toPlot);

numFaces = unique(toPlot(:,1));
data_grid = zeros(length(numFaces),2);
for i = 1:length(numFaces)
    mask = (toPlot(:,1) == numFaces(i));
    candidates = toPlot(mask,:);
    data_grid(i,1) = numFaces(i);
    data_grid(i,2) = sum(candidates(:,2));
    if i == 1
        data_grid(i,3) = data_grid(i,2);
    else
        data_grid(i,3) = data_grid(i,2) + data_grid(i-1, 3);
    end
end
data_grid(:,3) = data_grid(:,3)/sum(data_grid(:,2));


%%
plot(data_grid(:,1),data_grid(:,3),'k','LineWidth',2)
% line([zeroG_F, zeroG_F], [0,1], 'LineStyle','--', 'color', [0.5, 0.5, 0.5])
line([0, 100], [0.5,0.5], 'LineStyle','--', 'color', [0.5, 0.5, 0.5])
xlabel('F','FontSize',21);
ylabel('Cumulative Volume Fraction','FontSize',21);
set(gca,'fontsize',19)

[~, idx] = min(abs(data_grid(:,3)- 0.5));
text(max(data_grid(:,1))-20, 0.53, ['cdfV=0.5 at ', num2str(data_grid(idx, 1))],'FontSize',18);
xlim([0, 50])
title([keyWord, criterion]);
figureName = [keyWord, criterion];

print(figureName, '-dpng','-r300')

% scatter(data_grid(:,1), data_grid(:,2));
% set(gca, 'YScale', 'log')