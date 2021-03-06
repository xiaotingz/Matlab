%% ##### get position of a SixCVoxels #####
% """
% The ID calced for SixCNodes can be used directly for SixCVoxels. No need
% to minus 1 because both of them starts at 1.
% """

load('180625_Ni_TopologyResult.mat');
clearvars -except sixCoordNList_An4 sixCoordNList_An5 fiveCoordNList_An4 fiveCoordNList_An5

sixCoordNList_An4 = sixCoordNList_An4(~any(sixCoordNList_An4==0, 2), :);

file_An4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6.dream3d';
CI_An4 = squeeze(double(h5read(file_An4, '/DataContainers/ImageDataContainer/CellData/Confidence Index')));
dims_An4 = double(h5read(file_An4, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS'));

file_An5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6.dream3d';
CI_An5 = squeeze(double(h5read(file_An5, '/DataContainers/ImageDataContainer/CellData/Confidence Index')));
dims_An5 = double(h5read(file_An5, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS'));

%%
SixCV_CI_An4 = zeros(length(sixCoordNList_An4), 1);
SixCV_CI_An5 = zeros(length(sixCoordNList_An5), 1);

[X, Y, Z] = ind2sub(dims_An4 + 1, sixCoordNList_An4(:,1));
SixCV_pos_An4 = [X, Y, Z];
for i = 1:length(sixCoordNList_An4)
    SixCV_CI_An4(i) = CI_An4(SixCV_pos_An4(i,1), SixCV_pos_An4(i,2), SixCV_pos_An4(i,3));
end
[X, Y, Z] = ind2sub(dims_An5 + 1, sixCoordNList_An5(:,1));
SixCV_pos_An5 = [X, Y, Z];
for i = 1:length(sixCoordNList_An5)
    SixCV_CI_An5(i) = CI_An5(SixCV_pos_An5(i,1), SixCV_pos_An5(i,2), SixCV_pos_An5(i,3));
end

fiveCoordNList_An4 = fiveCoordNList_An4(~any(fiveCoordNList_An4==0, 2), :);
fiveCoordNList_An5 = fiveCoordNList_An5(~any(fiveCoordNList_An5==0, 2), :);

FiveCV_CI_An4 = zeros(length(fiveCoordNList_An4), 1);
FiveCV_CI_An5 = zeros(length(fiveCoordNList_An5), 1);

[X, Y, Z] = ind2sub(dims_An4 + 1, fiveCoordNList_An4(:,1));
FiveCV_pos_An4 = [X, Y, Z];
for i = 1:length(fiveCoordNList_An4)
    FiveCV_CI_An4(i) = CI_An4(FiveCV_pos_An4(i,1), FiveCV_pos_An4(i,2), FiveCV_pos_An4(i,3));
end

[X, Y, Z] = ind2sub(dims_An5 + 1, fiveCoordNList_An5(:,1));
FiveCV_pos_An5 = [X, Y, Z];
for i = 1:length(fiveCoordNList_An5)
    FiveCV_CI_An5(i) = CI_An5(FiveCV_pos_An5(i,1), FiveCV_pos_An5(i,2), FiveCV_pos_An5(i,3));
end

%%
histogram(FiveCV_CI_An4)
hold on
histogram(FiveCV_CI_An5)
legend('An4','An5','Location','NW');
xlabel('Confidence of FCNs');
ylabel('Counts');
print('FCN_CIs', '-dpng','-r300')

figure()
histogram(SixCV_CI_An4)
hold on
histogram(SixCV_CI_An5)
legend('An4','An5','Location','NW');
xlabel('Confidence of SixCNs');
ylabel('Counts');
print('SixCN_CIs', '-dpng','-r300')


%% ##### Plot CI of SixCNs #####

CI_An4_flat = reshape(CI_An4,[],1);
CI_An4_flat(CI_An4_flat==0) = [];
CI_An5_flat = reshape(CI_An5,[],1);
CI_An5_flat(CI_An5_flat==0) = [];


cdfplot(CI_An4_flat);
hold on 
cdfplot(CI_An5_flat);
legend('An4','An5','Location','NW');
xlabel('Confidence Index');
ylabel('CDF');
print('CI_cdf', '-dpng','-r300')




%% ##### Change SixCN orientation in a cropped volume #####

file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/SixCN_cropVolume/An4_cropSixCN_11353436.dream3d';
writeFile = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/SixCN_cropVolume/An4_cropSixCN_11353436_changeID.dream3d';
CI = squeeze(double(h5read(file, '/DataContainers/ImageDataContainer/CellData/Confidence Index')));
IDs = double(h5read(file, '/DataContainers/ImageDataContainer/CellData/FeatureIds'));
size = double(h5read(file, '/DataContainers/ImageDataContainer/CellFeatureData/NumElements'));
% EAs = double(h5read(file, '/DataContainers/ImageDataContainer/CellData/EulerAngles'));
% Quats = double(h5read(file, '/DataContainers/ImageDataContainer/CellData/Quats'));
dims = double(h5read(file, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS'));

[result, QNList, fiveCoordNList, sixCoordNList] = findQuadNodes(file);
[X, Y, Z] = ind2sub(dims + 1, sixCoordNList(1));
maxCI = 0;
maxCI_idx = [0,0,0];
for i = [-1, 1]
    for j = [-1, 1]
        for k = [-1, 1]
            if CI(X+i, Y+j, Z+k) > maxCI
                maxCI = CI(X+i, Y+j, Z+k);
                maxCI_idx = [X+i, Y+j, Z+k];
            end
        end
    end
end
% EAs(:, X, Y, Z) = EAs(:, maxCI_idx(1), maxCI_idx(2), maxCI_idx(3));
IDs(1, X, Y, Z) = IDs(1, maxCI_idx(1), maxCI_idx(2), maxCI_idx(3));
h5write(writeFile, '/DataContainers/ImageDataContainer/CellData/FeatureIds', IDs);


% ----- find Quads in the cropped volume -----

[~, ~, ~, sixCoordNList] = findQuadNodes(file);
[~, ~, ~, sixCoordNList2] = findQuadNodes(writeFile);
ID = squeeze(double(h5read(file, '/DataContainers/ImageDataContainer/CellData/FeatureIds')));
ID2 = squeeze(double(h5read(writeFile, '/DataContainers/ImageDataContainer/CellData/FeatureIds')));

% ----- check size of the host grain of SixCN, before and after change -----
sum(sum(sum(ID ~= ID2)))
size(ID(X, Y, Z))
size(ID2(X, Y, Z))


