% load('180625_Ni_TopologyResult.mat');
% clearvars -except sixCoordNList_An4 sixCoordNList_An5

file_An4 = '/Users/xiaotingzhong/Desktop/data/180510/An4new6_fixedOrigin_mesh.dream3d';
CI_An4 = squeeze(double(h5read(file_An4, '/DataContainers/ImageDataContainer/CellData/Confidence Index')));
dims_An4 = double(h5read(file_An4, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS'));

file_An5 = '/Users/xiaotingzhong/Desktop/data/180510/An5new6_mesh.dream3d';
CI_An5 = squeeze(double(h5read(file_An5, '/DataContainers/ImageDataContainer/CellData/Confidence Index')));
dims_An5 = double(h5read(file_An5, '/DataContainers/ImageDataContainer/_SIMPL_GEOMETRY/DIMENSIONS'));


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

CI_An4_flat = reshape(CI_An4,[],1);
CI_An4_flat(CI_An4_flat==0) = [];
CI_An5_flat = reshape(CI_An5,[],1);
CI_An5_flat(CI_An5_flat==0) = [];


% cdfplot(CI_An4_flat);
% hold on 
% cdfplot(CI_An5_flat);
% legend('An4','An5','Location','NW');
% title('Ni CI distributions');

%%
file = '/Users/xiaotingzhong/Desktop/SixCN_cropVolume/An4_cropSixCN_886184.dream3d';
CI = squeeze(double(h5read(file, '/DataContainers/ImageDataContainer/CellData/Confidence Index')));
EAs = squeeze(double(h5read(file, '/DataContainers/ImageDataContainer/CellData/EulerAngles')));
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
EAs(:,X, Y, Z) = EAs(:,maxCI_idx(1), maxCI_idx(2), maxCI_idx(3));


