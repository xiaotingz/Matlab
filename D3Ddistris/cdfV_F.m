% file = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub1_GBCD.dream3d');
% num_of_neigh_1 = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors')).';
% surfGrain_1 = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures')).';
% volumes_1 = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumElements')).';
% 
% num_of_neigh = [num_of_neigh_1; num_of_neigh_2];
% surfGrain = [surfGrain_1; surfGrain_2];
% volumes = [volumes_1; volumes_2];

innerGrain = ~logical(surfGrain);
toPlot = [num_of_neigh(innerGrain), volumes(innerGrain)];
toPlot = sortrows(toPlot);

numFaces = unique(num_of_neigh);
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

plot(data_grid(:,1),data_grid(:,3),'k','LineWidth',2)
line([17, 17], [0,1], 'LineStyle','--', 'color', [0.5, 0.5, 0.5])
xlabel('F','FontSize',21);
ylabel('Volume Fraction','FontSize',21);
set(gca,'fontsize',19)

% scatter(data_grid(:,1), data_grid(:,2));
% set(gca, 'YScale', 'log')