% % clear
% % file = '/Users/xiaotingzhong/Desktop/Datas/STO_1470/180311/180311_STO1470sub2_GBCD.dream3d';
file  = ('/Users/xiaotingzhong/Desktop/Datas/STO_1470/V6_MNK/subset1_fullrecon_MNK.dream3d');
% % % % --------------------------- V6 structure ---------------------------
numNeigh_1 = h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors').';
D_1 = h5read(file, '/DataContainers/ImageDataContainer/CellFeatureData/EquivalentDiameters').';
surfGrain_1 = double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/SurfaceFeatures')).';
% % centroids = roundn(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/Centroids'),-5).';
% % D = roundn(h5read(file,'/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
% % numNeigh = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
numNeigh_1(1) = [];
D_1(1) = [];
surfGrain_1(1) = [];
% centroids(1,:) = [];
% 
% mask = ~logical(surfGrain);
% D = D(mask);
% numNeigh = numNeigh(mask);

% Dave = sum(D)/length(D)
% Fave = sum(numNeigh)/length(numNeigh)
%%
scatter(data_grain(:,4),data_grain(:,3),50,'filled','k')
xlabel('E','FontSize',21);
ylabel('F','FontSize',21);
set(gca,'fontsize',19)
box on

%%
% % h = histogram(numNeigh,'Normalization','probability','FaceColor',[0.5,0.5,0.5]);
% h = histogram(log(data_grain(:,2)/sum(data_grain(:,2))*length(data_grain(:,2))),'Normalization','probability','FaceColor',[0.5,0.5,0.5]);
% h.NumBins =  60;
% xlabel('log(D/<D>)','FontSize',21);
% ylabel('Frequency','FontSize',21);
% set(gca,'fontsize',19)


figure(3)
x = D(~logical(surfGrain));
y = double(numNeigh(~logical(surfGrain)));
scatter(x,y,10,'filled','k')
% N = hist3([y,x], 'Nbins',[50,50]);
% colormap(flipud(hot))
% set(gca,'YDir','normal');
% % pcolor(N)
% p = imagesc(N);
% axis('xy')
%%

figure(1)
scatter(data_grain(:,2),data_grain(:,3),10,'filled','k')
figure(2)
x = data_grain(:,2);
y = data_grain(:,3);
xbins = 30; xMin = 1; xMax = 10;
ybins = 30; yMax = 55; yMin = 0;
% N = hist3([y,x], 'edges', {min(y):(max(y)-min(y))/Nbins:max(y) min(x):(max(x)-min(x))/Nbins:max(x)} );
N = hist3([y,x], 'edges', {yMin:(yMax-yMin)/xbins:yMax xMin:(xMax-xMin)/ybins:xMax} );
colormap(flipud(hot))
set(gca,'YDir','normal');
% pcolor(N)
p = imagesc(N);
hold on
for i = 1:size(N,1)
    for j = 1:size(N,2)
        if N(i,j) ~= 0
            line([j-0.5, j-0.5],[i-0.5, i+0.5],'color', [0.8, 0.8, 0.8])
            line([j+0.5, j+0.5],[i-0.5, i+0.5],'color', [0.8, 0.8, 0.8])
            line([j-0.5, j+0.5],[i-0.5, i-0.5],'color', [0.8, 0.8, 0.8])
            line([j-0.5, j+0.5],[i+0.5, i+0.5],'color', [0.8, 0.8, 0.8])
            hold on
        end
    end
end
yTicks = [0 : 10*ybins/55 : ybins*50/55];
xTicks = [0, xbins/9: 2*xbins/9: 50];
yLabels = ({'0', '10', '20', '30', '40', '50'});
xLabels = ({'1', '2', '4', '6', '8', '10'});
set(gca, 'YTick', yTicks, 'XTick', xTicks, 'XTickLabels',xLabels, 'YTickLabels',yLabels,  'XMinorTick','on','YMinorTick','on')
% grid minor
pbaspect([1 1 1])

axis('xy')
xlabel('D (\mum)','FontSize',21);
ylabel('F','FontSize',21);
set(gca,'fontsize',19)
colorbar



% xb = linspace(min(x),max(x),size(N,1));
% yb = linspace(min(y),max(y),size(N,2));
% 
% hist3([x,y],N);
% h = pcolor(xb,yb,N);
% colormap('hot') % Change color scheme 
% colorbar % Display colorbar
% h.ZData = ones(size(N))*-max(max(N));
% ax = gca;
% ax.ZTick(ax.ZTick < 0) = [];

% xlabel('D (\mum)','FontSize',21);
% ylabel('F','FontSize',21);
% set(gca,'fontsize',19)
% box on


% % disable tick on top and right of the box
%     % get handle to current axes
% a = gca;
%     % set box property to off and remove background color
% set(a,'box','off','color','none')
%     % create new, empty axes with box but without ticks
% b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'FontWeight','bold','linewidth',2);
%     % set original axes as active
% axes(a)
%     % link axes in case of zooming
% linkaxes([a b])
