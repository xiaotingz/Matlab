file = ('/Users/xiaotingzhong/Desktop/072416 validation/setTo0/A_stats_RingCount2_setTo0.dream3d');
% load data
facelabel = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
curvature_of_triangle = abs(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'));

data_raw = [facelabel; curvature_of_triangle.'];

tmp1 = data_raw(1,:);
tmp2 = data_raw(2,:);
tmp3 = data_raw(3,:);

    % get rid of bad datas: 
        % 1.facelabel <= 0; 2.extreme curvature value; 3.NaN(sometimes)
boolean1 = data_raw(1,:) > 0 & data_raw(2,:) > 0 & ~isnan(data_raw(3,:)) & data_raw(3,:) < 100 & data_raw(3,:) > -100;
data_cleared(1,:) = tmp1(boolean1);
data_cleared(2,:) = tmp2(boolean1);
data_cleared(3,:) = tmp3(boolean1);

% plot
h = histogram(data_cleared(3,:));
h.Normalization = 'Probability';
h.BinWidth = 0.1000;
h.FaceColor = [0.5,0.5,0.5];

% set label
ax = gca;
set(ax,'fontsize',19)
    % usual
% xlabel('Triangle Curvature (\mum^{-1})','FontSize',21);
% ylabel('Frequency','FontSize',21);
    % upper axis
ax.XTickLabel = {'+\infty','2.7','1.8','1.3','0.9','0.6'};
xlabel('Resolution as Log(voxels)','FontSize',21);
ylabel('Frequency','FontSize',21);

% set axis
axis([0,5,0,0.16])
ticks = [0:1:5];

ax.XTick = ticks;
set(ax,'XDir','reverse')
    % disable tick on top and right of the box
        % get handle to current axes
a = gca;
        % set box property to off and remove background color
set(a,'box','off','color','none')
        % create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',ticks,'XTickLabel','','ytick',[]);
        % set original axes as active
axes(a)
        % link axes in case of zooming
linkaxes([a b])



