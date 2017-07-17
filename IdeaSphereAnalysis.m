%% distribution data in .dat file: mean value and SD
clear
gbcd_graph = textread('/Users/xiaotingzhong/Desktop/Datas/Mar18 tripleLineSphereCheck/ts9/ts9_curvDistri2_gmt_1.dat');

% delete the first line and lines of MRD=nan
gbcd_graph(1,:) = [];
MRDs_tmp = gbcd_graph(:,3);
bool = ~isnan(MRDs_tmp);
MRDs = MRDs_tmp(bool);

MRD_mean = sum(MRDs)/length(MRDs)

total_MRD = 0;
for i = 1 : length(MRDs)
    total_MRD = total_MRD + (MRD_mean - MRDs(i))^2;
end
MRD_SD = sqrt(total_MRD/length(MRDs))

% num_NaN = length(gbcd_graph) - length(MRDs)
Max = max(gbcd_graph(:,3));
Min = min(gbcd_graph(:,3));

%% Character distribution data in DREAM3D file: mean value and SD
clear

% % version4
CurvDistri = h5read('/Users/xiaotingzhong/Desktop/ts6/ts9_curvDistri2.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
% % version6
% CurvDistri = h5read('/Volumes/RESEARCH/Oct.7 No1-10/No1/No1_CurvDistri.dream3d','/DataContainers/TriangleDataContainer/FaceEnsembleData/GBCD');

j = 1;
for i = 1 : length(CurvDistri)
    if CurvDistri(i,1) ~= 0 || CurvDistri(i,2) ~= 0;
        data(j,1) = i;
        data(j,2) = CurvDistri(i,1);
        data(j,3) = CurvDistri(i,2);
        j = j + 1;
    end
end

mean = sum(data(:,2))/length(data)
total = 0;
for i = 1 : length(data)
    total = total + (mean - data(i,2))^2;
end
SD = sqrt(total/length(data))

% MaxValue = max(data(:,2));
% MinValue = min(data(:,2));



%% triangle curvature data in DREAM3D file
%     mean1 is using #Triangles
%     mean2 is using TriangleAreas
clear

% INPUT: file & sphereID
file = ('/Users/xiaotingzhong/Desktop/Datas/Mar18 tripleLineSphereCheck/ts9/ts9_stats4.dream3d');
sphereID = 9;
facelabel = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
curvature_of_triangle = abs(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshMeanCurvatures'));
triangle_area = roundn(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-5);

data_raw = [facelabel; curvature_of_triangle.'; triangle_area.'];
tmp1 = data_raw(1,:);
tmp2 = data_raw(2,:);
tmp3 = data_raw(3,:);
tmp4 = data_raw(4,:);

boolean = data_raw(1,:) > 0 & data_raw(2,:) > 0 ;
data_cleared(1,:) = tmp1(boolean);
data_cleared(2,:) = tmp2(boolean);
data_cleared(3,:) = tmp3(boolean);
data_cleared(4,:) = tmp4(boolean);

boolean2 = data_cleared(1,:) == sphereID | data_cleared(2,:) == sphereID;
tmp5 = data_cleared(1,:);
tmp6 = data_cleared(2,:);
tmp7 = data_cleared(3,:);
tmp8 = data_cleared(4,:);

data_sphere(1,:) = tmp5(boolean2);
data_sphere(2,:) = tmp6(boolean2);
data_sphere(3,:) = tmp7(boolean2);
data_sphere(4,:) = tmp8(boolean2);

data_sphere = data_sphere.';

total = 0;
for i = 1:length(data_sphere)
    total = total + data_sphere(i,3) * data_sphere(i,4);
end
mean = total/sum(data_sphere(:,4))

total = 0;
for i = 1:length(data_sphere)
    total = total + (mean - data_sphere(i,3))^2;
end
SD = sqrt(total/length(data_sphere))



%% plot for paper --- laplacian comparison
file = 'Dec5.xlsx';
sheet = 1;
ideaCurv_range = 'B5:K5';      %dash line
laplacian1_range = 'B9:K9';    %scatter
laplacian2_range = 'B12:K12';
SD2_range = 'B13:K13';
laplacian3_range = 'B15:K15';
laplacian4_range = 'B18:K18';
ideaCurv = xlsread(file, sheet, ideaCurv_range);
laplacian1 = xlsread(file, sheet, laplacian1_range);
laplacian2 = xlsread(file, sheet, laplacian2_range);
laplacian3 = xlsread(file, sheet, laplacian3_range);
laplacian4 = xlsread(file, sheet, laplacian4_range);
SD2 = xlsread(file, sheet, SD2_range);
xAxis = [0.8,2:10];
% xAxis = [2.3,2.66,3.01,3.33,3.66,3.99,4.32,4.66,4.99,5.32];
% xAxis_range = 'C65:L65'
% xAxis = xlsread(file, sheet, xAxis_range);

figure
plot(xAxis,ideaCurv,'--','Color',[0.5 0.5 0.5]);
hold on
scatter1 = scatter(xAxis,laplacian1,200,[1,0,0],'s','filled');
errorbar(xAxis,laplacian2,SD2,'.','marker','o','markersize',14,'color',[0,0,0],'markerfacecolor',[0,0,0],'markeredgecolor',[0,0,0],'linewidth',1.2);
scatter3 = scatter(xAxis,laplacian3,200,[0,1,0],'d','filled');
scatter4 = scatter(xAxis,laplacian4,200,[0,0,1],'^','filled');

% Legend elsewhere

% Axes setting axis, ect...

axis([0.3,10.5,-0.4,1.4]);
ax = gca;
ax.XTick = xAxis;
% ax.XTickLabel = {'2.4','3.0','3.5','4.0','4.5','5.0','5.5','6.0','6.5','7.0',''};
ax.XTickLabel = {'4','6','9','13','19','28','41','60','88','129',''};

set(gca,'fontsize',19)
xlabel('Sphere Radius (voxels)','FontSize',21);
ylabel('Symmetry Averaged Curvature (\mum^{-1})','FontSize',21);

% % Resolution Range: Austenite maxSize=9.3706; Ferrite maxSize=12.9730
% maxRes = log10((4/3*pi*((1/0.05)^3)/(0.15*0.15*0.2)));
% line('XData',[maxRes,maxRes],'YData',[-0.1,0],'Color',[0.5 0.5 0.5]);


% disable tick on top and right of the box
    % get handle to current axes
a = gca;
    % set box property to off and remove background color
set(a,'box','off','color','none')
    % create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',xAxis,'XTickLabel','','ytick',[]);
    % set original axes as active
axes(a)
    % link axes in case of zooming
linkaxes([a b])


%% Axis for plot --- laplacian comparison
file = 'Dec5.xlsx';
sheet = 1;
ideaCurv_range = 'B5:K5';      %dash line
laplacian1_range = 'B9:K9';    %scatter
laplacian2_range = 'B12:K12';
SD2_range = 'B13:K13';
laplacian3_range = 'B15:K15';
laplacian4_range = 'B18:K18';
ideaCurv = xlsread(file, sheet, ideaCurv_range);
laplacian1 = xlsread(file, sheet, laplacian1_range);
laplacian2 = xlsread(file, sheet, laplacian2_range);
laplacian3 = xlsread(file, sheet, laplacian3_range);
laplacian4 = xlsread(file, sheet, laplacian4_range);
SD2 = xlsread(file, sheet, SD2_range);
xAxis = [0.8,2:10];

figure
plot(xAxis,ideaCurv,'--','Color',[0.5 0.5 0.5]);
hold on
scatter1 = errorbar(xAxis,laplacian2,SD2,'.','marker','s','markersize',13,'color',[1,0,0],'markerfacecolor',[1,0,0],'markeredgecolor',[1,0,0],'linewidth',1.2);
scatter2 = errorbar(xAxis,laplacian2,SD2,'.','marker','o','markersize',13,'color',[0,0,0],'markerfacecolor',[0,0,0],'markeredgecolor',[0,0,0],'linewidth',1.2);
scatter3 = errorbar(xAxis,laplacian2,SD2,'.','marker','d','markersize',13,'color',[0,1,0],'markerfacecolor',[0,1,0],'markeredgecolor',[0,1,0],'linewidth',1.2);
scatter4 = errorbar(xAxis,laplacian2,SD2,'.','marker','^','markersize',13,'color',[0,0,1],'markerfacecolor',[0,0,1],'markeredgecolor',[0,0,1],'linewidth',1.2);


Legend = legend('ideal curvature','smoothing I','smoothing II','smoothing III','smoothing IV');
set(Legend,'box','off','FontSize',16,'position',[0.7,0.62,0.1,0.27]);


axis([0.3,10.5,-0.4,1.4]);
ax = gca;
ax.XTick = xAxis;
ax.XTickLabel = {'S4','S6','S9','S13','S19','S28','S41','S60','S88','S129',''};
set(gca,'fontsize',19)
xlabel('SphereID','FontSize',21);
ylabel('Average Measured Curvature (\mum^{-1})','FontSize',21);

% disable tick on top and right of the box
    % get handle to current axes
a = gca;
    % set box property to off and remove background color
set(a,'box','off','color','none')
    % create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',xAxis,'XTickLabel','','ytick',[]);
    % set original axes as active
axes(a)
    % link axes in case of zooming
linkaxes([a b])
%% plot for paper --- CurvDistri for Smoothing2
% CHANGE xAxis tick label needs if file is modified

file = 'Dec5.xlsx';
sheet = 1;
xAxis_range = 'B48:K48';
ideaCurv_range = 'B5:K5'; 
ideaCurv = xlsread(file, sheet, ideaCurv_range);
laplacian2_range = 'B12:K12';
laplacian2 = xlsread(file, sheet, laplacian2_range);
SD2_range = 'B13:K13';
SD2 = xlsread(file, sheet, SD2_range);
gmt_Rang = 'B43:K44';
gmt = xlsread(file, sheet, gmt_Rang);
xAxis = [0.8,2:10];

plot(xAxis,ideaCurv,'--','Color',[0.5 0.5 0.5]); hold on;
errorbar(xAxis,laplacian2,SD2,'.','marker','o','markersize',14,'color','k','markerfacecolor','k','markeredgecolor','k','linewidth',1.2);
errorbar(xAxis,gmt(1,:),gmt(2,:),'.','marker','o','markersize',14,'color','r','markerfacecolor','r','markeredgecolor','r','linewidth',1.2);
set(gca,'fontsize',19)
xlabel('Log(Number of Voxels per Sphere)','FontSize',21);
ylabel('Symmetry Averaged Curvature (\mum^{-1})','FontSize',21);
axis([0.3,10.5,-0.4,1.4]);
ax = gca;
ax.XTick = xAxis;
ax.XTickLabel = {'2.4','3.0','3.5','4.0','4.5','5.0','5.5','6.0','6.5','7.0'};
% ax.XTickLabel = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10'};

Legend = legend('ideal curvature','triangle curvature',['symmetry averaged' char(10) 'curvature']);
set(Legend,'box','off','FontSize',16,'position',[0.65,0.69,0.1,0.2]);


% disable tick on top and right of the box
    % get handle to current axes
a = gca;
    % set box property to off and remove background color
set(a,'box','off','color','none')
    % create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',xAxis,'XTickLabel','','ytick',[]);
    % set original axes as active
axes(a)
    % link axes in case of zooming
linkaxes([a b])

%% plot for paper --- CurvDistri for large range
% CHANGE xAxis tick label needs if file is modified

file = 'Dec5.xlsx';
sheet = 1;
xAxis_range = 'B48:K48';
% xlRange2 = 'B7:G8';
xlRange3 = 'B50:K51';
xAxis = xlsread(file, sheet, xAxis_range);
% d3d = xlsread(file, sheet, xlRange2);
gmt = xlsread(file, sheet, xlRange3);


% scatter(,data(2,:));
% errorbar(xAxis,d3d(1,:),d3d(2,:),'.','marker','s','markersize',10,'markerfacecolor','auto','linewidth',1.2);
% set(gca,'fontsize',16)
% xlabel('Resolution curvature, um^{-1}','FontSize',17);
% ylabel('Normalized curvature','FontSize',17);
% hold on;

errorbar(xAxis,gmt(1,:),gmt(2,:),'.','marker','o','markersize',10,'markerfacecolor','auto','linewidth',1.2);
set(gca,'fontsize',14)
xlabel('sphereID','FontSize',17);
ylabel('Normalized curvature','FontSize',17);
axis([2.1,7.2,-0.5,2.5]);
ax = gca;
ax.XTick = xAxis;
ax.XTickLabel = {'2.4','3.0','3.5','4.0','4.5','5.0','5.5','6.0','6.5','7.0'};
% ax.XTickLabel = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10'};
line([2,8],[1,1],'LineStyle','--', 'Color',[0.5 0.5 0.5])


% MAX GRAIN SIZE
% Austenite maxSize=9.3706; Ferrite maxSize=12.9730
% maxRes = log10((4/3*pi*(12.9730^3)/(0.15*0.15*0.2)));

% disable tick on top and right of the box
    % get handle to current axes
a = gca;
    % set box property to off and remove background color
set(a,'box','off','color','none')
    % create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',xAxis,'XTickLabel','','ytick',[]);
    % set original axes as active
axes(a)
    % link axes in case of zooming
linkaxes([a b])



%% plot for paper --- curv1.2-0.2 corresponding
% CHANGE xAxis tick label needs if file is modified

file = 'Dec5.xlsx';
sheet = 1;
xAxis_range = 'B57:G57';
mean_range = 'B60:G60';
SD_range = 'B61:G61';
xAxis = xlsread(file, sheet, xAxis_range);
mean = xlsread(file, sheet, mean_range);
SD = xlsread(file, sheet, SD_range);

figure
errorbar(xAxis,mean,SD,'.','marker','o','markersize',10,'markerfacecolor','auto','linewidth',1.2);
line([0,1.4],[1,1],'LineStyle','--', 'Color',[0.5 0.5 0.5])

ax = gca;
set (ax, 'xdir', 'reverse')
% ax.XTickLabel = {'S1','S2','S3','S4','S5','S6','S7','S8','S9','S10'};
ax.XTickLabel = {'','0.2','0.4','0.6','0.8','1.0','1.2',''};
set(gca,'fontsize',14)
xlabel('Resolution Curvature, um^{-1}','FontSize',17);
ylabel('Normalized Curvature','FontSize',17);

% disable tick on top and right of the box
    % get handle to current axes
a = gca;
    % set box property to off and remove background color
set(a,'box','off','color','none')
    % create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
    % set original axes as active
axes(a)
    % link axes in case of zooming
linkaxes([a b])
