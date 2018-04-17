% ###############################################################
% data from the correlation program, read by import data
% ###############################################################
% There are 36 copies for each bin in the reduced zone. 
% Set default thres = 3*36
obj = GBCD_GBCurvD;
objName = 'GBCD_GBCurvD';
thres = 36*3;

x = obj(:,1);
y = obj(:,2);
mask = (obj(:,3) > thres);
x = x(mask);
y = y(mask);

[xlabelName, ylabelName] = setCorrAxis(objName);
scatter(x, y, 120, 's', 'filled', 'k')
xlabel(xlabelName, 'FontSize', 21)
ylabel(ylabelName, 'FontSize', 21)
% set(gca,'fontsize',19 ,'FontWeight','bold','linewidth',2)
set(gca,'fontsize',19)
box on 





%% All data, SD
% file = textread('Jan31_A_Curv10_v_energy_TWIP.txt');
file = textread('Jan31_A_Curv10_v_Jan31_A_gbcd10.txt');

cnt = 1;
figure(2)
for i = 1:length(file)
    if file(i,3) >= 1000
        scatter(file(i,1),file(i,2),100,'filled','k');hold on
        line([file(i,1),file(i,1)],[(file(i,2)-file(i,5)),(file(i,2)+file(i,5))],'color','k');
        line([file(i,1)-0.003,file(i,1)+0.003],[(file(i,2)-file(i,5)),(file(i,2)-file(i,5))],'color','k')
        line([file(i,1)-0.003,file(i,1)+0.003],[(file(i,2)+file(i,5)),(file(i,2)+file(i,5))],'color','k')
        F_energyPlot(cnt,1) = file(i,2);
        cnt = cnt + 1;
    end
end
set(gca,'fontsize',18)

% AXIS LABEL
xlabel('Nearest Neighbor Topology Difference (F - mF)','FontSize',15);
ylabel('Average Energy, a.u.','FontSize',15);
% xlabel('Population, MRD','FontSize',17);



% ax = gca;
% Ferrite, no Errorbar
% ax.XTick = [0.5:0.1:0.8];
% axis([0.47,0.8,0.75,1.05]); 
% Austenite, no Errorbar
% ax.YTick = [0.59:0.02:0.69];
% axis([0.45,0.67,0.59,0.69]);
% box on

% Austenite, Errorbar
% axis([0.45,0.67,0.46,0.82]);
% Ferrite, Errorbar
axis([0.47,0.79,0.57,1.15]);


% disable tick on top and right of the box
a = gca;
    % set box property to off and remove background color
set(a,'box','off','color','none')
    % create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
    % set original axes as active
axes(a)
    % link axes in case of zooming
linkaxes([a b])


%% All data, noSD
% file = textread('Jan31_A_Curv10_v_energy_TWIP.txt');
file = textread('Jan31_A_Curv10_v_Jan31_A_gbcd10.txt');

cnt = 1;
figure(2)
for i = 1:length(file)
    if file(i,3) >= 1000
        scatter(file(i,1),log(file(i,2)),100,'filled','k');hold on
    end
end
set(gca,'fontsize',18)

% AXIS LABEL
xlabel('Mean Curvature Range, \mum^{-1}','FontSize',18);
ylabel('Ln(Average Population, MRD)','FontSize',18);



% ax = gca;
% Ferrite, no Errorbar
% ax.XTick = [0.5:0.1:0.8];
% axis([0.47,0.8,0.75,1.05]); 
% Austenite, no Errorbar
% ax.YTick = [0.59:0.02:0.69];
% axis([0.45,0.67,0.59,0.69]);
% box on

% Austenite, Errorbar
% axis([0.45,0.67,0.46,0.82]);
% Ferrite, Errorbar
% axis([0.47,0.79,0.57,1.15]);


% disable tick on top and right of the box
a = gca;
    % set box property to off and remove background color
set(a,'box','off','color','none')
    % create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
    % set original axes as active
axes(a)
    % link axes in case of zooming
linkaxes([a b])

%% Zone data
clc
clear
% read in data
curvature = textread('/Users/xiaotingzhong/Desktop/programs/ZoneComparison/extracter/Jan31_A_Curv10_111.txt');
curvature = curvature(:,2);
curvature = curvature(1:31);
energy = textread('/Users/xiaotingzhong/Desktop/programs/ZoneComparison/extracter/energy_TWIP_111.txt');
energy = energy(:,2);
energy = energy(1:31);
angles = [0:2:60].';

% plot with secondary axis
[Ax,Line1,Line2] = plotyy(angles,curvature,angles,energy,'plot');
Line1.Color = 'k';
Line1.Marker = 'o';
Line1.MarkerEdgeColor = 'k';
Line1.MarkerSize  = 10;
Line2.Color = 'k';
Line2.Marker = 'o';
Line2.MarkerEdgeColor = 'k';
Line2.MarkerFaceColor  = 'k';
Line2.MarkerSize  = 10;

% set axes details: font, color, tick...
set(Ax,{'ycolor'},{'k';'k'})
set(Ax(1),'fontsize',19);
ylim(Ax(1),[1.15,1.23])
Ax(1).YTick = (1.15:0.02:1.23);
Ax(1).XTick = (0:15:90);
box(Ax(1),'off')
set(Ax(2),'fontsize',19)
ylim(Ax(2),[0.4,0.7])
Ax(2).YTick = (0.4:0.1:0.7);
Ax(2).XTick = (0:15:90);
upperLine = ones(length(angles),1);
upperLine(:) = 0.7;
h3 = line(angles,upperLine,'Parent',Ax(2),'Color','k');


% set legend and axes labels
legend('Curvature','Energy');
legend boxoff
xlabel('Angle from (110) along [111] Zone (°)','FontSize',21);
ylabel(Ax(1),'Grain Bounary Curvature (\mum^{-1})','FontSize',21);
ylabel(Ax(2),'Grain Boundary Energy (a.u.)','FontSize',21);

%% for Legend
figure
X = [0,19.1,30,40.9,60];
Y = [1,1,1,1,1];
scatter(X,Y)

axis([0,60,0,2]);
ax = gca;
ax.XTick = X;
ax.XTickLabel = {'(110)','(321)','(211)','(312)','(101)'};
set(gca,'fontsize',19)
xlabel('Grain Boundary Planes','FontSize',21);

box on

