% important data:
    % data_grain [grainId, grainDiameter, #Faces, grainCurvature]
    % data_grid [bin_low, bin_high, #Faces, AveCurvature, #Grains, SMD] 
%  notice
    % change step between 40&100 to make different scale plot
%% get data_grid 
start = 0;
width = 3;
step = max(data_grain(:,3))+1;
data_grid = zeros(step,6);
s_cur = 0;
Numfaces = 0;
cnt = 0;
s_sqdiff = 0;

for i = 1:step
% -- mean Curv for each Bin
    for j = 1:length(data_grain)
        if data_grain(j,3) >= start && data_grain(j,3) < start+width
            s_cur = s_cur + data_grain(j,4);
            Numfaces = Numfaces + data_grain(j,3);
            cnt = cnt + 1;
        end
    end
    data_grid(i,1) = start;
    data_grid(i,2) = start + width;
    data_grid(i,3) = Numfaces/cnt;
    data_grid(i,4) = s_cur/cnt;
    data_grid(i,5) = cnt;
% -- standard mean deviation
    for j = 1:length(data_grain)
        if data_grain(j,3) >= start && data_grain(j,3) < start+width
            s_sqdiff = s_sqdiff + (data_grain(j,4) - data_grid(i,4))^2;
        end
    end
    data_grid(i,6) = sqrt(s_sqdiff/cnt);
    
    start = start + width;
    cnt = 0;
    s_cur = 0;
    Numfaces = 0;
    s_sqdiff = 0;
end

% -- what percent of grains are considered in this plot
considered = double(sum(data_grid(:,5))/length(data_grain));

%% Plot 
% -- plot average curvature of bin
scatter(data_grid(:,3),data_grid(:,4),50,'filled','o','k');
ax = gca;
% ax.XTick = [0:5:40];

% -- plot standard Standard Deviation
range = zeros(length(data_grid),2);
for i = 1:length(data_grid)
    range(i,1) = data_grid(i,4) + data_grid(i,6); % up value
    range(i,2) = data_grid(i,4) - data_grid(i,6); % down value
    line([data_grid(i,3),data_grid(i,3)],[range(i,1),range(i,2)],'color','k');
    hold on
end

% -- axis labels
xlabel('F','FontSize',21);
ylabel('M_{S}, micrometers','FontSize',21);
% -- axis options
set(ax,'fontsize',19)
line([0,max(data_grain(:,3))+0.2],[0,0],'LineStyle','--', 'Color',[0.5 0.5 0.5],'linewidth',2)
set(gca,'FontWeight','bold','linewidth',2)
% xlim([0,40])
% ylim([-25,5])
% text(1.5,-22,'(b)','FontWeight','bold','FontSize',30)
% disable tick on top and right of the box
    % get handle to current axes
a = gca;
    % set box property to off and remove background color
set(a,'box','off','color','none')
    % create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[],'FontWeight','bold','linewidth',2);
    % set original axes as active
axes(a)
    % link axes in case of zooming
linkaxes([a b])
%% for Legend
% one = [1:1:15];
% two = [2:1:16];
% SD = zeros(15,1).';
% scatter1 = errorbar(one,one,SD,'.','marker','s','markersize',6,'color','k','markerfacecolor','k','markeredgecolor','k','linewidth',1.2);
% hold on
% scatter2 = errorbar(one,two,SD,'.','marker','o','markersize',6,'color','r','markerfacecolor','r','markeredgecolor','r','linewidth',1.2);
% 
% 
% Legend = legend('Ferrite','Austenite');
% set(Legend,'box','off','FontSize',16,'position',[0.7,0.4,0.1,0.15]);
