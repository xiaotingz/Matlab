% important data:
    % data_grain [grainId, grainDiameter, #Faces, grainCurvature]
    % data_grid [bin_low, bin_high, #Faces, AveCurvature, #Grains, SMD] 
%  notice
    % change step between 40&100 to make different scale plot
%% get data_grid 
start = 0;
width = 1;
step = 40;
data_grid = zeros(step,6);
s_cur = 0;
Numfaces = 0;
cnt = 0;
s_sqdiff = 0;

for i = 1:step
% mean Curv for each Bin
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
% standard mean deviation
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

% what percent of grains are considered in this plot
considered = double(sum(data_grid(:,5))/length(data_grain));

%% Plot 
% plot average curvature of bin
scatter(data_grid(:,3),data_grid(:,4),30,'filled','o','r');
set(gca,'fontsize',13)
xlabel('Grain Faces (F)','FontSize',15);
ylabel('Integral Mean Curvature (M), \mum','FontSize',15);
% plot standard Standard Deviation
range = zeros(length(data_grid),2);
for i = 1:length(data_grid)
    range(i,1) = data_grid(i,4) + data_grid(i,6); % up value
    range(i,2) = data_grid(i,4) - data_grid(i,6); % down value
    line([data_grid(i,3),data_grid(i,3)],[range(i,1),range(i,2)],'color','r');
    hold on
end


line([0,40],[0,0],'LineStyle','--', 'Color',[0.5 0.5 0.5])

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
