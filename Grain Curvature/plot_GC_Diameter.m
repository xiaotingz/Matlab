%  important data:
    % data_grain [grainId, grainDiameter, #Faces, grainCurvature]
    % data_grid [bin_low, bin_high, AveDiameter, AveCurvature, bin_#Grains, SMD] 

%% get data_grid 
% Plot GrainCurvature vs GrainDiameter 
start = 1.0; % notice, Austenite & Ferrite Sample start with 0.9 
width = 0.2;
%  notice
    % step=100 for fullsize ; Austenite&Ferrite=35, Mg=44,
step = 25;
data_grid = zeros(step,6);
s_cur = 0;
s_diameter = 0;
cnt = 0;
s_sqdiff = 0;

for i = 1:step
% mean Curv for each Bin
    for j = 1:length(data_grain)
        if data_grain(j,2) >= start && data_grain(j,2) < start+width
            s_cur = s_cur + data_grain(j,4);
            s_diameter = s_diameter + data_grain(j,2);
            cnt = cnt + 1;
        end
    end
    data_grid(i,1) = start;
    data_grid(i,2) = start + width;
    data_grid(i,3) = s_diameter/cnt;
    data_grid(i,4) = s_cur/cnt;
    data_grid(i,5) = cnt;
% standard mean deviation
    for j = 1:length(data_grain)
        if data_grain(j,2) >= start && data_grain(j,2) < start+width
            s_sqdiff = s_sqdiff + (data_grain(j,4) - data_grid(i,4))^2;
        end
    end
    data_grid(i,6) = sqrt(s_sqdiff/cnt);
    
    start = start + width;
    cnt = 0;
    s_cur = 0;
    s_diameter = 0;
    s_sqdiff = 0;
end

% what percent of grains are considered in this plot
considered = sum(data_grid(:,5))/length(data_grain);

%% Plot
% plot point: average curvature of bin
scatter(data_grid(:,3),data_grid(:,4),30,'filled');
set(gca,'fontsize',13)
xlabel('Diameter of Grains, \mu','FontSize',15);
ylabel('Integral Mean Curvature, \mu^{-1}','FontSize',15);

% axis([1,14,-0.4,1.2]);
% plot line: standard mean deviation
range = zeros(length(data_grid),2);
for i = 1:length(data_grid)
    range(i,1) = data_grid(i,4) + data_grid(i,6); % up value
    range(i,2) = data_grid(i,4) - data_grid(i,6); % down value
    line([data_grid(i,3),data_grid(i,3)],[range(i,1),range(i,2)]);
    hold on
end

line([1,6],[0,0],'LineStyle','--', 'Color',[0.5 0.5 0.5])

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
