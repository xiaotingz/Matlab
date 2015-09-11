%  important data:
    % data_grain [grainId, grainDiameter, #Faces, grainCurvature]
    % data_grid [bin_low, bin_high, AveDiameter, AveCurvature, bin_#Grains, SMD] 
%  notice
    % change step between 47&100 to make different scale plot

% Plot GrainCurvature vs GrainDiameter 
start = 0.9; % notice, Ferrite Sample need to start with 0.9 to count out a few extreme grain
width = 0.2;
step = 100;
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


% plot point: average curvature of bin
figure('name','Grain Curvature vs Diameter')
scatter(data_grid(:,3),data_grid(:,4),30,'filled');
xlabel('Grain Diameter','FontSize',13);
ylabel('Grain Curvature','FontSize',13);
grid on
% plot line: standard mean deviation
range = zeros(step,2);
for i = 1:step
    range(i,1) = data_grid(i,4) + data_grid(i,6); % up value
    range(i,2) = data_grid(i,4) - data_grid(i,6); % down value
    line([data_grid(i,3),data_grid(i,3)],[range(i,1),range(i,2)]);
    hold on
end

