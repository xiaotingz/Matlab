% important data:
    % data_grain [grainId, grainDiameter, grainCurvature]
    % data_grid [bin_low, bin_high, #Faces, AveCurvature, #Grains, SMD] 
%  notice
    % change step between 40&100 to make different scale plot

start = 0;
width = 1;
step = 100;
data_grid = zeros(step,6);
s_cur = 0;
s_faces = 0;
cnt = 0;
s_sqdiff = 0;

for i = 1:step
% mean Curv for each Bin
    for j = 1:length(data_grain)
        if data_grain(j,3) >= start && data_grain(j,3) < start+width
            s_cur = s_cur + data_grain(j,4);
            s_faces = s_faces + data_grain(j,3);
            cnt = cnt + 1;
        end
    end
    data_grid(i,1) = start;
    data_grid(i,2) = start + width;
    data_grid(i,3) = s_faces/cnt;
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
    s_faces = 0;
    s_sqdiff = 0;
end

% what percent of grains are considered in this plot
considered = double(sum(data_grid(:,5))/length(data_grain));


% plot average curvature of bin
figure('name','Grain Curvature vs #Faces')
scatter(data_grid(:,3),data_grid(:,4),30,'filled');
xlabel('#Faces','FontSize',13);
ylabel('Grain Curvature','FontSize',13);
grid on
% plot standard mean deviation
range = zeros(step,2);
for i = 1:step
    range(i,1) = data_grid(i,4) + data_grid(i,6); % up value
    range(i,2) = data_grid(i,4) - data_grid(i,6); % down value
    line([data_grid(i,3),data_grid(i,3)],[range(i,1),range(i,2)]);
    hold on
end
axis([0,70,-0.5,3]);
% % % Plot Grain Diameter vs Number of Faces
% % figure('name','Grain Diameter vs Number of Faces')
% % scatter(grain_diameter_raw,num_of_neigh,30,'filled');
% % xlabel('Grain Diameter','FontSize',13);
% % ylabel('Number of Neighbor','FontSize',13);
% % axis([0.9,10,0,100]);
% % grid on
