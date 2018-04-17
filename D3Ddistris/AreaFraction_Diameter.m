%% for each subset, get a table=[grainDiameter, grainArea]
clear

file = ('/Users/xiaotingzhong/Desktop/091616_STO_1470/Modified/subset1/subset1_misA_FindNagain.dream3d');
dia = h5read(file,'/VoxelDataContainer/FIELD_DATA/EquivalentDiameters');
facelabel_raw = double(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceLabels'));
tri_area_raw = roundn(h5read(file,'/SurfaceMeshDataContainer/FACE_DATA/SurfaceMeshFaceAreas'),-8).';

dia(1,:) = [];
tmp1 = facelabel_raw(1,:);
tmp2 = facelabel_raw(2,:);
boolean = facelabel_raw(1,:) > 0 & facelabel_raw(2,:) > 0;
facelabel(1,:) = tmp1(boolean);
facelabel(2,:) = tmp2(boolean);
tri_area = tri_area_raw(boolean);
tri_area = tri_area.';

GrainArea = zeros(size(dia));

for i = 1:length(tri_area)
    g1 = facelabel(1,i);
    g2 = facelabel(2,i);
    GrainArea(g1) = GrainArea(g1) + tri_area(i);
    GrainArea(g2) = GrainArea(g2) + tri_area(i);
end

table = [dia,GrainArea];

%%
load('sub1_dia_area.mat')
load('sub2_dia_area.mat')
grain_table = [GrainTable_sub1; GrainTable_sub2];
grain_table = sortrows(grain_table, 1);

start = 0;
stepsize = 0.5;
range = ceil(max(GrainTable(:,1)));
plot_table = zeros(range,2);
next = 1;
cnt = 0;
sum_d = 0;
sum_area = 0;

for i = 1:length(grain_table)
    if grain_table(i,1) <= start + stepsize;
        sum_d = sum_d + grain_table(i,1);
        sum_area = sum_area + grain_table(i,2);
        cnt = cnt + 1;
    else
        plot_table(next,1) = sum_d/cnt;
        plot_table(next,2) = sum_area;
        
        start = start + stepsize;
        sum_d = grain_table(i,1);
        sum_area = grain_table(i,2);
        cnt = 1;
        next = next + 1;
    end
end
plot_table(next,1) = sum_d/cnt;
plot_table(next,2) = sum_area;

plot_table(:,2) = plot_table(:,2)/sum(plot_table(:,2));
plot(plot_table(:,1),plot_table(:,2),'o-','linewidth',3,'markersize',10,'markerfacecolor','g')
set(gca,'fontsize',19,'XMinorTick','on');
xlabel('Grain Diameter','FontSize',21);
ylabel('Area Fraction','FontSize',21);
    
    
    
    
    
    
    
    
