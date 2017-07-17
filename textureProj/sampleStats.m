%% histogram
h = histogram(DihedralAngle);
h.Normalization = 'probability';
h.BinWidth = 2;

xlim([0,180])
ax = gca;
ax.XTick = [0:30:180];
xlabel('Dihedral Angle, °','FontSize',20);
ylabel('Frequency','FontSize',20);
set(gca,'fontsize',18)
title('Dihedral Angle Distritbution, Nickel','FontSize',18);




%% ECDF
[fA,xA] = ecdf(DA_A);
plot(xA,fA);
hold on

[fF,xF] = ecdf(DA_F);
plot(xF,fF);

xlim([0,180])
ax = gca;
ax.XTick = [0:30:180];
legend('Austenite','Ferrite','Location','northwest');
xlabel('Dihedral Angle, °','FontSize',20);
ylabel('Cumulative Probability','FontSize',20);
set(gca,'fontsize',18)
title('ECDF','FontSize',18);
hold off



%%
% Austenite = '/Volumes/DATAS/Datas/Jan.31 Austenite/Jan31_A_stats.dream3d'
% Ferrite = '/Volumes/DATAS/Datas/Ferrite/Jan.31 Ferrite/Jan31_Fc_stats.dream3d'
clear
file = ('/Users/xiaotingzhong/Desktop/DAdistri data/Jan31_Fc_stats.dream3d');
centroids = roundn(h5read(file,'/VoxelDataContainer/FIELD_DATA/Centroids'),-5).';
grain_diameter_raw = roundn(h5read(file,'/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
num_of_neigh = double(h5read(file,'/VoxelDataContainer/FIELD_DATA/NumNeighbors'));

%% Normalized Size Distri, removed & not remove
aveD_raw = sum(grain_diameter_raw) / length(grain_diameter_raw);
% Austenite --- X=434*0.15; Y=267*0.15; Z=100*0.2;
% Ferrite   --- X=234*0.15; Y=267*0.15; Z=68*0.2;
X=434*0.15; Y=267*0.15; Z=100*0.2;

delete_bool = ones(length(centroids),1);
for i = 1:length(centroids)
    if centroids(i,1) < aveD_raw || X - centroids(i,1) < aveD_raw
        delete_bool(i) = 0;
    elseif centroids(i,2) < aveD_raw || Y - centroids(i,2) < aveD_raw
        delete_bool(i) = 0;
    elseif centroids(i,3) < aveD_raw || Z - centroids(i,3) < aveD_raw
        delete_bool(i) = 0;
    end
end
delete_bool = logical(delete_bool);


x_raw = (grain_diameter_raw/aveD_raw);
h1 = histogram(x_raw(2:length(x_raw)));
h1.Normalization = 'probability';
h1.BinWidth = aveD_raw/20;


hold on
xA = grain_diameter_raw(delete_bool);
aveD = sum(xA)/length(xA);
xA = (xA/aveD);
h2 = histogram(xA);
h2.Normalization = 'probability';
h2.BinWidth = aveD/20;
Legend = legend('not removed','removed');
set(Legend,'FontSize',12);

xlabel('D/<D>','FontSize',20);
ylabel('Frequency','FontSize',20);
set(gca,'fontsize',18)
title('Normalized Grain Size Distribution, Ferrite','FontSize',18);
hold off


%% Log Size Distri
xA = grain_diameter_raw;
xA(1) = [];
aveD = sum(xA)/length(xA);
xA = log(xA/aveD);
h2 = histogram(xA);
h2.Normalization = 'probability';

xlabel('Log( D/<D> )','FontSize',20);
ylabel('Frequency','FontSize',20);
set(gca,'fontsize',18)
title('Log Size Disribution, Ferrite','FontSize',18);
hold off

%% #Faces distribution
xA = num_of_neigh;
xA(1) = [];
h2 = histogram(xA);
h2.Normalization = 'probability';
h2.BinWidth = 2;

xlabel('Number of Faces','FontSize',20);
ylabel('Frequency','FontSize',20);
set(gca,'fontsize',18)
title('Topology Distribution, Austenite','FontSize',18);
hold off

%% Topology vs Size
grain_diameter_raw(1) = [];
num_of_neigh(1) = [];
scatter(grain_diameter_raw,num_of_neigh,30,'filled');

xlabel('Grain Diameter, \mum^{-1}','FontSize',20);
ylabel('Number of Faces','FontSize',20);
set(gca,'fontsize',18)
title('Ferrite','FontSize',18);
hold off
xlim([0.8,14])


%% write txt file
fileID = fopen('DA_A.txt','w');
fprintf(fileID,'%12.8f\n',DA_A);
fclose(fileID);

fileID = fopen('DA_F.txt','w');
fprintf(fileID,'%12.8f\n',DA_F);
fclose(fileID);