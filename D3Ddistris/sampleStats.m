% Contents

% histogram of Dihedral angles
% histogram of DAD: Dihedral angles in Character distribution
% Combine two Ferrite volume, Dihedral angles
% ECDF
% Normalized Size Distri, removed & not remove
% Log Size Distri
% #Faces distribution
% Topology vs Size
%% histogram of Dihedral angles
clear
data = textread('060216_DAlist_Ni.txt');
data(1,:) = [];
DihedralAngle = data(:,2);
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

%% histogram of DAD: Dihedral angles in Character distribution
clear
DAD = h5read('/Users/xiaotingzhong/Desktop/for_FindDAD/060216_Ni_DAdistri.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
DAD = DAD(:,1);
h = histogram(DAD);
h.Normalization = 'probability';
h.BinWidth = 1;

xlim([0,180])
ax = gca;
ax.XTick = [60:30:180];
xlabel('Dihedral Angle(Bined), °','FontSize',20);
ylabel('Frequency','FontSize',20);
set(gca,'fontsize',18)
title('Dihedral Angle Distritbution(Bined), Nickel','FontSize',18);

%% histogram of DAD: Dihedral angles in Character distribution? weighted by #observations
DAD = h5read('/Users/xiaotingzhong/Desktop/for_FindDAD/Jan31_Frun4_DAdistri.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
DAD = DAD(:,1);
Counter = h5read('/Users/xiaotingzhong/Desktop/for_FindDAD/Jan31_Frun4_DAdistri.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCDCounters');
Counter = Counter(:,1);

weights = zeros(180,1);
for i = 1:length(DAD)
%     use upper limit
    tmp = floor(DAD(i)) + 1;
    weights(tmp) = weights(tmp) + Counter(i);
end
bar([1:1:180],weights);

ax = gca;
ax.XTick = [60:30:180];
xlabel('Binned Dihedral Angle, °','FontSize',20);
ylabel('Frequency','FontSize',20);
set(gca,'fontsize',18)
title('Binned Dihedral Angle Weighted by Counter, Ferrite run4','FontSize',18);







%% Combine two Ferrite volume, Dihedral angles
data_F = textread('DAlist_F.txt');
data_F(1,:) = [];
DA_F = data_F(:,2);
data_run4 = textread('DAlist_Frun4.txt');
data_run4(1,:) = [];
DA_run4 = data_run4(:,2);

DA_F = DA_F(:,1);
DA_run4 = DA_run4(:,1);

DihedralAngle = [DA_F;DA_run4];
%% Combine two Ferrite volume, Dihedral angles as in character distribution
DAD_F = h5read('/Users/xiaotingzhong/Desktop/for_FindDAD/Jan31_F_DAdistri.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
Counter_F = h5read('/Users/xiaotingzhong/Desktop/for_FindDAD/Jan31_F_DAdistri.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCDCounters');
DAD_run4 = h5read('/Users/xiaotingzhong/Desktop/for_FindDAD/Jan31_Frun4_DAdistri.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCD');
Counter_run4 = h5read('/Users/xiaotingzhong/Desktop/for_FindDAD/Jan31_Frun4_DAdistri.dream3d','/SurfaceMeshDataContainer/ENSEMBLE_DATA/GBCDCounters');
DAD_F = DAD_F(:,1);
Counter_F = Counter_F(:,1);
DAD_run4 = DAD_run4(:,1);
Counter_run4 = Counter_run4(:,1);

totalDAD = zeros(length(DAD_F),1);
for i = 1 : length(DAD_F)
    totalDAD(i) = (DAD_F(i)*Counter_F(i) + DAD_run4(i)*Counter_run4(i)) / (Counter_F(i)+Counter_run4(i));
end



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