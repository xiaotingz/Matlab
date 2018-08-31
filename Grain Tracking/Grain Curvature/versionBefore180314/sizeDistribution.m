% OUTLINE, size distribution plots
    % samples for plot, Austenite, Ferrite, SmallIN100
    % for every sample, two plots containing different set of grain are compared.
        % grainSet1 -- all grains
        % grainSet2 -- grains removed according to centroid
    % in every plot, result of two reconstruction settings are compared.
        % reconsSetting1 -- minSize = 100, minNeigh = 4
        % reconsSetting2 -- minSize = 16,  minNeigh = 2

% CHANEG --- CENTROID! plot title, aveSize line

%% Austenite -- all grains 
D1 = roundn(h5read('/Volumes/RESEARCH/Oct.26 MinNeigh=4/Oct.26 Austenite/Oct26_A_CurvDistri.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
D2 = roundn(h5read('/Volumes/RESEARCH/Nov.4 Austenite/Nov4_A_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
D1(1) = [];
D2(1) = [];

%% Austenite -- data for centroid
D1_all = roundn(h5read('/Volumes/RESEARCH/Oct.26 MinNeigh=4/Oct.26 Austenite/Oct26_A_CurvDistri.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
centroid1 = roundn(h5read('/Volumes/RESEARCH/Oct.26 MinNeigh=4/Oct.26 Austenite/Oct26_A_CurvDistri.dream3d','/VoxelDataContainer/FIELD_DATA/Centroids'),-5).';
D2_all = roundn(h5read('/Volumes/RESEARCH/Nov.4 Austenite/Nov4_A_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
centroid2 = roundn(h5read('/Volumes/RESEARCH/Nov.4 Austenite/Nov4_A_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/Centroids'),-5).';


%% SmallIN100 -- all grains
D1 = roundn(h5read('/Volumes/RESEARCH/Nov.4 SmallIN100/Nov.4 minSize100_minNeigh4/Nov4_nP_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
D2 = roundn(h5read('/Volumes/RESEARCH/Nov.4 SmallIN100/Nov.4 prebuild/Nov4_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
D1(1) = [];
D2(1) = [];

%% SmallIN100 -- data for centroid 
D1_all = roundn(h5read('/Volumes/RESEARCH/Nov.4 SmallIN100/Nov.4 minSize100_minNeigh4/Nov4_nP_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
centroid1 = roundn(h5read('/Volumes/RESEARCH/Nov.4 SmallIN100/Nov.4 minSize100_minNeigh4/Nov4_nP_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/Centroids'),-5).';
D2_all = roundn(h5read('/Volumes/RESEARCH/Nov.4 SmallIN100/Nov.4 prebuild/Nov4_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
centroid2 = roundn(h5read('/Volumes/RESEARCH/Nov.4 SmallIN100/Nov.4 prebuild/Nov4_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/Centroids'),-5).';


%% Ferrite -- all grains
D1 = roundn(h5read('/Volumes/RESEARCH/Oct.26 MinNeigh=4/Oct.26 Ferrite/Processing/Oct26_F_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
D2 = roundn(h5read('/Volumes/RESEARCH/Nov.4 Ferrite/Nov4_F_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
D1(1) = [];
D2(1) = [];

%% Ferrite -- data for centroid
D1_all = roundn(h5read('/Volumes/RESEARCH/Oct.26 MinNeigh=4/Oct.26 Ferrite/Processing/Oct26_F_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
centroid1 = roundn(h5read('/Volumes/RESEARCH/Oct.26 MinNeigh=4/Oct.26 Ferrite/Processing/Oct26_F_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/Centroids'),-5).';
D2_all = roundn(h5read('/Volumes/RESEARCH/Nov.4 Ferrite/Nov4_F_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
centroid2 = roundn(h5read('/Volumes/RESEARCH/Nov.4 Ferrite/Nov4_F_ForCurv.dream3d','/VoxelDataContainer/FIELD_DATA/Centroids'),-5).';


%% clean data according to centroid 
D1_all(1) = [];
D2_all(1) = [];
centroid1(1,:) = [];
centroid2(1,:) = [];
 
aveD1_all = sum(D1_all) / length(D1_all);
aveD2_all = sum(D2_all) / length(D2_all);

% Austenite  --- X=434*0.15; Y=267*0.15; Z=100*0.2;
% Ferrite    --- X=234*0.15; Y=267*0.15; Z=68*0.2;
% SmallIN100 --- X=189*0.25; Y=201*0.25; Z=117*0.25;
X=234*0.15; Y=267*0.15; Z=68*0.2;

delete_bool1 = ones(length(centroid1),1);
for i = 1:length(centroid1)
    if centroid1(i,1) < aveD1_all || X - centroid1(i,1) < aveD1_all
        delete_bool1(i) = 0;
    elseif centroid1(i,2) < aveD1_all || Y - centroid1(i,2) < aveD1_all
        delete_bool1(i) = 0;
    elseif centroid1(i,3) < aveD1_all || Z - centroid1(i,3) < aveD1_all
        delete_bool1(i) = 0;
    end
end
delete_bool1 = logical(delete_bool1);

delete_bool2 = ones(length(centroid2),1);
for i = 1:length(centroid2)
    if centroid2(i,1) < aveD1_all || X - centroid2(i,1) < aveD1_all
        delete_bool2(i) = 0;
    elseif centroid2(i,2) < aveD1_all || Y - centroid2(i,2) < aveD1_all
        delete_bool2(i) = 0;
    elseif centroid2(i,3) < aveD1_all || Z - centroid2(i,3) < aveD1_all
        delete_bool2(i) = 0;
    end
end
delete_bool2 = logical(delete_bool2);

D1 = D1_all(delete_bool1);
D2 = D2_all(delete_bool2);
%% overlap plot --- x axis not normalized
% remember to modify plot title

aveD1 = sum(D1) / length(D1);
aveD2 = sum(D2) / length(D2);

% histogram
h1 = histogram(D1);
h1.Normalization = 'probability';
h1.BinWidth = aveD1/10;
hold on;
h2 = histogram(D2);
h2.Normalization = 'probability';
h2.BinWidth = aveD2/10;
% add two lines for average size
line('XData',[aveD1,aveD1],'YData',[0,0.15],'color',[0,0.4470,0.7410]);
line('XData',[aveD2,aveD2],'YData',[0,0.15],'color',[0.8500,0.3250,0.0980]);
% setting details -- legend, label, title, ticksize
Legend = legend('minSize=100, minNeigh=4','minSize=16, minNeigh=2','aveSize1 -- 100,4','aveSize2 -- 16,2');
set(Legend,'FontSize',12);
xlabel('Grain Diameter, um','FontSize',15);
ylabel('Frequency','FontSize',15);
set(gca,'fontsize',12)
title('Ferrite size distribution -- removed as centroid','FontSize',15);
hold off


%% seperate plots --- x axis normalized
figure(1)
aveD1 = sum(D1) / length(D1);
h1 = histogram(D1);
h1.Normalization = 'probability';
h1.BinWidth = aveD1/10;
% setting details -- legend, label, title, ticksize
xlabel('R/<R>','FontSize',15);
ylabel('Frequency','FontSize',15);
set(gca,'fontsize',12)
% write data tick. write large range, matlab will adjust axis size automatically
ax = gca;
ax.XTick = 0:aveD1/2:20;
ax.XTickLabel = {'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5','5.5','6','6.5','7'};
title('Austenite, removed, 100&4','FontSize',15);

figure(2)
aveD2 = sum(D2) / length(D2);
h2 = histogram(D2);
h2.Normalization = 'probability';
h2.BinWidth = aveD2/10;
set(h2,'facecolor',[0.8500,0.3250,0.0980]);
% setting details -- legend, label, title, ticksize
xlabel('R/<R>','FontSize',15);
ylabel('Frequency','FontSize',15);
set(gca,'fontsize',12)
% write data tick. write large range, matlab will adjust axis size automatically
ax = gca;
ax.XTick = 0:aveD2/2:20;
ax.XTickLabel = {'0','0.5','1','1.5','2','2.5','3','3.5','4','4.5','5','5.5','6','6.5','7'};
title('Austenite, removed, 16&2','FontSize',15);


%% attempt 
% aveSize = sum(grain_diameter_raw)/length(grain_diameter_raw);
% maxSize = max(grain_diameter_raw);
% 
% % decide the #steps and stepSize
% steps = 0;
% stepSize = 0.2;
% while steps*stepSize < maxSize
%     steps = steps + 1;
% end
% 
% % prepare data for plot, first column---start size, second column---frequency
%     % first record frequency as #grains, then normalize it with total #grains to be percentage
% grid = zeros(steps,2);
% frequency = 0;
% start = 0;
% for i = 1:steps
%     for j = 1:length(grain_diameter_raw)
%         if grain_diameter_raw(j) > start && grain_diameter_raw(j) <= (start + stepSize)
%             frequency = frequency + 1;
%         end
%     end
%     frequency = frequency/(length(grain_diameter_raw)-1);
%     
%     grid(i,1) = start;
%     grid(i,2) = frequency;
%     
%     start = start + stepSize;
% end
