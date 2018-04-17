D_Araw = roundn(h5read('/Volumes/RESEARCH/Grain Curvature/MinSize=100/Jul.20 Austenite/processing/Jul20_ForCurv_centroid.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
centr_A = roundn(h5read('/Volumes/RESEARCH/Grain Curvature/MinSize=100/Jul.20 Austenite/processing/Jul20_ForCurv_centroid.dream3d','/VoxelDataContainer/FIELD_DATA/Centroids'),-5).';
nneigh_Araw = double(h5read('/Volumes/RESEARCH/Grain Curvature/MinSize=100/Jul.20 Austenite/processing/Jul20_ForCurv_centroid.dream3d','/VoxelDataContainer/FIELD_DATA/NumNeighbors'));
D_Fraw = roundn(h5read('/Volumes/RESEARCH/Grain Curvature/MinSize=100/Jun.24 Ferrite/Jun.24 ForCurv_centroid.dream3d','/VoxelDataContainer/FIELD_DATA/EquivalentDiameters'),-5);
centr_F = roundn(h5read('/Volumes/RESEARCH/Grain Curvature/MinSize=100/Jun.24 Ferrite/Jun.24 ForCurv_centroid.dream3d','/VoxelDataContainer/FIELD_DATA/Centroids'),-5).';
nneigh_Fraw = double(h5read('/Volumes/RESEARCH/Grain Curvature/MinSize=100/Jun.24 Ferrite/Jun.24 ForCurv_centroid.dream3d','/VoxelDataContainer/FIELD_DATA/NumNeighbors'));


% clean data according to centroid 
D_Araw(1) = [];
D_Fraw(1) = [];
centr_A(1,:) = [];
centr_F(1,:) = [];
nneigh_Araw(1,:) = [];
nneigh_Fraw(1,:) = [];
 
aveD_Aall = sum(D_Araw) / length(D_Araw);
aveD_Fall = sum(D_Fraw) / length(D_Fraw);

% Austenite  --- X=434*0.15; Y=267*0.15; Z=100*0.2;
% Ferrite    --- X=234*0.15; Y=267*0.15; Z=68*0.2;
% SmallIN100 --- X=189*0.25; Y=201*0.25; Z=117*0.25;
X_A=434*0.15; Y_A=267*0.15; Z_A=100*0.2;
X_F=234*0.15; Y_F=267*0.15; Z_F=68*0.2;

delete_bool1 = ones(length(centr_A),1);
for i = 1:length(centr_A)
    if centr_A(i,1) < aveD_Aall || X_A - centr_A(i,1) < aveD_Aall
        delete_bool1(i) = 0;
    elseif centr_A(i,2) < aveD_Aall || Y_A - centr_A(i,2) < aveD_Aall
        delete_bool1(i) = 0;
    elseif centr_A(i,3) < aveD_Aall || Z_A - centr_A(i,3) < aveD_Aall
        delete_bool1(i) = 0;
    end
end
delete_bool1 = logical(delete_bool1);

delete_bool2 = ones(length(centr_F),1);
for i = 1:length(centr_F)
    if centr_F(i,1) < aveD_Fall || X_F - centr_F(i,1) < aveD_Fall
        delete_bool2(i) = 0;
    elseif centr_F(i,2) < aveD_Fall || Y_F - centr_F(i,2) < aveD_Fall
        delete_bool2(i) = 0;
    elseif centr_F(i,3) < aveD_Fall || Z_F - centr_F(i,3) < aveD_Fall
        delete_bool2(i) = 0;
    end
end
delete_bool2 = logical(delete_bool2);

D_A = D_Araw(delete_bool1);
D_F = D_Fraw(delete_bool2);
nneigh_A = nneigh_Araw(delete_bool1);
nneigh_F = nneigh_Fraw(delete_bool2);

AveD_A = sum(D_A)/length(D_A);
AveD_F = sum(D_F)/length(D_F);
AveNN_A = sum(nneigh_A)/length(nneigh_A);
AveNN_F = sum(nneigh_F)/length(nneigh_F);

% bin data for plot
% faces = [#faces, A_freq, F_freq]
faces = zeros(40,3);
for i = 1:40
    faces(i,1) = i;
    cntA = 0;
    for j_A = 1:length(nneigh_A)
        if nneigh_A(j_A) == i
            cntA = cntA + 1;
        end
    end
    faces(i,2) = cntA/length(nneigh_A);
    cntF = 0;
    for j_F = 1:length(nneigh_F)
        if nneigh_F(j_F) == i
            cntF = cntF + 1;
        end
    end
    faces(i,3) = cntF/length(nneigh_F);    
end

% D = [norm_D, A_freq, F_freq];  normD = BinSize*D/aveD, BinSize = 1/13 to cover 3*D range
Ds = zeros(40,3);
bin_A = AveD_A/13;
bin_F =  AveD_F/13;
for i = 1:40
    Ds(i,1) = i;
    cntA = 0;
    for j_A = 1:length(D_A)
        if D_A(j_A) >= (i-1)*bin_A && D_A(j_A) < i*bin_A
            cntA = cntA + 1;
        end
    end
    Ds(i,2) = cntA/length(nneigh_A);
    cntF = 0;
    for j_F = 1:length(D_F)
        if D_F(j_F) >= (i-1)*bin_F && D_F(j_F) < i*bin_F
            cntF = cntF + 1;
        end
    end
    Ds(i,3) = cntF/length(nneigh_F);
end

%% plot bin data
plot(faces(:,1),faces(:,2),'-o','linewidth',2,'color','k','markersize',5,'markerfacecolor','k')
hold on
plot(faces(:,1),faces(:,3),'-d','linewidth',2,'color','r','markersize',5,'markerfacecolor','r')
plot(Ds(:,1),Ds(:,2),'-s','linewidth',2,'color','b','markersize',5,'markerfacecolor','b')
plot(Ds(:,1),Ds(:,3),'-^','linewidth',2,'color','g','markersize',5,'markerfacecolor','g')
Legend = legend('A, #faces','F, #faces','A, normalized D(1/13*40)','F, normalized D(1/13*40)');
set(Legend,'FontSize',12);

%% histogram
h1 = histogram(nneigh_A);
h1.Normalization = 'probability';
h1.BinWidth = 1;
hold on;
axis([0,40,0,0.1]);
h1 = histogram(nneigh_F);
h1.Normalization = 'probability';
h1.BinWidth = 1;
Legend = legend('Austenite','Ferrite');
set(Legend,'FontSize',12);


