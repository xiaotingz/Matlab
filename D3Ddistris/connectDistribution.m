% ############################################################################
% Method
% - From the misorientation of interestget get the objective grain faces. 
%     or equivalently, the two grainIDs [A,B] defining that face.  
% - Identify the grain faces of interest, thats the faces connecting the objective faces
%     this is achieved by searching the common neighbors (a list C) of grainA and
%     grainB. The faces in touch with [A,B] are then defined by [B,C_i] and [A,C_i]
% ############################################################################
clear
O = CrysSym();

file = '/Users/xiaotingzhong/Desktop/Datas/SteelFinal_setTo0/Jan31_Ac_setTo0.dream3d';
thresMisA = 3;
thresRFvec = 0.001;
objMisO = [10, 1, 0, 0];
materialName = 'Austenite ';
currentMisO = [num2str(objMisO(1)), '°@[', sprintf('%d', objMisO(2:4)),'], thres=', num2str(thresRFvec)];
g1 = eye(3);
g2 = AAToG(objMisO(1), objMisO(2:4));
[objRFvec, ~] = dgInFZ(g1, g2, O);


% --------------------------- V6 structure ---------------------------
% numNeigh = h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NumNeighbors');
% NeighborList = h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/NeighborList');
% --------------------------- V4 structure ---------------------------
numNeigh = h5read(file, '/VoxelDataContainer/FIELD_DATA/NumNeighbors');
NeighborList = h5read(file, '/VoxelDataContainer/FIELD_DATA/NeighborList');
EAs = h5read(file, '/VoxelDataContainer/FIELD_DATA/EulerAngles').';
misAs = h5read(file, '/VoxelDataContainer/FIELD_DATA/MisorientationList');
numNeigh(1) = [];
EAs(1,:) = [];
EAs = rad2deg(EAs);


% --------------------------- find the faces of interest by misorientation ---------------------------
% -- first filter all faces by misorientation angle, then check for the axis
faces_tmp = [];
for i = 1:length(numNeigh)
    faces_tmp = vertcat(faces_tmp,ones(numNeigh(i),1)*i);
end
faces = [faces_tmp, NeighborList];

mask = (abs(misAs - objMisO(1)) < thresMisA);
faces_candidate = faces(mask,:);

mask2 = logical(zeros(length(faces_candidate),1));
RFvecs_candidate = zeros(length(faces_candidate),3);
diff_RFs = zeros(length(faces_candidate),3);
for i = 1:length(faces_candidate)
    g1 = EAtoG(EAs(faces_candidate(i,1),:));
    g2 = EAtoG(EAs(faces_candidate(i,2),:));
    [RFvecs_candidate(i,:), ~] = dgInFZ(g1, g2, O);
    diff_RFs(i,:) = RFvecs_candidate(i,:) - objRFvec;
end

% diff_RFs = RFvecs_candidate - repmat(objRFvec,length(RFvecs),1);
diff_norms = sum(diff_RFs.*diff_RFs,2);
mask2 = (diff_norms < thresRFvec);
faces_obj = faces_candidate(mask2, :);


%% --------------------------- identify the grain faces connecting the objective faces ---------------------------
faces_tmp2 = [];
for i = 1:length(faces_obj)
    grainA = faces_obj(i,1);
    grainB = faces_obj(i,2);
    listA = getNeighList(faces_obj(i,1), numNeigh, NeighborList);
    listB = getNeighList(faces_obj(i,2), numNeigh, NeighborList);
    intersection = intersect(listA, listB);
    interestA = [int32(ones(length(intersection),1))*grainA, intersection];
    interestB = [int32(ones(length(intersection),1))*grainB, intersection];
    faces_tmp2 = vertcat(faces_tmp2, interestA, interestB);
end
faces_tmp2 = sort(faces_tmp2,2,'ascend');
faces_interest = unique(faces_tmp2,'rows');

% --------------------------- get misorientation for the faces of interest ---------------------------
RFvecs = zeros(length(faces_interest),3);
misAs = zeros(length(faces_interest),1);  % misA to check calculation
for i = 1:length(faces_interest)
    g1 = EAtoG(EAs(faces_interest(i,1),:));
    g2 = EAtoG(EAs(faces_interest(i,2),:));
    [RFvecs(i,:), misAs(i)] = dgInFZ(g1, g2, O);
end


fig1 = scatter3(RFvecs(:,1),RFvecs(:,2),RFvecs(:,3), ones(length(RFvecs),1));
fig1.MarkerFaceAlpha = .2;
fig1.MarkerEdgeAlpha = .2;
hold on 
scatter3(1/3, 1/3, 1/3, 40, 'r','filled')
scatter3(0.25,0.25,0.0, 40, 'r','filled')
rotate3d on
Name1 = [materialName,', ', currentMisO];
title(Name1,'FontSize',21)
% tmp = repmat([0.25 0.25 0.0],length(RFvecs),1);
% diff_objTmp = (RFvecs - tmp);
% diff_obj = sum(diff_objTmp .* diff_objTmp, 2);
% mask3 = (diff_obj < 0.03);
% sum(mask3)

RFvecNorm = sum(RFvecs .* RFvecs, 2);
figure(2)
fig2 = histogram(misAs,65);
title(Name1,'FontSize',21)
xlabel('misA(objGB, theConnectedGB)', 'Interpreter', 'none','FontSize',21)
ylabel('Count','FontSize',21)
set(gca,'fontsize',19)
set(gca,'FontWeight','bold','linewidth',2)

figure(3)
fig3 = histogram(RFvecNorm,65)
title(Name1,'FontSize',21)
xlabel('norm(RFvecs)', 'Interpreter', 'none','FontSize',21)
ylabel('Count','FontSize',21)
set(gca,'fontsize',19)
set(gca,'FontWeight','bold','linewidth',2)

Info = [currentMisO, ', #objFaces=', num2str(length(faces_obj)), ', #connectionsToObjGB=', num2str(length(RFvecs)), ', #totalFaces=', num2str(length(faces))];
Info