data = [faceCorresp, FA_diff, FCurv_diff];
data = sortrows(data, [4, 3]);

faceMap_An4 = containers.Map(faceCorresp(:,1), num2cell(faces_An4(faceCorresp(:,1),:),2));
faceMap_An5 = containers.Map(faceCorresp(:,2), num2cell(faces_An5(faceCorresp(:,2),:),2));

%% ##### Plot Triangle Curvature Distribution #####
set(0,'defaultAxesLabelFontSize',1.1)
set(0,'defaultAxesFontSize',19)

FL_An4 = double(h5read(file_An4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
triCurves_An4 =  abs(roundn(h5read(file_An4,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)).';
mask_An4 = ~(any(FL_An4 <= 0, 2) | triCurves_An4 > 100);
triCurves_An4 = triCurves_An4(mask_An4);
FL_An5 = double(h5read(file_An5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
triCurves_An5 =  abs(roundn(h5read(file_An5,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5)).';
mask_An5 = ~(any(FL_An5 <= 0, 2) | triCurves_An5 > 100);
triCurves_An5 = triCurves_An5(mask_An5);

histogram(triCurves_An4,100)
hold on
histogram(triCurves_An5,100)
set(gca, 'YScale', 'log')
aveTriCurv_An4 = sum(triCurves_An4)/length(triCurves_An4);
aveTriCurv_An5 = sum(triCurves_An5)/length(triCurves_An4);
legend(['aveTriCurv, An4 = ', num2str(aveTriCurv_An4)], ['aveTriCurv, An5 = ', num2str(aveTriCurv_An5)])
xlabel('Triangle Curvature, \mum^{-1}')
ylabel('# Triangles')


%% ##### Plot Triangle Curvature Distribution #####
% get the coords of one face pair

idx = 1;
idx_An4 = faceCorresp(idx, 1);
idx_An5 = faceCorresp(idx, 2);
label_An4 = faces_An4(idx_An4, :);
label_An5 = faces_An5(idx_An5, :);

FL_An4 = double(h5read(file_An4,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
triNodes_An4 = 1 + double(h5read(file_An4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
NCoords_An4 = double(h5read(file_An4,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';
FL_An5 = double(h5read(file_An5,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
triNodes_An5 = 1 + double(h5read(file_An5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
NCoords_An5 = double(h5read(file_An5,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';

mask_An4 = ((FL_An4(:,1) == label_An4(1) & FL_An4(:,2) == label_An4(2)) | (FL_An4(:,1) == label_An4(2) & FL_An4(:,2) == label_An4(1)));
mask_An5 = ((FL_An5(:,1) == label_An5(1) & FL_An5(:,2) == label_An5(2)) | (FL_An5(:,1) == label_An5(2) & FL_An5(:,2) == label_An5(1)));

thisFaceNodes_An4 = triNodes_An4(mask_An4,:);
thisFaceNodes_An5 = triNodes_An5(mask_An5,:);
thisFaceNodes_An4 = reshape(thisFaceNodes_An4, [], 1);
thisFaceNodes_An5 = reshape(thisFaceNodes_An5, [], 1);
triNodeCoords_An4 = NCoords_An4(thisFaceNodes_An4,:);
triNodeCoords_An5 = NCoords_An5(thisFaceNodes_An5,:);

tmpNodes_An4 = triNodes_An4(mask_An4,:);
tmpNodes_An5 = triNodes_An5(mask_An5,:);
%%
trisurf( tmpNodes_An4, NCoords_An4(:,1), NCoords_An4(:,2), NCoords_An4(:,3) );
rotate3d on
hold on
trisurf( tmpNodes_An5, NCoords_An5(:,1), NCoords_An5(:,2), NCoords_An5(:,3),'Facecolor',[218/255 168/255 32/255] );
rotate3d on









