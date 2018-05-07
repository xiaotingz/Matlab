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