[~, ~, faceCorresp_AllF] = TrackFace(file_An4, file_An5, lookUp, false);
[~, ~, faceCorresp_CF] = TrackFace(file_An4, file_An5, lookUp, true);
% 
% FAcomplete_diff = FItgCurvs_An5(faceCorresp_full(:,2),1) - FItgCurvs_An4(faceCorresp_full(:,1),1);
% FAall_diff = FItgCurvs_An5(faceCorresp_full(:,2),1) - FItgCurvs_An4(faceCorresp_full(:,1),1);
% sum(FAcomplete_diff)
% 
% CInfo = [C_An4, faces_An4(faceCorresp_full(:,1),:), zeros(length(C_An4),1), C_An5, faces_An5(faceCorresp_full(:,2),:)];
% CInfo = sortrows(CInfo, 1);

faceCorresp = faceCorresp_CF;
% FItgCurvs_An4 = faceCurvatureForTrack(file_An4, faces_An4);
% FItgCurvs_An5 = faceCurvatureForTrack(file_An5, faces_An5);
diff_tmp = FItgCurvs_An5(faceCorresp(:,2),:) - FItgCurvs_An4(faceCorresp(:,1),:);
FA_diff = diff_tmp(:,1);


plotBinData(FItgCurvs_An4(faceCorresp(:,1),2), FA_diff, [0, 200], ylim_bin, 5, 'Integral |Face Curvature| in An4, \mum', 'Face Area Difference, \mum^{2}', true, false)
print('FItgCurv_FADiff_completeFaces', '-dpng','-r300')



