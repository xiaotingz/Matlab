set(0,'defaultAxesFontSize',16)
% load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/data_matlab/190425_Hsmooth_geo_topo_an5crop2.mat', ...
%     'tracked_uniqueface_an4', 'face_area_diff', 'face_itg_abscurv_diff');

% --------------------------- full data ---------------------------
load('/Volumes/XIAOTING/Ni/working/190621_tracked_faces_full.mat')
dists = dists_an4_full;
% --------------------------- data ---------------------------
% load('/Volumes/XIAOTING/Ni/working/190425_Hsmooth_geo_topo_an5crop2.mat', ...
%     'tracked_uniqueface_an4', 'tracked_uniqueface_an5','face_area_diff', ...
%     'face_area_an4', 'face_itg_abscurv_diff');
% load('/Volumes/XIAOTING/Ni/working/190425_dist_from_misorientations.mat')
% face_area_an5 = face_area_an4 + face_area_diff;
% --------------------------- clean data ---------------------------
% load('/Volumes/XIAOTING/Ni/working/An4new6_fixOrigin3_GBCDtris_noSmallFace.mat', ...
%     'face_area_diff', 'face_itg_abscurv_diff', 'tracked_uniqueface_an4')
% load('/Volumes/XIAOTING/Ni/working/190425_Hsmooth_geo_topo_an5crop2.mat', 'face_area_an4');
% tracked_uniqueface_an4_clean = tracked_uniqueface_an4;
% load('/Volumes/XIAOTING/Ni/working/190425_dist_from_misorientations.mat')
% mask = ismember(tracked_uniqueface_an4, tracked_uniqueface_an4_clean, 'rows');
% face_area_an4 = face_area_an4(mask);
% face_area_an5 = face_area_an4 + face_area_diff;
% dists_an4 = dists_an4(mask, :);
% dists_an5 = dists_an5(mask, :);
% dists = dists_an4;

% --------------------------- box plot ---------------------------
data = [face_area_diff, face_itg_abscurv_diff];
data_sigma3 = data(dists(:,1) < 5, :);
data_sigma5 = data(dists(:,2) < 5, :);
data_sigma7 = data(dists(:,3) < 5, :);
data_sigma9 = data(dists(:,4) < 5, :);
data_sigma27a = data(dists(:,5) < 5, :);
data_sigma27b = data(dists(:,6) < 5, :);
data_general = data(all(dists>5, 2), :);

labels = {['#(sigma3) = ', num2str(length(data_sigma3))], ['#(sigma5) = ', num2str(length(data_sigma5))], ...
        ['#(sigma7) = ', num2str(length(data_sigma7))], ['#(sigma9) = ', num2str(length(data_sigma9))], ...
        ['#(sigma27a) = ', num2str(length(data_sigma27a))], ['#(sigma27b) = ', num2str(length(data_sigma27b))],...
        ['#general = ', num2str(length(data_general))]};
 
idx_plot = 1;
area_diff = [data_sigma3(:,idx_plot); data_sigma5(:,idx_plot); data_sigma7(:,idx_plot); data_sigma9(:,idx_plot); ...
            data_sigma27a(:,idx_plot); data_sigma27b(:,idx_plot); data_general(:,idx_plot)];
groups = [ones(size(data_sigma3(:,1))); ones(size(data_sigma5(:,1)))*2; ...
            ones(size(data_sigma7(:,1)))*3; ones(size(data_sigma9(:,1)))*4; ...
            ones(size(data_sigma27a(:,1)))*5; ones(size(data_sigma27b(:,1)))*6; ...
            ones(size(data_general(:,1)))*7];
boxplot(area_diff,groups, 'Symbol', '');
% ylim([-650, 650])
set(gca,'XTickLabel', labels)
xtickangle(45)
line ([0, 8], [0, 0], 'Color', [0.5, 0.5, 0.5], 'LineStyle',':')
ylabel('\DeltaA')
% title('all tracked faces')
% title('no small faces')
title('full tracked faces')
print('dA_byMisorientation_box_full_tracked','-dtiff','-r300')

% disp([num2str(sum(face_area_an4)), ',  ', num2str(sum(face_area_an5)), ',  ', ...
%     num2str((sum(face_area_an5) - sum(face_area_an4)) / sum(face_area_an4))])
% disp([num2str((sum(face_area_an4(dists(:,1) < 5)))), ',  ', num2str((sum(face_area_an5(dists(:,1) < 5)))), ',  ', ...
%     num2str((sum(face_area_an5(dists(:,1) < 5)) - sum(face_area_an4(dists(:,1) < 5))) / sum(face_area_an4(dists(:,1) < 5)))])
%%
% --------------------------- cdf plot ---------------------------
figure()
h(1) = cdfplot(data_sigma3(:,1));
hold on
h(2) = cdfplot(data_sigma5(:,1));
h(3) = cdfplot(data_sigma7(:,1));
h(4) = cdfplot(data_sigma9(:,1));
h(5) = cdfplot(data_sigma27a(:,1));
h(6) = cdfplot(data_sigma27b(:,1));
h(7) = cdfplot(data_general(:,1));
set( h(:), 'LineWidth',2);
xlim([-1000, 1000])
legend(labels, 'Location', 'southeast')
xlabel('dA')
% title('all tracked faces')
% title('no small faces')
title('full tracked faces')
print('dA_byMisorientation_cdf_full_tracked','-dtiff','-r300')

% csvwrite('no_small_face_mask.txt', mask)

