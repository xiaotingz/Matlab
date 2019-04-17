load('181028_migration.mat')
load('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/181025_mig_input.mat', 'face_to_calc')

mig_svm_abs = mig_svm_proj(:,2);
mig_normal_an4_abs = mig_normal_proj(:,3);
mig_normal_an5_abs = mig_normal_proj(:,4);
mig_normal_avg_abs = (mig_normal_proj(:,3) + mig_normal_proj(:,4))/2;

% figure
% scatter(mig_normal_an4_abs, mig_normal_an5_abs, 20, 'filled')
% xlabel('mig\_abs\_normal\_an4')
% ylabel('mig\_abs\_normal\_an5')
% set(gca,'fontsize',18)
% daspect([1 1 1])
% % print('mig_abs_normal_an4_&_an5', '-dpng','-r300')
% figure
% scatter(mig_normal_aveproj, mig_svm_abs, 20, 'filled')
% xlabel('mig\_abs\_normal\_avg')
% ylabel('mig\_abs\_normal\_svm')
% set(gca,'fontsize',18)
% daspect([1 1 1])
% print('Migration_CornerDiff', '-dpng','-r300')
% % print('mig_abs_normal_avg_&_mig_abs_svm', '-dpng','-r300')


sum(mig_normal_an4_abs > mig_normal_an5_abs);
sum(mig_normal_an4_abs < mig_normal_an5_abs);
sum(mig_svm_abs > mig_normal_an4_abs);

mask_diff_4_5 = (mig_normal_an4_abs > 3*mig_normal_an5_abs);
mask_diff_5_4 = (mig_normal_an5_abs > 3*mig_normal_an4_abs);
mask_diff_svm_avg = (mig_svm_abs > 5*mig_normal_avg_abs);
mask_svm_bad = (mig_svm_proj(:,3) == 1);
mask_svm_nan = (isnan(mig_svm_proj(:,1)));

idx_array = (1:length(mig_svm_proj))';
idx_diff_svm_avg = idx_array(mask_diff_svm_avg);
idx_diff_4_5 = idx_array(mask_diff_4_5);
idx_diff_5_4 = idx_array(mask_diff_5_4);
idx_svm_bad = idx_array(mask_svm_bad);
idx_svm_nan = idx_array(mask_svm_nan);


%% ##### Compare mig_sign & mig_abs #####
cdf1 = cdfplot(mig_svm_proj(:,2) - abs(mig_svm_proj(:,1)));
hold on
h2 = cdfplot(mig_normal_proj(:,3) - abs(mig_normal_proj(:,1)));
h3 = cdfplot(mig_normal_proj(:,4) - abs(mig_normal_proj(:,2)));
set(cdf1, 'LineWidth', 3);
set(h2, 'LineWidth', 3);
set(h3, 'LineWidth', 3);
set(gca,'FontSize', 18, 'LineWidth', 2)
xlim([0, 5])
xlabel('mig\_abs - abs(mig\_sign)')
legend('svm','normal\_an4','normal\_an5')


%% ##### Sign Problem #####
% idx_array = (1:length(mig_svm_proj))';
% mask_sign_diff = (mig_normal_proj(:,3) - abs(mig_normal_proj(:,1))) > 2 | ...
%     (mig_normal_proj(:,4) - abs(mig_normal_proj(:,2))) > 2;
% idx_sign_diff = idx_array(mask_sign_diff);

idx = idx_sign_diff(randi(312, 1))
x_to_y = X_to_Y_onepiece{idx};
face_node_info = getSingleFaceNodes(tracked_uniqueface_an4(face_to_calc(idx),:), tracked_uniqueface_an5(face_to_calc(idx),:));
visualizeFace(face_node_info, x_to_y)
mig_normal_proj(idx, :)
mig_svm_proj(idx, :)

% plotSingleFaceWithNormal(file_an4, tracked_uniqueface_an4(face_to_calc(idx),:), 1)
% hold on
% plotSingleFaceWithNormal(file_an5, tracked_uniqueface_an5(face_to_calc(idx),:), 1)


%% ##### Migration of Twins #####
% [rfvecs_an4]  = getFaceRFvecs(file_an4, tracked_uniqueface_an4);
% [rfvecs_an5]  = getFaceRFvecs(file_an5, tracked_uniqueface_an5);
% 
axis = [1,1,1];
axis = axis/norm(axis);
angle = 60;
rfvec_ctwin = axis*tand(angle/2);
rfvec_ctwin = repmat(rfvec_ctwin, length(rfvecs_an4), 1);

mask_an4_twin = vecnorm(rfvecs_an4 - rfvec_ctwin, 2, 2) < 0.06;
mask_an5_twin = vecnorm(rfvecs_an5 - rfvec_ctwin, 2, 2) < 0.06;
if sum(mask_an4_twin == mask_an5_twin) ~= length(mask_an4_twin)
    warning('twins are not the same in the two states')
else 
    mask_twin = mask_an4_twin;
    clear mask_an4_ctwin mask_an5_ctwin axis angle 
end


mig_svm_twin = mig_svm_proj(mask_twin, :);
mig_svm_not_twin = mig_svm_proj(~mask_twin, :);

mig_normal_twin = (mig_normal_proj(mask_twin, 3) + mig_normal_proj(mask_twin, 4))/2;
mig_normal_not_twin = (mig_normal_proj(~mask_twin, 3) + mig_normal_proj(~mask_twin, 4))/2;

bin_edges = 0:0.5:50;

figure
% histogram(mig_svm_not_twin(:,2), bin_edges)
% histogram(mig_svm_twin(:,2), bin_edges)
cdf1 = cdfplot(mig_svm_not_twin(:,2));
hold on
cdf2 = cdfplot(mig_svm_twin(:,2));
cdf3 = cdfplot(mig_normal_not_twin);
cdf4 = cdfplot(mig_normal_twin);
cdf1.LineWidth = 3; cdf1.LineStyle = '--';
cdf2.LineWidth = 3; cdf2.LineStyle = '--';
cdf3.LineWidth = 3;
cdf4.LineWidth = 3;
xlabel('Migration Distance, \mum', 'FontSize', 20)
ylabel('Fraction of Data', 'FontSize', 20)
set(gca,'FontSize', 18, 'LineWidth', 2)
legend('Proj\_svm\_notTwin', 'Proj\_svm\_Twin', 'Proj\_normal\_notTwin', 'Proj\_normal\_Twin',...
    'Location', 'southeast')
xlim([0, 20])
print('migration_PieceCorrespAllData_OptCover_CDF','-dpng','-r300')

%% ##### Mobile Twins #####
% load('180822_FaceCorresp.mat')

mask_mobile_twin = (mask_twin & (mig_normal_proj(:,3) > 5));
idx_array = (1:length(mig_normal_proj))';
idx_mobile_twin = idx_array(mask_mobile_twin);

idx = idx_mobile_twin(randi(length(idx_mobile_twin), 1))
x_to_y = X_to_Y_piece_corresp{idx};
if iscell(x_to_y)
    x_to_y = cell2mat(X_to_Y_piece_corresp{idx});
end
face_node_info = getSingleFaceNodes(tracked_uniqueface_an4(idx,:), tracked_uniqueface_an5(idx,:));
visualizeFace(face_node_info, x_to_y)
mig_normal_proj(idx, :)
mig_svm_proj(idx, :)


% file_name = ['mobile_twin_pair_', num2str(face_to_calc(idx))];





