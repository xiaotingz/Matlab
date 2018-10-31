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







