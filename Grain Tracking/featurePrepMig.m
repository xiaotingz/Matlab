% ##########################################################################
% * Notes
%     - sign = 1    if face moving towards the grain with lower data_grain; 
%       sign = -1   if face moving toward the grain with higher data_grain;
%       sign = 0    if the two grains swapped data_grain level.
%     - Dependence: calcFaceToGrainCentroidDist.m
% ##########################################################################
load('/Users/xiaotingzhong/Documents/Matlab/Grain Tracking/projection/190404_proj_dists.mat')
mig_abs_nearest = zeros(size(mig_localnorm_nearest, 1), 1);
mig_abs_ot = zeros(size(mig_localnorm_nearest, 1), 1);
for i = 1:size(mig_localnorm_nearest, 1)
    cnt = 2;
    nearest = sum(mig_localnorm_nearest(i, 1:2));
    if mig_lr_nearest(i, 1) > 0 
        nearest = nearest + mig_lr_nearest(i, 1);
        cnt = cnt + 1;
    end
    if mig_svm_nearest(i, 1) > 0 
        nearest = nearest + mig_svm_nearest(i, 1);
        cnt = cnt + 1;
    end
    mig_abs_nearest(i) = nearest / cnt;
    
    cnt = 2;
    ot = sum(mig_localnorm_ot(i, 1:2));
    if mig_lr_ot(i, 1) > 0 
        ot = ot + mig_lr_ot(i, 1);
        cnt = cnt + 1;
    end
    if mig_svm_ot(i, 1) > 0 
        ot = ot + mig_svm_ot(i, 1);
        cnt = cnt + 1;
    end
    mig_abs_ot(i) = ot / cnt;
end
    

%% #################################### Write txt file ####################################
%%
fileID = fopen('190408_features_mig.txt','w');
fprintf(fileID,'%s, %s\n', 'mig_abs_nearest', 'mig_abs_ot');
for i = 1:length(mig_localnorm_nearest)
    fprintf(fileID, '%6.3f, %6.3f\n', mig_abs_nearest(i), mig_abs_ot(i));
end
fclose(fileID);


%%
fileID = fopen('190408_features_mig.txt','w');
fprintf(fileID,'%s, %s, %s, %s, %s, %s, %s, %s, %s\n', ...
    'move_left', 'svm_abs', 'lr_abs', 'normal_abs', 'pillar_abs', ...
                 'svm_sign', 'lr_sign', 'normal_sign', 'pillar_sign');
for i = 1:length(move_left)
    fprintf(fileID, '%6d, %6.3f\n',  move_left(i));
end
fclose(fileID);









