clc

[mig_sign, mig_abs, bad_fit, features_lr]  = calcFaceMigByLinearRegPlaneProj(features, x_to_y, 0.1);
sum(mig_sign) / length(mig_sign)
sum(mig_abs) / length(mig_abs)


[mig_sign, mig_abs, bad_fit, features_svm]  = calcFaceMigBySVMPlaneProj(features, x_to_y, 0.1);
sum(mig_sign) / length(mig_sign)
sum(mig_abs) / length(mig_abs)



