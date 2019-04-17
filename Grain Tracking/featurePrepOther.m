% ##########################################################################
% * Notes
%     - Misorientaion distance is calculated from misorientation angle. 
%     - Distance to coherent twin is calculated as defined by Moraweic.
% ##########################################################################

load('/Users/xiaotingzhong/Documents/Matlab/Grain Tracking/data/181108_rfvec.mat', 'rfvecs_an4', 'rfvecs_an5');
load('/Users/xiaotingzhong/Documents/Matlab/Grain Tracking/data/181107_mig_piececorresp_comb',...
        'face_piecewise', 'tracked_uniqueface_an4')
file_an4 = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixOrigin2_smooth.dream3d');
    
% ------------------ identify broken faces ------------------
is_onepiece = ones(size(tracked_uniqueface_an4, 1), 1);
is_onepiece(face_piecewise) = 0;


% ------------------ identify twins ------------------
not_twin = ones(size(tracked_uniqueface_an4, 1), 1);
rfvec_twin = [1,1,1]/norm([1,1,1]) * tand(60/2);
eps_twin = 0.05;
for i = 1:length(rfvecs_an4)
    if norm(rfvecs_an4(i,:) - rfvec_twin) < eps_twin || norm(rfvecs_an5(i,:) - rfvec_twin) < eps_twin
        not_twin(i) = 0;
    end
end


% ------------------ distance from twins ------------------
dist_twin = vecnorm(rfvecs_an4 - rfvec_twin, 2, 2);





%% #################################### Write txt file ####################################
fileID = fopen('190408_features_otherinfo.txt','w');
% fprintf(fileID,'%s, %s, %s, %s\n','is_onepiece', 'not_twin', 'dist_twin', 'weighteddist_ctwin');
fprintf(fileID,'%s, %s, %s\n','is_onepiece', 'not_twin', 'dist_twin');
for i = 1:length(is_onepiece)
    fprintf(fileID, '%6d, %6d, %6.3f\n', ...
        is_onepiece(i), not_twin(i), dist_twin(i)) ;
end
fclose(fileID);

