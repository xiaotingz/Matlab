% ##########################################################################
% * Notes
%     - This file calculates 
% ##########################################################################

% ------------------ identify twins and broken faces ------------------
load('181108_rfvec', 'rfvecs_an4', 'rfvecs_an5');
load('181107_mig_piececorresp_comb','face_piecewise')

is_onepiece = ones(size(tracked_uniqueface_an4, 1), 1);
is_onepiece(face_piecewise) = 0;

not_twin = ones(size(tracked_uniqueface_an4, 1), 1);
rfvec_twin = [1,1,1]/norm([1,1,1]) * tand(60/2);
eps_twin = 0.05;
for i = 1:length(rfvecs_an4)
    if norm(rfvecs_an4(i,:) - rfvec_twin) < eps_twin || norm(rfvecs_an5(i,:) - rfvec_twin) < eps_twin
        not_twin(i) = 0;
    end
end





%% #################################### Write txt file ####################################
fileID = fopen('190408_features_otherinfo.txt','w');
fprintf(fileID,'%s ,%s\n','is_onepiece', 'not_twin', 'dist_twin', 'weighteddist_ctwin');
for i = 1:length(move_left)
    fprintf(fileID, '%6d, %6d\n', ...
        is_onepiece(i), not_twin(i)) ;
end
fclose(fileID);

