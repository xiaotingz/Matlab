% ##### plot the centroids #####
i = randi(5878, 1);
idx = face_onepiece(i);

disp(' ');
disp(['Pair ', num2str(idx), '  ;  one_piece_id = ', num2str(i)]);
dispFacePairInfo(file_an4, file_an5, tracked_uniqueface_an4, tracked_uniqueface_an5, idx)
disp(['mig_sign_norm_nearest = ',num2str(mig_localnorm_nearest(i, 3:4))])
disp(['mig_sign_norm_ot = ',num2str(mig_localnorm_ot(i, 3:4))])

figure
plotCentroids(file_an4, file_an5, tracked_uniqueface_an4, tracked_uniqueface_an5, idx)
% if move_left(idx) > 0
%     disp('move left')
%     title('move left', 'fontSize', 18)
% elseif move_left(idx) < 0
%     disp('move right')
%     title('move right', 'fontSize', 18)
% else
%     disp('stable')
%     title('stable', 'fontSize', 18)
% end

%% ##### save plot #####
f_name = ['centroids_pair_', num2str(idx)];
print(f_name, '-dtiff', '-r300')



