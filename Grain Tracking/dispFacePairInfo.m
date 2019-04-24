function dispFacePairInfo(file_1, file_2, tracked_uniqueface_1, tracked_uniqueface_2, idx)
% ############################################################################
% * Input 
%     - idx, the pairID of tracked_uniqueface_
% ############################################################################
% ------------------ load data for debug --------------------
% file_1 = file_an4;
% file_2 = file_an5;
% tracked_uniqueface_1 = tracked_uniqueface_an4;
% tracked_uniqueface_2 = tracked_uniqueface_an5;
% % -----------------------------------------------------------
face_id_1 = h5read(file_1, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
face_id_2 = h5read(file_2, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
face_id_1(1,:) = [];
face_id_2(1,:) = [];
face_id_1 = sort(face_id_1, 2);
face_id_2 = sort(face_id_2, 2);


obj_label_1 = sort(tracked_uniqueface_1(idx, :));
[idx_1, ~] = find(face_id_1(:,1) == obj_label_1(1) & face_id_1(:,2) == obj_label_1(2));
obj_label_2 = sort(tracked_uniqueface_2(idx, :));
[idx_2, ~] = find((face_id_2(:,1) == obj_label_2(1) & face_id_2(:,2) == obj_label_2(2)) | (face_id_2(:,1) == obj_label_2(2) & face_id_2(:,2) == obj_label_2(1)));

% disp(' ');
% disp(['Pair ', num2str(idx)]);
% disp(['FaceLabel in an4:  [', num2str(face_id_1(idx_1, :)), ']', ]);
disp(['face_label in an4:  [', num2str(tracked_uniqueface_1(idx, :)), ']', ]);
disp(['face_id in an4:   ', num2str(idx_1)])
% disp(['FaceLabel in an5:  [', num2str(face_id_2(idx_2, :)), ']', ]);
disp(['face_label in an5:  [', num2str(tracked_uniqueface_2(idx, :)), ']', ]);
disp(['face_id in an5:   ', num2str(idx_2)])

end






