file_an4 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
file_an5 = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An5new6_smooth.dream3d';

% face_an4 = h5read(file_an4, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
% face_an5 = h5read(file_an5, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
% face_an4 = face_an4(all(face_an4>0, 2), :);
% face_an5 = face_an5(all(face_an5>0, 2), :);
% face_an4 = sortrows(face_an4);
% face_an5 = sortrows(face_an5);
% 
% 
% % [result_an4, ~, ~, ~] = findQuadNodes(file_an4);
% % tl_an4 = findTripleLines(file_an4);
% [num_corners_an4, num_edges_an4] = getFaceCharacter(face_an4, tl_an4, result_an4{1}, result_an4{2});
% 
% % [result_an5, ~, ~, ~] = findQuadNodes(file_an5);
% % tl_an5 = findTripleLines(file_an5);
% [num_corners_an5, num_edges_an5] = getFaceCharacter(face_an5, tl_an5, result_an5{1}, result_an5{2});

[faces_an5, num_neigh_face_an5, neigh_list_faceid_an5] = findFaceConnection(file_an5);
num_nnface_avgcorner_an5 = findFaceNNAvgCorner(file_an5, num_corners_an5);






