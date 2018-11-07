% function idx_ctwins = getCoherentTwin(file, twin_faces)
% ###########################################################################
% * Inputs
%     - twin_faces = [n, 2]
%         facelabels for twins, returned probably by getFaceRFvecs.m
% * Notes
%     - This scripts take in the face_labels of twins and checks the
%     normals of the twins to see if they are coherent.
% ###########################################################################
% ------------------ load data for debug --------------------
% file = ('/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d');
% obj_facelabel = [332, 2275];
% reverse_winding = 1;
% -----------------------------------------------------------
