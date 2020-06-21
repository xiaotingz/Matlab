function rf_miso = calcMisorientationAsRFvec(obj_faces, file)
% ##########################################################################
% * Input
%     - obj_faces = [n , 2]
%           n = #obj_faces
% * Output
%     - rf_miso = [n, 3]
%           Misorientation in Rodrigues Vectors
% * Note
%     - Dependency: EAtoG.m, dgInFZ.
%     - For debug data, see OneNote/Readings/Notations/CLSs
% ##########################################################################
% ----------------------- load debug data -----------------------
% file = '/Users/xiaotingzhong/Desktop/Datas/Ni_an4_5/An4new6_fixedOrigin_smooth.dream3d';
% obj_faces = obj_faces_an4;
% ---------------------------------------------------------------
ea = rad2deg(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles'))';
ea(1,:) = [];
O = CrysSym;

rf_miso = zeros(size(obj_faces, 1), 3);
for i = 1:length(obj_faces)
    g_1 = EAtoG(ea(obj_faces(i,1), :));
    g_2 = EAtoG(ea(obj_faces(i,2), :));
    [rf_miso(i, :), ~, ~] = dgInFZ(g_1, g_2, O);
end

end







