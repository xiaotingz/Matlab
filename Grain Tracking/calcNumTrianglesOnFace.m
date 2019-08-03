function num_tris = calcNumTrianglesOnFace(file, faces, eps_curv, eps_area, eps_min_ang)
% ##########################################################################
% * Input
%     - faces([n,2]) and faceCorresp([n,1]) 
%         is returned by the function trackFace.m or trackUniqueFace.m
%     - eps_
%         thresholds for good quality triangles. Especially important for data by Hsmooth.
% * Output
%     - num_tris = [n, 1]
%         number of triangles on the given face
% * NOTE
%     - This function is used with writeSeudoGBCDTris.m or calcFaceItgCurv.m
% ##########################################################################
% ----------------------- load debug data -----------------------
% file = file_an4;
% faces = faces_an4;
% ---------------------------------------------------------------
% -------- Read Data --------
facelabel = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels'))';
tri_curv =  abs(roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5))';
tri_area = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5)';
tri_min_ang = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceDihedralAngles'),-5)';

% -------- Clean Data --------
mask = all(facelabel > 0, 2) & abs(tri_curv) < eps_curv & tri_area < eps_area & tri_min_ang > eps_min_ang;
facelabel = facelabel(mask, :);
facelabel = sort(facelabel, 2);
faces = sort(faces, 2);

num_tris = zeros(size(faces, 1), 1);
for i = 1:size(faces, 1)
    mask = facelabel(:,1) == faces(i, 1) & facelabel(:,2) == faces(i, 2);
    num_tris(i) = sum(mask);
end


end


% ####################### Checks #######################
% see featurePrepGeoAndTopo.m, section 5.1















