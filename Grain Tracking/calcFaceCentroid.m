function face_centroids = calcFaceCentroid(file, faces)
% ##########################################################################
% * Input
%     - faces([n,2]) and faceCorresp([n,1]) 
%         is returned by the function trackFace.m or trackUniqueFace.m
% * Output
%     - face_centroids = [n,1] 
%         the data was in the same order as faces so the correspodence 
%         returned by 'TrackFace' can be applied directly. 
% * NOTE
%     - this file is basically designed for anomaly detection. The
%     assumption is abnormal faces cluster together. 
% ##########################################################################
% % ----------------------- load debug data -----------------------
% file = '/Volumes/XIAOTING/Ni/An4new6_fixOrigin3_Hsmooth.dream3d';
% load('/Volumes/XIAOTING/Ni/190425_Hsmooth_geo_topo_an5crop2.mat', 'tracked_uniqueface_an4');
% faces = tracked_uniqueface_an4;
% % ---------------------------------------------------------------

% ----------------------- load data -----------------------
fl = h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
tri_nodes = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
node_coords = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';

mask = all(fl>0, 2);
fl = fl(mask, :);
tri_nodes = tri_nodes(mask, :);
fl = sort(fl, 2);
faces = sort(faces, 2);

% ----------------------- calculate centroid for each given face -----------------------
face_centroids = zeros(size(faces, 1), 3);
for i = 1:length(faces)
    mask = (fl(:, 1) == faces(i, 1) & fl(:, 2) == faces(i, 2));
    nodes = unique(tri_nodes(mask, :));
    coords = node_coords(nodes, :);
    face_centroids(i, :) = sum(coords) / length(coords);
end

end
% csvwrite('An4new6_fixOrigin3_Hsmooth_face_centroids.csv',face_centroids);





