function faceCentroids = findFaceCentorids(file, faces)
% ##########################################################################
% * Input
%     - faces([N,2]) and faceCorresp([N,1]) is returned by the function 'TrackFace'
% * Output
%     - faceCentroids = [N,3] = [x, y, z]
%         the data was in the same order as faces so the correspodence 
%         returned by 'TrackFace' can be applied directly. 
% * NOTE
%     - make sure samples are aligned by changing ORIGIN.
%     - the triangles with extreme curvature aren't excluded
% ##########################################################################
% ----------------------- load debug data -----------------------
% file = file_An5;
% faces = faces_An5;
% ---------------------------------------------------------------

FL = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
    % -- NOTE triNodes are indexes starting from zero 
triNodes = 1 + double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedTriList'))';
NCoords = double(h5read(file,'/DataContainers/TriangleDataContainer/_SIMPL_GEOMETRY/SharedVertexList'))';

% ##### filter bad data #####
mask = all(FL > 0, 2);
FL = FL(mask,:);
triNodes = triNodes(mask,:);

% ##### sort triangle data (node coords) by facelabel. #####
FL = sort(FL, 2);
data = [FL, zeros(length(FL), 1), triNodes];
data = sortrows(data);

% ##### prepare a holder the same length as faces: add ID then sort. #####
faces = sort(faces, 2);
faces = [(1:length(faces))', faces];
faces = sortrows(faces, [2,3]);
faceCentroids = [faces, zeros(length(faces), 3)];

% ##### loop list to calc faceCentorid #####
% ##### add result into the corresponding place in holder. Copy its equivalence if the face didn't exist #####
cnt = 0;
idx = 1;
tmp_centroid = zeros(1,3);
for i = 1:length(data)-1
    if data(i, 1) == data(i+1, 1) && data(i, 2) == data(i+1, 2)
        triCoords = NCoords(data(i, 4:6),:);
        tmp_centroid = tmp_centroid + sum(triCoords);
        cnt = cnt + 3;
    else
        triCoords = NCoords(data(i, 4:6),:);
        tmp_centroid = tmp_centroid + sum(triCoords);
        cnt = cnt + 3;
        if data(i,1) == faceCentroids(idx, 2) && data(i,2) == faceCentroids(idx, 3)
            faceCentroids(idx, 4:6) = tmp_centroid / cnt;
            faceCentroids(idx+1, 4:6) = faceCentroids(idx, 4:6);
            idx = idx + 2;
        end
        cnt = 0;
        tmp_centroid = zeros(1,3);
    end
end
% -- the last triangle/face wasn't really processed
triCoords = NCoords(data(i, 4:6),:);
tmp_centroid = tmp_centroid + sum(triCoords);
cnt = cnt + 3;
if data(i,1) == faceCentroids(idx, 2) && data(i,2) == faceCentroids(idx, 3)
    faceCentroids(idx, 4:6) = tmp_centroid / cnt;
    faceCentroids(idx+1, 4:6) = faceCentroids(idx, 4:6);
end 

faceCentroids = sortrows(faceCentroids, 1);
faceCentroids = faceCentroids(:,4:6);
end













