function data_face = calcFaceCurvature(data_raw)
% -- calc integral curvature for each face
% -- get rid of bad datas
    % 1.facelabel <= 0; 2.extreme curvature value; 3.NaN(sometimes)
% ----------------------- load debug data -----------------------
% file = file_An4;
% FL = double(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceLabels')).';
% triCurves =  roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/MeanCurvatures'),-5).';
% triAreas = roundn(h5read(file,'/DataContainers/TriangleDataContainer/FaceData/FaceAreas'),-5).';
% data_raw = [FL.'; triCurves.'; triAreas.'];
% ---------------------------------------------------------------

% -- prepare for face curvature calculation: sort the data by facelabel
bool1 = data_raw(1,:) > 0 & data_raw(2,:) > 0 & ~isnan(data_raw(3,:)) & data_raw(3,:) < 100 & data_raw(3,:) > -100;
data_cleared  = data_raw(:,bool1);
data_sorted = sortrows(data_cleared.');
clear bool1 data_raw data_cleared

% -- calculate face curvature
tmp5 = unique(data_sorted(:,1:2),'rows');
data_face = zeros(length(tmp5),4);
j = 0;
k = 1;
counter = 2;
s_area = data_sorted(j+counter-1,4);
s_cur = data_sorted(j+counter-1,3)*data_sorted(j+counter-1,4);
for i = 2:length(data_sorted)
%   -- if next trianlge still in the same face
    if data_sorted(j+counter-1,1) == data_sorted(j+counter,1) && data_sorted(j+counter-1,2) == data_sorted(j+counter,2)
       s_cur = s_cur + data_sorted(j+counter,3)*data_sorted(j+counter,4);
       s_area = s_area + data_sorted(j+counter,4);
       counter = counter + 1;
%   -- if next trianlge in another face: record former one and restart counting
    else
        data_face(k,1) = data_sorted(j+counter-1,1);
        data_face(k,2) = data_sorted(j+counter-1,2);
        data_face(k,3) = s_area;
        data_face(k,4) = s_cur/s_area;
        k = k+1;   % check k = length(x)
        j = j + counter -1;
        counter = 2;
        s_area = data_sorted(j+counter-1,4);
        s_cur = data_sorted(j+counter-1,3)*data_sorted(j+counter-1,4);
    end
end
% -- looped all but no comparison to execute else for the last face ==> write manually
data_face(k,1) = data_sorted(j+counter-1,1);
data_face(k,2) = data_sorted(j+counter-1,2);
data_face(k,3) = s_area;
data_face(k,4) = s_cur/s_area;

end
