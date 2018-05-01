function [RFvecs, misAs]  = getFaceRFvecs(file, faces)
% ###########################################################################
% get the Rodrigues vector associated with each grain face
% * Output
%     - RFvecs = [N,3], in the same order as faces. 
% ###########################################################################
% ------------------ load data for debug --------------------
% file = file_An4;
% faces = faces_An4;
% -----------------------------------------------------------
    EAs = roundn(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles'),-5).';
    misAs = roundn(double(h5read(file,'/DataContainers/ImageDataContainer/CellFeatureData/MisorientationList')),-5);
    EAs(1,:) = [];
    %   ### 'EAtoG' and 'dgInFZ' takes data as degrees ### 
    EAs = rad2deg(EAs);
    O = CrysSym();

    %   ### RF caculation is expensive, for each unique face, do it once ###
    faces = [(1:length(faces))', faces];
    faces(:,2:3) = sort(faces(:,2:3), 2);
    faces = sortrows(faces, [2,3]);

    % tmp = zeros(length(faces),1);
    RFvecs = zeros(length(faces),4);
    RFvecs(:,1) = faces(:,1);
    for i = 1:2:length(faces)
        g1 = EAtoG(EAs(faces(i,2),:));
        g2 = EAtoG(EAs(faces(i,3),:));
        [RFvecs(i,2:4), ~] = dgInFZ(g1, g2, O);
        RFvecs(i+1,2:4) = RFvecs(i,2:4);
    %     tmp(i+1) = tmp(i);
    end
    %   ### sort RFvecs back into the order of faces ###
    RFvecs = sortrows(RFvecs,1);
    RFvecs = RFvecs(:,2:4);

    % %   ### validate dgInFZ by check the returned misA(named tmp) with the misA given by D3D ###
    % tmp2 = [RFvecs(:,1), tmp];
    % tmp2 = sortrows(tmp2,1);
    % sum(roundn(tmp2(:,2),-2)==roundn(misAs,-2))
end





