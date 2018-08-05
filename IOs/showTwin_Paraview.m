function showTwin_Paraview(file, Angle_thres, Axis_thres)
% ##################################################################
% Instruction
%   * This script is to help visualize grain faces associated with two
%        unique misorientations: 57° / [11-20] and 90° / [10-10] with HCP symmetry
%   * The comparison with the twin orientation is naively speciafied by
%        comparing the misorientation angleDifference and norm(axisAifference)
%         - Angles are all in degrees
%   * In order to run this script, one have run the following D3D filters
%         - 'Find Feature Average Orientations' (make sure to name average Euler Angles as?AvgEulerAngles)
%         - 'Generate Triangle Face Ids'. 
%                 !!! NOTE !!! --> this filter has to be run twice.
%                       The first time name the variable 'Feature Face Ids' to be 'TwinIds'.
%                       The second time, use 'Delete data' filter to delete
%                           'TriangleDataContainer/FaceFeatureData', then run 
%   * There are two dependency: which should be put in the same directory as this script
%           - generateHexSym.m
%           - EAtoG.m
%   * After run this script, read and write the file in D3D again to update the .Xdmf file
%   * In paraview, chose the filter FeatureFaceId
%        0 --- not designed twin
%        1 --- 57° / [11-20]
%        2 --- 90° / [10-10]
% ##################################################################
% ------------------------------------------------------
% Input for Debug
% file = ('/Users/xiaotingzhong/Desktop/Datas/HCP twin/labeled_mesh2_forParaview.dream3d');
% A_thres = 5;
% Axis_thres = 0.3;
% ------------------------------------------------------
TwinId = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/TwinId')';
FL_unique = h5read(file, '/DataContainers/TriangleDataContainer/FaceFeatureData/FaceLabels')';
FL_tris = h5read(file, '/DataContainers/TriangleDataContainer/FaceData/FaceLabels')';
EAs = h5read(file, '/DataContainers/ImageDataContainer/CellFeatureData/AvgEulerAngles')';
% goodMisA = h5read(file, '/DataContainers/ImageDataContainer/CellFeatureData/MisorientationList')';
EAs(1, :) = [];
EAs = rad2deg(EAs);

% --- the objective twins ---
twin_A1 = 57;
twin_Axis1 = [1,1,0];
twin_Axis1 = twin_Axis1/norm(twin_Axis1);
twin_A2 = 90;
twin_Axis2 = [2,1,0];
twin_Axis2 = twin_Axis2/norm(twin_Axis2);

% --- prepare symmetry operator ---
O = generateHexSym();
% --- convert Euler Angles to orientation matrix ---
Gs = zeros(3,3,length(EAs));
for i = 1:length(EAs)
    Gs(:,:,i) = EAtoG(EAs(i,:));
end

% --- get misorientation for all inner grain faces as axis angle ---
mask_innerFace = all(FL_unique > 0, 2);
FL_unique = FL_unique(mask_innerFace, :);
misA = zeros(length(FL_unique), 1);
misAxis = zeros(length(FL_unique), 3);
for i = 1:length(FL_unique)
    g1 = Gs(:, :, FL_unique(i,1));
    g2 = Gs(:, :, FL_unique(i,2));
    
    minA = 180;
    minAxis = [0,0,0];
    for j = 1:size(O,3)
        for k = 1:size(O,3)
            gg1 = O(:,:,j)*g1;
            gg2 = O(:,:,k)*g2;
            
            dg = gg1*gg2';
            tmpA = acosd(0.5*(trace(dg)-1));
            tmpAxis = [dg(2,3)-dg(3,2), dg(3,1)-dg(1,3), dg(1,2)-dg(2,1)];
            if tmpA < minA && all(tmpAxis > 0, 2)
                minA = tmpA;
                minAxis = tmpAxis;
            end
            
            dg = gg2*gg1';
            tmpA = acosd(0.5*(trace(dg)-1));
            tmpAxis = [dg(2,3)-dg(3,2), dg(3,1)-dg(1,3), dg(1,2)-dg(2,1)];
            if tmpA < minA && all(tmpAxis > 0, 2)
                minA = tmpA;
                minAxis = tmpAxis;
            end           
            
        end
    end
    misA(i) = minA;
    misAxis(i, :) = minAxis;
end

% --- normalize the misorientation axes to be unit vector ---
misAxis = misAxis ./ repmat(sqrt(sum(misAxis .* misAxis, 2)), 1,3);

% histogram(misA, 'normalization', 'probability')
% hold on
% histogram(goodMisA, 'normalization', 'probability')

diff_Axis1 = [misAxis(:,1) - twin_Axis1(1), misAxis(:,2) - twin_Axis1(2), misAxis(:,3) - twin_Axis1(3)];
diff_Axis1 = sum(diff_Axis1 .* diff_Axis1, 2);
mask_twin1 = (abs(misA - twin_A1) < Angle_thres  & diff_Axis1 < Axis_thres);
diff_Axis2 = [misAxis(:,1) - twin_Axis2(1), misAxis(:,2) - twin_Axis2(2), misAxis(:,3) - twin_Axis2(3)];
diff_Axis2 = sum(diff_Axis2 .* diff_Axis2, 2);
mask_twin2 = (abs(misA - twin_A2) < Angle_thres  & diff_Axis2 < Axis_thres);

Twin1s = FL_unique(mask_twin1, :);
Twin2s = FL_unique(mask_twin2, :);

mask_T1 = ismember(FL_tris, [Twin1s; [Twin1s(:,2), Twin1s(:,1)]], 'rows');
mask_T2 = ismember(FL_tris, [Twin2s; [Twin2s(:,2), Twin2s(:,1)]], 'rows');
mask_none = ~(mask_T1 | mask_T2);

TwinId(mask_none) = 0;
TwinId(mask_T1) = 1;
TwinId(mask_T2) = 2;

h5write(file, '/DataContainers/TriangleDataContainer/FaceData/TwinId', TwinId');

disp('FeatureFaceId=0 are the faces which are not twins.');
disp(['FeatureFaceId=1 are the 57°@[11-20] twins, for which there are   ', num2str(sum(mask_twin1)), '   faces']);
disp(['FeatureFaceId=2 are the 90°@[10-10] twins., for which there are   ', num2str(sum(mask_twin2)), '   faces']);
end

