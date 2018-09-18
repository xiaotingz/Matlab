% function result = isTwin(obj_facelabel_an4, orientation_an4, obj_facelabel_an5, orienation_an5)
% ##################################################################
% * Input
%     - obj_facelabel_ = [2, 1], label of the objective grain face
%     - orientation_ = [n, 3], orientations of the grains
% * Output
%     - result = 0 | 1, 0 --- not twin, 1 --- twin
% * Note
%     - !!! Make sure the first row of EA has been deleted !!!
% ##################################################################
% ------------------------------------------------------
% Input for Debug
% file = ('/Users/xiaotingzhong/Desktop/Datas/HCP twin/labeled_mesh2_forParaview.dream3d');
% A_thres = 5;
% Axis_thres = 0.3;
% ------------------------------------------------------

EAs = rad2deg(EAs);

% --- the objective twins ---
twin_ang = 60;
twin_axis = [1,1,1];
twin_axis = twin_axis/norm(twin_axis);

% --- prepare symmetry operator ---
O = generateHexSym();
% --- convert Euler Angles to orientation matrix ---
G = zeros(3,3,length(EAs));
for i = 1:length(EAs)
    G = EAtoG(EAs(i,:));


% --- get misorientation for all inner grain faces as axis angle ---
mask_innerFace = all(FL_unique > 0, 2);
FL_unique = FL_unique(mask_innerFace, :);
misA = zeros(length(FL_unique), 1);
misAxis = zeros(length(FL_unique), 3);
for i = 1:length(FL_unique)
    g1 = G(:, :, FL_unique(i,1));
    g2 = G(:, :, FL_unique(i,2));
    
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
diff_Axis2 = [misAxis(:,1) - twin_axis(1), misAxis(:,2) - twin_axis(2), misAxis(:,3) - twin_axis(3)];
diff_Axis2 = sum(diff_Axis2 .* diff_Axis2, 2);
mask_twin2 = (abs(misA - twin_ang) < Angle_thres  & diff_Axis2 < Axis_thres);

Twin1s = FL_unique(mask_twin1, :);
Twin2s = FL_unique(mask_twin2, :);

mask_T1 = ismember(FL_tris, [Twin1s; [Twin1s(:,2), Twin1s(:,1)]], 'rows');
mask_T2 = ismember(FL_tris, [Twin2s; [Twin2s(:,2), Twin2s(:,1)]], 'rows');
mask_none = ~(mask_T1 | mask_T2);

TwinId(mask_none) = 0;
TwinId(mask_T1) = 1;
TwinId(mask_T2) = 2;

