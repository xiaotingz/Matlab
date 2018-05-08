function UFcorresp = TrackUniqueFace(UF_1, UF_2, lookUp)
% ###########################################################################
% * Very similar to TrackFace, except that the input is the uniqueFaces.
% * Input
%     - lookUp = [N,2] = [ID_s1, ID_s2] 
%         is a table which contains only the tracked faces and was sorted in the order of state2.
%     - UF_1,2
%         output of makeFaceCoords, labels of the unique faces. 
% * Output
%     - UFcorresp = [faceID_state1, faceID_state2]
%         this correspondence can be applied to faceCoords directly.
% ###########################################################################
% ------------------ load data for debug --------------------
% UF_1 = UFlabels_An4;
% UF_2 = UFlabels_An5;
% -----------------------------------------------------------
    lookUp = sortrows(lookUp, 2);  %  somehow the look up table is not really sorted well 
    lookUp_2to1 = zeros(max(lookUp(:,2)),2);
    %   --- implicitly, grainSizeDiff is for the grains in state_An5 ---
    idx = 1;
    for i = 1 : max(lookUp(:,2))
        if lookUp(idx,2) == i
            lookUp_2to1(i,1) = lookUp(idx,1);
            lookUp_2to1(i,2) = lookUp(idx,2);
            idx = idx + 1;
        else
            lookUp_2to1(i,1) = NaN;
            lookUp_2to1(i,2) = i;
        end
    end
    
    UF_1 = [(1:length(UF_1)).', UF_1];
    UF_2_inID1 = [(1:length(UF_2)).', lookUp_2to1(UF_2)];
    mask = ~(any(isnan(UF_2_inID1),2));
    UF_2_inID1 = UF_2_inID1(mask,:);
    UF_2_inID1(:,2:3) = sort(UF_2_inID1(:,2:3),2);
    UF_2_inID1 = sortrows(UF_2_inID1,[2,3]);

% --- compare the common element thus track the faces, when FULL list ---
    length_1 = length(UF_1);
    length_2 = length(UF_2_inID1);
    i = 1;
    j = 1;
    UFcorresp = [];
    while j <= length_1 && i <= length_2
        if UF_2_inID1(i,2) == UF_1(j,2)
            if UF_2_inID1(i,3) == UF_1(j,3)
                UFcorresp = vertcat(UFcorresp, [UF_1(j,1), UF_2_inID1(i,1)]);
                i = i + 1;
                j = j + 1;
            elseif UF_2_inID1(i,3) > UF_1(j,3)
                j = j + 1;
            else
                i = i + 1;
            end
        elseif UF_2_inID1(i,2) > UF_1(j,2)
            j = j + 1;
        else
            i = i + 1;
        end
    end
    UFcorresp = sortrows(UFcorresp, 2);
end
