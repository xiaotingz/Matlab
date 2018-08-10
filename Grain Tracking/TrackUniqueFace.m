function unique_facelabel_corresp = trackUniqueFace(unique_facelabel_1, unique_facelabel_2, look_up_table)
% ###########################################################################
% * Very similar to TrackFace, except that the input is the uniqueFaces.
% * Input
%     - lookUp = [N,2] = [ID_s1, ID_s2] 
%         is a table which contains only the tracked faces and was sorted in the order of state2.
%     - unique_facelabel_1,2
%         output of makeFaceCoords, labels of the unique faces. 
% * Output
%     - UFcorresp = [faceID_state1, faceID_state2]
%         this correspondence can be applied to faceCoords directly.
% ###########################################################################
% ------------------ load data for debug --------------------
% UF_1 = UFlabels_An4;
% UF_2 = UFlabels_An5;
% -----------------------------------------------------------
    look_up_table = sortrows(look_up_table, 2);  %  somehow the look up table is not really sorted well 
    look_up_2to1 = zeros(max(look_up_table(:,2)),2);
    %   --- implicitly, grainSizeDiff is for the grains in state_An5 ---
    idx = 1;
    for i = 1 : max(look_up_table(:,2))
        if look_up_table(idx,2) == i
            look_up_2to1(i,1) = look_up_table(idx,1);
            look_up_2to1(i,2) = look_up_table(idx,2);
            idx = idx + 1;
        else
            look_up_2to1(i,1) = NaN;
            look_up_2to1(i,2) = i;
        end
    end
    
    unique_facelabel_1 = [(1:length(unique_facelabel_1)).', unique_facelabel_1];
    unique_facelabel_2in1 = [(1:length(unique_facelabel_2)).', look_up_2to1(unique_facelabel_2)];
    mask = ~(any(isnan(unique_facelabel_2in1),2));
    unique_facelabel_2in1 = unique_facelabel_2in1(mask,:);
    unique_facelabel_2in1(:,2:3) = sort(unique_facelabel_2in1(:,2:3),2);
    unique_facelabel_2in1 = sortrows(unique_facelabel_2in1,[2,3]);

% --- compare the common element thus track the faces, when FULL list ---
    length_1 = length(unique_facelabel_1);
    length_2 = length(unique_facelabel_2in1);
    i = 1;
    j = 1;
    unique_facelabel_corresp = [];
    while j <= length_1 && i <= length_2
        if unique_facelabel_2in1(i,2) == unique_facelabel_1(j,2)
            if unique_facelabel_2in1(i,3) == unique_facelabel_1(j,3)
                unique_facelabel_corresp = vertcat(unique_facelabel_corresp, [unique_facelabel_1(j,1), unique_facelabel_2in1(i,1)]);
                i = i + 1;
                j = j + 1;
            elseif unique_facelabel_2in1(i,3) > unique_facelabel_1(j,3)
                j = j + 1;
            else
                i = i + 1;
            end
        elseif unique_facelabel_2in1(i,2) > unique_facelabel_1(j,2)
            j = j + 1;
        else
            i = i + 1;
        end
    end
    unique_facelabel_corresp = sortrows(unique_facelabel_corresp, 2);
end
