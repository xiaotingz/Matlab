function [faces_1_raw, faces_2_raw, corresp] = TrackFace(numNeigh_1, neighList_1, numNeigh_2, neighList_2, lookUp)
% ###########################################################################
% * Input
%     - lookUp = [N,2] = [ID_s1, ID_s2] 
%         is a table which contains only the tracked faces and was sorted in the order of state2.
% * Output
%     - faces = [N,2] = [labelA, labelB] are raw lists including all possible faces. Note just the ones that are tracked.
%     - corresp = [n,2] = [faceID_state1, faceID_state2]
% * Note the facesList is duplicated by [A,B] and [B,A]. This code only gives half correspondence. 
% ###########################################################################
% ------------------ load data for debug --------------------
% numNeigh_1 = numNeigh_An4;
% numNeigh_2 = numNeigh_An5;
% neighList_1 = neighList_An4;
% neighList_2 = neighList_An5;
% -----------------------------------------------------------

    % ##### Find Face Correspondence #####
    %	--- make a look up list that converts ID_An5 to ID_An4. ---
    %	---	if no correspondence, assign NaN ---
    %    *  lookUp_2to1 is made the following way:
    %      -  It was a table of size [N,2], in which N is the number of grains in state_2
    %      -  The first column is the grainID in state_1 and the second column is the grainID in state_2.
    %      -  The rows of the table was sorted so the second column is 1:N. 
    %               Note if a grain had ID_state_2=i but had no correspondence in state_1, the rows is [NaN, i]
    %      -  Now ID_state_2 is implicit. Take the first column to be lookUp_2to1, then it will give the ID_state_1 for grains in state_2
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

    faces_1_raw = [];
    for i = 1:length(numNeigh_1)
        neighbors = getNeighList(i, numNeigh_1, neighList_1);
        faces_1_raw = [faces_1_raw; [i*ones(size(neighbors)), neighbors]];
    end
    faces_1 = [(1:length(faces_1_raw)).', faces_1_raw];
    faces_1 = sortrows(faces_1,[2,3]);
    
    faces_2_raw = [];
    for i = 1:length(numNeigh_2)
        neighbors = getNeighList(i, numNeigh_2, neighList_2);
        faces_2_raw = [faces_2_raw; [i*ones(size(neighbors)), neighbors]];
    end
    
    faces_2_inID1 = [(1:length(faces_2_raw)).', lookUp_2to1(faces_2_raw)];
    mask = ~(any(isnan(faces_2_inID1),2));
    faces_2_inID1 = faces_2_inID1(mask,:);
    faces_2_inID1 = sortrows(faces_2_inID1,[2,3]);

% --- compare the common element thus track the faces, when FULL list ---
    length_1 = length(faces_1);
    length_2 = length(faces_2_inID1);
    i = 1;
    j = 1;
    corresp = [];
    while j <= length_1 && i <= length_2
        if faces_2_inID1(i,2) == faces_1(j,2)
            if faces_2_inID1(i,3) == faces_1(j,3)
                corresp = vertcat(corresp, [faces_1(j,1), faces_2_inID1(i,1)]);
                i = i + 1;
                j = j + 1;
            elseif faces_2_inID1(i,3) > faces_1(j,3)
                j = j + 1;
            else
                i = i + 1;
            end
        elseif faces_2_inID1(i,2) > faces_1(j,2)
            j = j + 1;
        else
            i = i + 1;
        end
    end
    corresp = sortrows(corresp, 2);
end