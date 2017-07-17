clc
clear
% read data
    % moved GBs
load('MGB_indexes');
    % The GBs in current map
FID=fopen('/Users/xiaotingzhong/Desktop/SEMs/053117_Cu/053117_hole/053117_bottom.txt');
datacell = textscan(FID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 11, 'Delimiter',',');
fclose(FID);
IDs = 1:length(datacell{1});
    % GBInMap[GBInMap_index, angle, axis, length]
GBInMap = [IDs.', datacell{7}, datacell{8},datacell{9},datacell{10},datacell{14}];
GBInMap = sortrows(GBInMap,2);
    % Get rid of the many sigma3s and the artificial duplicates 
i = 1;
while i <= length(GBInMap)-1
    angle = GBInMap(i,2);
    axis = abs(GBInMap(i,3:5)/norm(GBInMap(i,3:5)));
    IsSigma3 = (angle - 60)<3 && norm(axis - [1,1,1]/norm([1,1,1]))<0.05;
    short = GBInMap(i,6)<200;
    if IsSigma3 || short
        GBInMap(i,:) = [];
    else
        sameGB = (abs(GBInMap(i,2) - GBInMap(i+1,2)) < 1) && (norm(GBInMap(i,3:5)/norm(GBInMap(i,3:5)) - GBInMap(i+1,3:5)/norm(GBInMap(i+1,3:5))) < 0.05);
        if sameGB || short
            GBInMap(i+1,:) = [];
        end
        i = i + 1;
    end
end
    % convert the angles in degrees to radians
GBInMap(:,2) = deg2rad(GBInMap(:,2));

% get misorientation matrix
CD = 9.0;
O = CrysSym;
cnt = 1;

box_ind_full = zeros(length(GBInMap), 24*24);

for i = 1:length(GBInMap)
    angle = GBInMap(i,2);
    axis = GBInMap(i,3:5)/norm(GBInMap(i,3:5));
    delg = AAToG(angle, axis);    
    
    % apply symmetry to get symmetrically equivalent representations
    for j = 1:24
        for k = 1:24
            delg_sym = O(:,:,j)*delg*O(:,:,k).';
            EA_tmp = GToE(delg_sym);
            
            % get box_index in Euler Space
            if (EA_tmp(1) < pi/2) && (EA_tmp(2) < pi/2) && (EA_tmp(3) < pi/2)
                c1 = int32(floor(CD*2.0*EA_tmp(1)/pi) + 1);
                c2 = int32(floor(CD*cos(EA_tmp(2))) + 1);
                c3 = int32(floor(CD*2.0*EA_tmp(3)/pi) + 1);
                if c1 == CD+1
                    c1 = CD;
                elseif c2 == CD+1
                    c2 = CD;
                elseif c3 == CD+1
                    c3 = CD;
                end
                box_ind = c1 + CD*c2 + CD*CD*c3;
                box_ind_full(i, j + 24*(k-1)) = box_ind;
                
                
                % if the box_index is in the list of moved GBs, record this GB and its box_index
                % According to the box_index would be able to find out
                % which moved GB it is equivalent to.
                if ismember(box_ind, MGB_indexes(:,1))
                    matched_full(cnt,1) = GBInMap(i,1);
                    matched_full(cnt,2) = box_ind;
                    cnt = cnt+1;
                end
                
            end
        end
    end
end
matched = unique(matched_full,'rows');                
matchedGB_ind = unique(matched(:,1));

length(datacell{7})
length(GBInMap)
length(matchedGB_ind)