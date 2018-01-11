clear
clc
% load data
load('Cu_GBseed.mat')

FID=fopen('/Users/xiaotingzhong/Desktop/SEMs/091217_Cu/091217_fullEBSD/091217_8_top_GB.txt');
datacell = textscan(FID, '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f', 'HeaderLines', 11, 'Delimiter',',');
fclose(FID);
IDs = 1:length(datacell{1});
%     GBInMap[GBInMap_index, angle, axis, length]
InMapGB = [IDs.', datacell{7}, datacell{8},datacell{9},datacell{10},datacell{14}];
InMapGB = sortrows(InMapGB,2);
%     The threshold for telling the same GB



% Get rid of the many sigma3s and the artificial duplicates 
i = 1;
while i <= size(InMapGB,1)-1
    angle = InMapGB(i,2);
    axis = abs(InMapGB(i,3:5)/norm(InMapGB(i,3:5)));
    IsSigma3 = (angle - 60)<3 && norm(axis - [1,1,1]/norm([1,1,1]))<0.05;
    short = InMapGB(i,6)<150;
    if IsSigma3 || short
        InMapGB(i,:) = [];
    else
        sameGB = (abs(InMapGB(i,2) - InMapGB(i+1,2)) < 3) && (norm(InMapGB(i,3:5)/norm(InMapGB(i,3:5)) - InMapGB(i+1,3:5)/norm(InMapGB(i+1,3:5))) < 0.05);
        if sameGB || short
            InMapGB(i+1,:) = [];
        end
        i = i + 1;
    end
end

% For the moved GBs, construct the matrix from the axis angle 
InMapGB_mat = zeros(3,3,size(InMapGB,1));
for i = 1:size(InMapGB,1)
    angle = deg2rad(InMapGB(i,2));
    axis = InMapGB(i,3:5)/norm(InMapGB(i,3:5));
    InMapGB_mat(:,:,i) = AAToG(angle, axis); 
end
% For the moved GBs, transfer their matrix into FZ
O = CrysSym;
for i = 1:size(InMapGB_mat,3)
    g1 = InMapGB_mat(:,:,i);
    disA = pi/2;
    for j = 1:24
        for k = 1:24
            gg1 = O(:,:,j) * g1;
            gg2 = O(:,:,k);

            dg = gg1 * gg2.';
            misA = acos(0.5*(trace(dg)-1));
            RFvec = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA)) *tan(misA/2);       
            axisInFZ = ( (sum(RFvec >= 0) == 3) && (RFvec(1) <= (sqrt(2) - 1)) && ...
                        (RFvec(1) >= RFvec(2)) && (RFvec(2) >= RFvec(3)) && (sum(RFvec) <= 1));
            if ( misA < disA && axisInFZ )
                InMapGB_mat(:,:,i) = dg;
                disA = misA;
            end
            dg = gg2 * gg1.';
            misA = acos(0.5*(trace(dg)-1));
            RFvec = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA)) *tan(misA/2);
            axisInFZ = ( (sum(RFvec >= 0) == 3) && (RFvec(1) <= (sqrt(2) - 1)) && ...
                        (RFvec(1) >= RFvec(2)) && (RFvec(2) >= RFvec(3)) && (sum(RFvec) <= 1));
            if ( misA < disA && axisInFZ )
                InMapGB_mat(:,:,i) = dg;
                disA = misA;
            end       
        end
    end
end



% Compare the filtered inMapGB to the GB_seeds
sameGB_thres = 5;

GB_comp = zeros(size(InMapGB_mat,3),size(GB_seed,3));
matchedGB = [];
cnt = 1;
for i = 1:size(InMapGB_mat,3)
    for j = 1:size(GB_seed,3)
        disA = 100000;
        MatchID = 0;
        for k = 1:24
            for l = 1:24
                gg1 = O(:,:,k) * InMapGB_mat(:,:,i);
                gg2 = O(:,:,l) * GB_seed(:,:,j);
                delg = gg1*gg2.';
                angle = acosd(0.5*(trace(delg)-1));
                if angle < disA
                    disA = angle;
%                     disAxis = [(delg(2,3)-delg(3,2)),(delg(3,1)-delg(1,3)),(delg(1,2)-delg(2,1))] /(2*sin(disA));
                end
                delg = gg2*gg1.';
                angle = acosd(0.5*(trace(delg)-1));
                if angle < disA
                    disA = angle;
%                     disAxis = [(delg(2,3)-delg(3,2)),(delg(3,1)-delg(1,3)),(delg(1,2)-delg(2,1))] /(2*sin(disA));
                end
            end
        end
        GB_comp(i,j) = disA;
        if disA < sameGB_thres
            matchedGB(cnt, 1) = InMapGB(i,1);
            matchedGB(cnt, 2) = j;
            matchedGB(cnt, 3) = disA;
            matchedGB(cnt, 4) = InMapGB(i,2);
            matchedGB(cnt, 5:7) = InMapGB(i,3:5);
            matchedGB(cnt, 8) = InMapGB(i,6);
            cnt = cnt + 1;
        end
    end
end

InMapGB_ID = matchedGB(:,1);
SeedGB_ID = matchedGB(:,2);
matched_disA = matchedGB(:,3);
InMapGB_misA = matchedGB(:, 4);
InMapGB_misAxis = matchedGB(:, 5:7);
InMapGB_length = matchedGB(:, 8);
matchGB_table = table(InMapGB_ID, SeedGB_ID, matched_disA, InMapGB_misA, InMapGB_misAxis, InMapGB_length)


