% % Most interesting GB
% MGB = zeros(6,4);
% MGB(1,:) = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','C3:F3');
% MGB(2,:) = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','C6:F6');
% MGB(3,:) = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','C13:F13');
% MGB(4,:) = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','C18:F18');
% MGB(5,:) = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','C23:F23');
% MGB(6,:) = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','C28:F28');

clear

% The old data which don't have the back info
% old data
% MGB = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','C2:F19');
% new data
MGB = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','C22:F32');
% convert misorientations to matrix
M_matrix = zeros(3,3,size(MGB,1));
for i = 1:size(MGB,1)
    angle = deg2rad(MGB(i,1));
    axis = MGB(i,2:4)/norm(MGB(i,2:4));
    M_matrix(:,:,i) = AAToG(angle, axis); 
end


%% convert the misorientation into FZ
O = CrysSym;
for i = 1:size(M_matrix,3)
    g1 = M_matrix(:,:,i);
    disA = pi/2;
    for j = 1:24
        for k = 1:24
            gg1 = O(:,:,j) * g1;
            gg2 = O(:,:,k);

            dg = gg1 * gg2.';
            misA = acos(0.5*(trace(dg)-1));
            RFvec = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA)) *tan(misA/2);
    %         if ( misA < disA && sum(misAxis > 0) == 3 )        
            axisInFZ = ( (sum(RFvec >= 0) == 3) && (RFvec(1) <= (sqrt(2) - 1)) && ...
                        (RFvec(1) >= RFvec(2)) && (RFvec(2) >= RFvec(3)) && (sum(RFvec) <= 1));
            if ( misA < disA && axisInFZ )
                M_matrix(:,:,i) = dg;
            end
            dg = gg2 * gg1.';
            misA = acos(0.5*(trace(dg)-1));
            RFvec = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA)) *tan(misA/2);
    %         if ( misA < disA && sum(misAxis > 0) == 3 )
            axisInFZ = ( (sum(RFvec >= 0) == 3) && (RFvec(1) <= (sqrt(2) - 1)) && ...
                        (RFvec(1) >= RFvec(2)) && (RFvec(2) >= RFvec(3)) && (sum(RFvec) <= 1));
            if ( misA < disA && axisInFZ )
                M_matrix(:,:,i) = dg;
            end       
        end
    end
end

            
% calculate misA between misAs
M_comp = zeros(size(M_matrix,3));


for i = 1:size(M_matrix,3)
    for j = 1:size(M_matrix,3)
%         GG = M_matrix(:,:,i) * M_matrix(:,:,j).';
        minA = 100000;
        for k = 1:24
            for l = 1:24
                gg1 = O(:,:,k) * M_matrix(:,:,i);
                gg2 = O(:,:,l) * M_matrix(:,:,j);
%                 delg_sym = O(:,:,k)*delg*O(:,:,l).';
                delg = gg1*gg2.';
                angle = acosd(0.5*(trace(delg)-1));
                if angle < minA
                    minA = angle;
                end
                delg = gg2*gg1.';
                angle = acosd(0.5*(trace(delg)-1));
                if angle < minA
                    minA = angle;
                end
            end
        end
        M_comp(i,j) = minA;
    end
end
M_comp = round(M_comp,3);

lowerMat = tril(M_comp);

imagesc(lowerMat);
colormap jet;
colorbar
