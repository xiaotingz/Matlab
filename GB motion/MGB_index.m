clear

CD = 9;
MGB = xlsread('/Users/xiaotingzhong/Dropbox/Paper\ \&\ Documentation/Documentations/Cu_moved\ GB.xlsx','A2:D19');
MGB(:,1) = deg2rad(MGB(:,1));

O = CrysSym;
cnt = 1;
for k = 1:length(MGB)
    % matrix from Axis-Angle
    theta = MGB(k,1);
    u = MGB(k,2:4)/norm(MGB(k,2:4));
    delg = AAToG(theta, u);

    % apply symmetry on matrix and get box ID from BEAs
    for i = 1:24
        for j = 1:24
            delg_sym = O(:,:,i)*delg*O(:,:,j);
            EA_tmp = GToE(delg_sym);
            if (EA_tmp(1) < pi/2 && EA_tmp(2) < pi/2 && EA_tmp(3) < pi/2)
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
           
                MGB_indexFull(cnt,1) = c1 + CD*c2 + CD*CD*c3;
                MGB_indexFull(cnt,2) = k;
                cnt = cnt + 1;
            end
        end
    end
end
MGB_indexes = sortrows(MGB_indexFull,1);
MGB_indexes = unique(MGB_indexes,'rows');