% The old data which don't have the back info
MGB_old = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','C2:F19');
% The new data which have the back info
MGB_new = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','C22:F32');
cs = crystalSymmetry('m-3m');


old_comp = zeros(size(MGB_old,1));
new_comp = zeros(size(MGB_new,1));

for i = 1:size(MGB_old,1)
    for j = 1:size(MGB_old,1)
        v1 = vector3d(MGB_old(i,2:4));
        o1 = rotation('axis',v1,'angle',MGB_old(i,1)*degree, cs);
        v2 = vector3d(MGB_old(j,2:4));
        o2 = rotation('axis',v2,'angle',MGB_old(j,1)*degree, cs);
        old_comp(i,j) = angle(o1,o2);
    end
end

v1 = vector3d(-1,1,1);   
o1 = rotation('axis',v1,'angle',60*degree, cs);
v2 = vector3d(1,1,1); 
o2 = rotation('axis',v2,'angle',60*degree, cs);
rad2deg(angle(o1,o2))