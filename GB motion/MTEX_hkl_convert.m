cs = crystalSymmetry('m-3m');
axis = [-12,19,-4];
tmp = Miller(axis(1),axis(2),axis(3),cs,'hkl');
tmp2 = round(tmp);

axis_rounded = tmp2.hkl
