function b = calcMorawiecBMat(dg, n_1)
% ##########################################################################
% * Input
%     - dg = [3, 3]
%           Orientation matrix
%     - n_1 = [3, 1]
%           Plane normal, in reference frame of the first grain, NOT sample frame. 
% * Output
%     - b = [4, 4]
%           A matrix for grain boundary.
%           Note Yufeng's formula is different from Kryz's formula. Using Yufeng's.
% ##########################################################################
b = zeros(4,4);
b(1:3, 1:3) = dg;
b(1:3, 4) = n_1;
b(4, 1:3) = - n_1' * dg;
% b(4, 1:3) = - dg * n_1;
end




















