% #################################################################
% MATLAB optimization example
%   https://www.mathworks.com/help/optim/ug/example-linear-programming-via-problem.html
% #################################################################

% ##### Variables #####
P1 = optimvar('P1','LowerBound',2500,'UpperBound',6250);
P2 = optimvar('P2','LowerBound',3000,'UpperBound',9000);
I1 = optimvar('I1','LowerBound',0,'UpperBound',192000);
I2 = optimvar('I2','LowerBound',0,'UpperBound',244000);
C = optimvar('C','LowerBound',0,'UpperBound',62000);
LE1 = optimvar('LE1','LowerBound',0);
LE2 = optimvar('LE2','LowerBound',0,'UpperBound',142000);
HE1 = optimvar('HE1','LowerBound',0);
HE2 = optimvar('HE2','LowerBound',0);
HPS = optimvar('HPS','LowerBound',0);
MPS = optimvar('MPS','LowerBound',271536);
LPS = optimvar('LPS','LowerBound',100623);
BF1 = optimvar('BF1','LowerBound',0);
BF2 = optimvar('BF2','LowerBound',0);
EP = optimvar('EP','LowerBound',0);
PP = optimvar('PP','LowerBound',0);

% ##### Objective Function #####
linprob = optimproblem('Objective',0.002614*HPS + 0.0239*PP + 0.009825*EP);
% showproblem(linprob)

% ##### Constraints #####
linprob.Constraints.cons1 = I1 - HE1 <= 132000;
linprob.Constraints.cons2 = EP + PP >= 12000;
linprob.Constraints.cons3 = P1 + P2 + PP >= 24550;

linprob.Constraints.econs1 = LE2 + HE2 == I2;
linprob.Constraints.econs2 = LE1 + LE2 + BF2 == LPS;
linprob.Constraints.econs3 = I1 + I2 + BF1 == HPS;
linprob.Constraints.econs4 = C + MPS + LPS == HPS;
linprob.Constraints.econs5 = LE1 + HE1 + C == I1;
linprob.Constraints.econs6 = HE1 + HE2 + BF1 == BF2 + MPS;
linprob.Constraints.econs7 = 1267.8*HE1 + 1251.4*LE1 + 192*C + 3413*P1 == 1359.8*I1;
linprob.Constraints.econs8 = 1267.8*HE2 + 1251.4*LE2 + 3413*P2 == 1359.8*I2;

% ##### solve the problem #####
[linsol, fval] = solve(linprob);

%%
% #################################################################
% Matching nodes: naive trial
% #################################################################
clear
clc
% --- 1D simple case ---
% X = (1:5)'
% Y = (2:2:10)';
% normDiff = X - Y;
% m = size(X,1);
% n = size(Y,1);

% --- 3D simple case ---
X = [0,0,0; 0,1,0; 1,0,0; 1,1,0; 2,2,0];
Y = X;
Y(:,3) = Y(:,3) + 1;
% Y =  [0,0.1,1; 0.2,1,1; 1,1.8,1];
m = size(X,1);
n = size(Y,1);

tmp1 = repmat(X, 1, 1, n);
tmp1 = permute(tmp1, [1,3,2]);
tmp2 = repmat(Y, 1, 1, m);
tmp2 = permute(tmp2, [3,1,2]);
normDiff = tmp1 - tmp2;
normDiff = sum(normDiff .* normDiff, 3);

coeff = optimvar('coeff', m, n, 'LowerBound', 0);

prob = optimproblem('Objective', sum(sum(coeff .* normDiff)));

prob.Constraints.cons1 = sum(coeff,1) == (ones(1,n) * m/n);
prob.Constraints.cons2 = sum(coeff,2) == ones(m,1);
[linsol, fval] = solve(prob);
linsol.coeff


%%
% ##### checks ##### 
% --- check the construction of diff matrix (xi-yj) ---
tmp = zeros(m,n,3);
for i = 1:m
    for j = 1:n
        tmp(i,j,:) = X(i,:) - Y(j,:);
    end
end
sum(sum(sum(tmp1 - tmp2 ~= tmp)));
% --- check the calculation of norm ---
tmp_norm = vecnorm(tmp, 2, 3);
tmp_norm = tmp_norm .* tmp_norm;
sum(normDiff ~= int32(tmp_norm));



