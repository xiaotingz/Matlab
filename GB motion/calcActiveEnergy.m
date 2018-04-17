% -------------------------------------------------------------------
% Notice
%     Somehow the literature plots like log10-scale, but activation energy
%     Q has to be aquired from ln-scale.
% Read in data and modify X_axis
%     Notice the 1000s are just to fill space, need to filter them out 
%     The values are a^2/(t*f(alpha)), with unit = 10^-9 cm^2/sec
% Kprime = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','Q73:S90');
Kprime = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','F73:H90');
x_500 = 1000/(500+273);
x_600 = 1000/(600+273);
x_650 = 1000/(650+273);
x_700 = 1000/(700+273);
x_750 = 1000/(750+273);
x_800 = 1000/(800+273);
X = [x_500*ones(4,1); x_600*ones(3,1); x_650*ones(2,1); x_700*ones(4,1); x_750*ones(3,1); x_800*ones(2,1)];
X_temp = [(273+500)*ones(4,1); (273+600)*ones(3,1); (273+650)*ones(2,1); (273+700)*ones(4,1); (273+750)*ones(3,1); (273+800)*ones(2,1)];


% -------------------------------------------------------------------
% Filter out the 1000s to get real data
seed1_raw = [X, Kprime(:,1), X_temp];
seed2_raw = [X, Kprime(:,2), X_temp];
seed3_raw = [X, Kprime(:,3), X_temp];

bool_seed1 = (seed1_raw(:,2) ~= 1000);
bool_seed2 = (seed2_raw(:,2) ~= 1000);
bool_seed3 = (seed3_raw(:,2) ~= 1000);

seed1_data = seed1_raw(bool_seed1,:);
seed2_data = seed2_raw(bool_seed2,:);
seed3_data = seed3_raw(bool_seed3,:);

% -------------------------------------------------------------------
% Linear regression
    % Y = log(Kprimt*T), unit doesn't matter because it's always canceled in y_diff. (using cm^2*K/sec)
load accidents
x1 = [ones(length(seed1_data),1), seed1_data(:,1)];
y1 = log(seed1_data(:,2).*seed1_data(:,3) * 10^(-9));
x2 = [ones(length(seed2_data),1), seed2_data(:,1)];
y2 = log(seed2_data(:,2).*seed2_data(:,3) * 10^(-9));
x3 = [ones(length(seed3_data),1), seed3_data(:,1)];
y3 = log(seed3_data(:,2).*seed3_data(:,3) * 10^(-9));

format long
b1 = x1\y1;
b2 = x2\y2;
b3 = x3\y3;

% -------------------------------------------------------------------
% Plot
color1 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980];
color3 = [0.9290, 0.6940, 0.1250;];

% figure(1)
scatter(x1(:,2), y1, 60, color1, 'filled');
hold on
yCalc1 = x1*b1;
l1 = plot(x1(:,2), yCalc1, 'LineWidth',2, 'Color', color1);
% figure(2)
scatter(x2(:,2), y2, 60, color2, 'filled');
yCalc2 = x2*b2;
l2 = plot(x2(:,2), yCalc2, 'LineWidth',2, 'Color', color2);
% figure(3)
scatter(x3(:,2), y3, 60, color3, 'filled');
yCalc3 = x3*b3;
l3 = plot(x3(:,2), yCalc3, 'LineWidth',2, 'Color', color3);

ax = gca;
set(ax,'fontsize',21, 'LineWidth', 2, 'box','on');
xlabel('1000/T, 1/K','FontWeight','bold','FontSize',21);
% ylabel('log(Reduced Mobility), 10^-9 cm^2/sec','FontWeight','bold','FontSize',21);
ylabel('log(K''T), cm^2/sec','FontWeight','bold','FontSize',21);
set(gcf, 'Position', [500, 500, 650, 700])

%% -------------------------------------------------------------------
% Calc active energy and diffusion coefficient & note in legend
% K' = K/f(alpha) = a^2/f(alpha)/t = 2*gamma*M = 2*gamma*D0/NkT*exp(-Q/kT)
% Boltzman constant (k) = 8.617×10^?5 eV/K = 1.38*10^-16 erg/K
%  plot ln(K'T)   =>   Q = -slope*k
%                 =>   D0 = intercept * N*k/(2*gamma) || exp(intercept) * N*k/(2*gamma)
%                      Viswanathan_1973 data: N=2*10^15 atom/cm^2; gamma=475 - 600 erg/cm^2

k = 8.617*10^(-5);
k_erg = 1.38*10^(-16);
N = 2*10^15;
eVToKcal_mol = 23.0605419453293;
gamma = 550;

Q1 = - b1(2)*1000 * k;
Q1_eV = roundn(Q1, -2);
Q1_Kcal_mol = roundn(Q1*eVToKcal_mol,-2);
% D0_1 = exp(b1(1)) * N*k_erg/gamma/2;

Q2 = - b2(2)*1000 * k;
Q2_eV = roundn(Q2, -2);
Q2_Kcal_mol = roundn(Q2*eVToKcal_mol,-2);
% D0_2 = exp(b2(1)) * N*k_erg/gamma/2;

Q3 = - b3(2)*1000 * k;
Q3_eV = roundn(Q3, -2);
Q3_Kcal_mol = roundn(Q3*eVToKcal_mol,-2);
% D0_3 = exp(b3(1)) * N*k_erg/gamma/2;


Rsq1 = 1 - sum((y1 - yCalc1).^2)/sum((y1 - mean(y1)).^2);
Rsq2 = 1 - sum((y2 - yCalc2).^2)/sum((y2 - mean(y2)).^2);
Rsq3 = 1 - sum((y3 - yCalc3).^2)/sum((y3 - mean(y3)).^2);

% l1info = ['seed1 slope=', num2str(b1(2)), ', interc=', num2str(b1(1))];
% l2info = ['seed2 slope=', num2str(b2(2)), ', interc=', num2str(b2(1))];
% l3info = ['seed3 slope=', num2str(b3(2)), ', interc=', num2str(b3(1))];
%%
l1info = ['seed1 Q1 = ', num2str(Q1_eV), ' eV = ', num2str(Q1_Kcal_mol), ' kcal/mol'];
l2info = ['seed2 Q2 = ', num2str(Q2_eV), ' eV = ', num2str(Q2_Kcal_mol), ' kcal/mol'];
l3info = ['seed3 Q3 = ', num2str(Q3_eV), ' eV = ', num2str(Q3_Kcal_mol), ' kcal/mol'];
legend([l1, l2, l3],{l1info, l2info, l3info}, 'Location','northoutside')

set(gcf, 'Position', [500, 500, 650, 700])





