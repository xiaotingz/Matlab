% MGB = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','Q73:S90');
MGB = xlsread('/Users/xiaotingzhong/Dropbox/Paper & Documentation/Documentations/Cu_moved GB.xlsx','F73:H90');

x_500 = 1000/(500+273);
x_600 = 1000/(600+273);
x_650 = 1000/(650+273);
x_700 = 1000/(700+273);
x_750 = 1000/(750+273);
x_800 = 1000/(800+273);
X = [x_500*ones(4,1); x_600*ones(3,1); x_650*ones(2,1); x_700*ones(4,1); x_750*ones(3,1); x_800*ones(2,1)];
% X_tick = [x_500, x_600, x_650, x_700];

% seed1
scatter(X, MGB(:,1), 80, 'filled');
hold on
scatter(X, MGB(:,2), 80, 'filled');
scatter(X, MGB(:,3), 80, 'filled');

% seed2


% seed3

% axis and labels
ax = gca;
set(ax, 'YScale', 'log')
set(ax,'fontsize',21, 'LineWidth', 2, 'box','on')
xlabel('1000/T, 1/K','FontWeight','bold','FontSize',21);
ylabel('Reduced Mobility, 10^-9 cm^2/sec','FontWeight','bold','FontSize',21);
legend('seed1=23°@[113]', 'seed2=39°@[110]', 'seed3=52°@[123]');
axis([0.9, 1.3, 0, 200])

