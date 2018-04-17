% Viswanathan data, 9deg, w-normal
Kprime = [22.3, 70].';
T = [875, 950].';
Kprime = Kprime * 10^(-9);
KprimeT = Kprime .* T;
T = 1./T;

% load accidents
x = [ones(length(T),1), T];
y = log(KprimeT);

% format long
b = x\y;

plot(x(:,2), y)

k = 8.617*10^(-5);
eVToKcal_mol = 23.0605419453293;
Q = - (b(2)*k)*eVToKcal_mol

