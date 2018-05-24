% ------------------- second derivative of twist[100] -------------------
% -- parameter, using Cu
[pars,AlCuparameter] = makeparvec();
a = pars(10);   % 100 twist maximum energy
c = pars(11); 
diff = deg2rad(0.1);
lowLim = deg2rad(0);
highLim = deg2rad(45);
% -- analytic
theta = (lowLim:diff:highLim).';
sins = sin(theta);
coss = cos(theta);
firstDev = a.*(coss - a.*coss.*log(sins) - a.*coss);
secondDev = a.*(- sins + c.*sins.*log(sins) - c.*coss.*1./sins.*coss + c.*sins);
figure(1)
plot(rad2deg(theta), firstDev,'LineWidth',2)
figure(2)
plot(rad2deg(theta), secondDev,'LineWidth',2)
% -- finite difference
sins = sin(theta) ;
xlogx = sins.*log(sins);
% xlogx(isnan(xlogx))=0; % Force the limit to zero as x -> 0.
en =  a.*(sins - c.*xlogx);
en_pos = en(3:end);
en_neg = en(1:end-2);
en_point = en(2:end-1);
firstDiff = (en_pos - en_point)./diff;
figure(3)
plot(rad2deg(theta(2:end-1)), firstDiff,'LineWidth',2)
hold on
secondDiff = (en_pos - 2.*en_point + en_neg)./(diff^2);
plot(rad2deg(theta(2:end-1)), secondDiff,'LineWidth',2)
legend('secondDev', 'secondDiff');


%% ------------------- second derivative of 60@[111] -------------------
% -- parameter, using Cu
[pars,AlCuparameter] = makeparvec();
diff = deg2rad(0.1);
lowLim = deg2rad(0);
highLim = deg2rad(60);


eta = [lowLim : diff : highLim].';
phi = eta;
en = zeros(length(eta), length(phi));

for i = 1:length(eta)
    for j = 1:length(phi)
        geom111(1) = pi/3;
        geom111(2) = eta(i);
        geom111(3) = phi(j);
        en(i,j) = set111_oneGB(geom111,pars);
    end
end

% -- assume i=eta, j=phi
en_iPos = en(3:end,:);
en_i = en(2:end-1,:);
en_iNeg = en(1:end-2,:);
enSecondDiff_11 = (en_iPos - 2.*en_i + en_iNeg)./(diff^2);

en_jPos = en(:,3:end);
en_j = en(:,2:end-1);
en_jNeg = en(:,1:end-2);
enSecondDiff_22 = (en_jPos - 2.*en_j + en_jNeg)./(diff^2);

en_ijPos = en(3:end,3:end);
en_iPosjNeg = en(3:end,1:end-2);
en_iNegjPos = en(1:end-2,3:end);
en_ijNeg = en(1:end-2,1:end-2);
enSecondDiff_12 = (en_ijPos - en_iPosjNeg - en_iNegjPos + en_ijNeg)/(4*diff^2);

eigenValues = zeros(length(enSecondDiff_12), 2);
for i = 1:length(en)-2
    for j = 1:length(en)-2
        mat = [enSecondDiff_11(i,i), enSecondDiff_12(i,j); enSecondDiff_12(j,i), enSecondDiff_22(j,j)];
        eigenValue= eig(mat);
        eigenValues(i,j,1) = eigenValue(1);
        eigenValues(i,j,2) = eigenValue(2);
    end
end

% -- histograms
% figure(1)
% histogram(enSecondDiff_11,100)
% set(gca, 'YScale', 'log')
% xlabel('$\frac{1}{{\partial}^2 \eta}$', 'Interpreter','latex')
% figure(2)
% histogram(enSecondDiff_22,100)
% set(gca, 'YScale', 'log')
% xlabel('$\frac{1}{{\partial}^2 \phi}$', 'Interpreter','latex') 
% figure(3)
% histogram(enSecondDiff_12,100)
% xlabel('$\frac{1}{\partial \eta \partial \phi}$', 'Interpreter','latex')
% figure(9)
% histogram(eigenValues(:,:,1),200)
% set(gca, 'YScale', 'log')
% hold on
% histogram(eigenValues(:,:,2),200)
% xlabel('eigenValues')
% legend('eigenValues 1','eigenValues 2','Location','northwest')


% -- color images
figure(4)
imagesc(rad2deg(phi(2:end-1)), rad2deg(eta(2:end-1)),enSecondDiff_11,[-0.1,0])
title('$\frac{1}{{\partial}^2 \eta}$ fine', 'Interpreter','latex')
colorbar
xlabel('Tilt Angle $\phi$, $^{\circ}$', 'Interpreter','latex')
ylabel('Tilt Assymetry Angle $\eta$, $^{\circ}$', 'Interpreter','latex')
figure(5)
imagesc(rad2deg(phi(2:end-1)), rad2deg(eta(2:end-1)), enSecondDiff_22)
title('$\frac{1}{{\partial}^2 \phi}$ fine', 'Interpreter','latex')
colorbar
xlabel('Tilt Angle $\phi$, $^{\circ}$', 'Interpreter','latex')
ylabel('Tilt Assymetry Angle $\eta$, $^{\circ}$', 'Interpreter','latex')
figure(6)
imagesc(rad2deg(phi(2:end-1)), rad2deg(eta(2:end-1)), enSecondDiff_12)
title('$\frac{1}{\partial \eta \partial \phi}$', 'Interpreter','latex')
colorbar
xlabel('Tilt Angle $\phi$, $^{\circ}$', 'Interpreter','latex')
ylabel('Tilt Assymetry Angle $\eta$, $^{\circ}$', 'Interpreter','latex')
figure(7)
imagesc(rad2deg(phi(2:end-1)), rad2deg(eta(2:end-1)), eigenValues(:,:,1))
title('eigenValue 1')
colorbar
xlabel('Tilt Angle $\phi$, $^{\circ}$', 'Interpreter','latex')
ylabel('Tilt Assymetry Angle $\eta$, $^{\circ}$', 'Interpreter','latex')
figure(8)
imagesc(rad2deg(phi(2:end-1)), rad2deg(eta(2:end-1)), eigenValues(:,:,2))
title('eigenValue 2')
colorbar
xlabel('Tilt Angle $\phi$, $^{\circ}$', 'Interpreter','latex')
ylabel('Tilt Assymetry Angle $\eta$, $^{\circ}$', 'Interpreter','latex')


%% ------------------- second derivative of 38.9@[110] -------------------
[pars,AlCuparameter] = makeparvec();
diff = deg2rad(0.1);
lowLim = deg2rad(0);
highLim = deg2rad(180);


eta = [lowLim : diff : highLim].';
phi = eta;
en = zeros(length(eta), length(phi));

for i = 1:length(eta)
    for j = 1:length(phi)
        geom110(1) = deg2rad(38.9);
        geom110(2) = eta(i);
        geom110(3) = phi(j);
        en(i,j) = set110_oneGB(geom110,pars);
    end
end

% -- assume i=eta, j=phi
en_iPos = en(3:end,:);
en_i = en(2:end-1,:);
en_iNeg = en(1:end-2,:);
enSecondDiff_11 = (en_iPos - 2.*en_i + en_iNeg)./(diff^2);

en_jPos = en(:,3:end);
en_j = en(:,2:end-1);
en_jNeg = en(:,1:end-2);
enSecondDiff_22 = (en_jPos - 2.*en_j + en_jNeg)./(diff^2);

en_ijPos = en(3:end,3:end);
en_iPosjNeg = en(3:end,1:end-2);
en_iNegjPos = en(1:end-2,3:end);
en_ijNeg = en(1:end-2,1:end-2);
enSecondDiff_12 = (en_ijPos - en_iPosjNeg - en_iNegjPos + en_ijNeg)/(4*diff^2);

eigenValues = zeros(length(enSecondDiff_12), 2);
for i = 1:length(en)-2
    for j = 1:length(en)-2
        mat = [enSecondDiff_11(i,i), enSecondDiff_12(i,j); enSecondDiff_12(j,i), enSecondDiff_22(j,j)];
        eigenValue= eig(mat);
        eigenValues(i,j,1) = eigenValue(1);
        eigenValues(i,j,2) = eigenValue(2);
    end
end

% -- histograms
figure(1)
histogram(enSecondDiff_11,100)
set(gca, 'YScale', 'log')
xlabel('$\frac{1}{{\partial}^2 \eta}$', 'Interpreter','latex')
figure(2)
histogram(enSecondDiff_22,100)
set(gca, 'YScale', 'log')
xlabel('$\frac{1}{{\partial}^2 \phi}$', 'Interpreter','latex') 
figure(3)
histogram(enSecondDiff_12,100)
xlabel('$\frac{1}{\partial \eta \partial \phi}$', 'Interpreter','latex')
figure(9)
histogram(eigenValues(:,:,1),200)
set(gca, 'YScale', 'log')
hold on
histogram(eigenValues(:,:,2),200)
xlabel('eigenValues')
legend('eigenValues 1','eigenValues 2','Location','northwest')


% -- color images
figure(4)
imagesc(rad2deg(phi(2:end-1)), rad2deg(eta(2:end-1)),enSecondDiff_11)
title('$\frac{1}{{\partial}^2 \eta}$ fine', 'Interpreter','latex')
colorbar
xlabel('Tilt Angle $\phi$, $^{\circ}$', 'Interpreter','latex')
ylabel('Tilt Assymetry Angle $\eta$, $^{\circ}$', 'Interpreter','latex')
figure(5)
imagesc(rad2deg(phi(2:end-1)), rad2deg(eta(2:end-1)), enSecondDiff_22)
title('$\frac{1}{{\partial}^2 \phi}$ fine', 'Interpreter','latex')
colorbar
xlabel('Tilt Angle $\phi$, $^{\circ}$', 'Interpreter','latex')
ylabel('Tilt Assymetry Angle $\eta$, $^{\circ}$', 'Interpreter','latex')
figure(6)
imagesc(rad2deg(phi(2:end-1)), rad2deg(eta(2:end-1)), enSecondDiff_12)
title('$\frac{1}{\partial \eta \partial \phi}$', 'Interpreter','latex')
colorbar
xlabel('Tilt Angle $\phi$, $^{\circ}$', 'Interpreter','latex')
ylabel('Tilt Assymetry Angle $\eta$, $^{\circ}$', 'Interpreter','latex')
figure(7)
imagesc(rad2deg(phi(2:end-1)), rad2deg(eta(2:end-1)), eigenValues(:,:,1))
title('eigenValue 1')
colorbar
xlabel('Tilt Angle $\phi$, $^{\circ}$', 'Interpreter','latex')
ylabel('Tilt Assymetry Angle $\eta$, $^{\circ}$', 'Interpreter','latex')
figure(8)
imagesc(rad2deg(phi(2:end-1)), rad2deg(eta(2:end-1)), eigenValues(:,:,2))
title('eigenValue 2')
colorbar
xlabel('Tilt Angle $\phi$, $^{\circ}$', 'Interpreter','latex')
ylabel('Tilt Assymetry Angle $\eta$, $^{\circ}$', 'Interpreter','latex')


