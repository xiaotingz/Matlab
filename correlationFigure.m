file1 = textread('A_curvEn_0.05.txt');
A_curvEn_0_05 = file1(:,1:2);
file2 = textread('A_enCurv_0.05.txt');
A_enCurv_0_05 = file2(:,1:2);
file3 = textread('A_gbCurv_10.txt');
A_gbCurv_10 = file3(:,1:2);
file4 = textread('F_curvEn_0.05.txt');
F_curvEn_0_05 = file4(:,1:2);
file5 = textread('F_enCurv_0.05.txt');
F_enCurv_0_05 = file5(:,1:2);
file6 = textread('F_gbCurv_0.02.txt');
F_gbCurv_0_02 = file6(:,1:2);
file7 = textread('F_gbCurv_0.04.txt');
F_gbCurv_0_04 = file7(:,1:2);

figure(1)
scatter(A_curvEn_0_05(:,1),A_curvEn_0_05(:,2),'filled');
title('TWIP Curvature v.s. Energy','FontSize',14);
xlabel('curvature','FontSize',14)
ylabel('energy','FontSize',14)

figure(2)
scatter(A_enCurv_0_05(:,1),A_enCurv_0_05(:,2),'filled');
title('A_enCurv_0_05','FontSize',14);
xlabel('energy','FontSize',14)
ylabel('curvature','FontSize',14)

figure(3)
scatter(A_gbCurv_10(:,1),A_gbCurv_10(:,2),'filled');
title('A_gbCurv_10','FontSize',14);
xlabel('population','FontSize',14)
ylabel('curvature','FontSize',14)

figure(4)
scatter(F_curvEn_0_05(:,1),F_curvEn_0_05(:,2),'filled');
title('F_curvEn_0_05','FontSize',14);
xlabel('curvature','FontSize',14)
ylabel('energy','FontSize',14)

figure(5)
scatter(F_enCurv_0_05(:,1),F_enCurv_0_05(:,2),'filled');
title('F_enCurv_0_05','FontSize',14);
xlabel('energy','FontSize',14)
ylabel('curvature','FontSize',14)

figure(6)
scatter(F_gbCurv_0_02(:,1),F_gbCurv_0_02(:,2),'filled');
title('F_gbCurv_0_02','FontSize',14);
xlabel('population','FontSize',14)
ylabel('curvature','FontSize',14)

figure(7)
scatter(F_gbCurv_0_04(:,1),F_gbCurv_0_04(:,2),'filled');
title('F_gbCurv_0_04','FontSize',14);
xlabel('poplution','FontSize',14)
ylabel('curvature','FontSize',14)

%%
% file = textread('Jan31_A_Curv10_v_energy_TWIP.txt');
file = textread('F_CurvDistri10Ave_v_energy_Ferrite.txt');

cnt = 1;
figure(2)
for i = 1:length(file)
    if file(i,3) >= 1000
        scatter(file(i,1),file(i,2),100,'filled','k');hold on
        line([file(i,1),file(i,1)],[(file(i,2)-file(i,5)),(file(i,2)+file(i,5))],'color','k');
        line([file(i,1)-0.003,file(i,1)+0.003],[(file(i,2)-file(i,5)),(file(i,2)-file(i,5))],'color','k')
        line([file(i,1)-0.003,file(i,1)+0.003],[(file(i,2)+file(i,5)),(file(i,2)+file(i,5))],'color','k')
        F_energyPlot(cnt,1) = file(i,2);
        cnt = cnt + 1;
    end
end
set(gca,'fontsize',18)

% AXIS LABEL
xlabel('Curvature Range, \mum^{-1}','FontSize',20);
ylabel('Average Energy, a.u.','FontSize',20);
% xlabel('Population, MRD','FontSize',17);



% ax = gca;
% Ferrite, no Errorbar
% ax.XTick = [0.5:0.1:0.8];
% axis([0.47,0.8,0.75,1.05]); 
% Austenite, no Errorbar
% ax.YTick = [0.59:0.02:0.69];
% axis([0.45,0.67,0.59,0.69]);
% box on

% Austenite, Errorbar
% axis([0.45,0.67,0.46,0.82]);
% Ferrite, Errorbar
axis([0.47,0.79,0.57,1.15]);


% disable tick on top and right of the box
a = gca;
    % set box property to off and remove background color
set(a,'box','off','color','none')
    % create new, empty axes with box but without ticks
b = axes('Position',get(a,'Position'),'box','on','xtick',[],'ytick',[]);
    % set original axes as active
axes(a)
    % link axes in case of zooming
linkaxes([a b])
