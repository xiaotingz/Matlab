%% Electrolyte 1
plot(x071017_1_scan1_C01(:,1), x071017_1_scan1_C01(:,2),'LineWidth',1.5)
hold on
plot(x071017_1_scan2_C01(:,1), x071017_1_scan2_C01(:,2),'LineWidth',1.5)
plot(x071017_1_scan3_C01(:,1), x071017_1_scan3_C01(:,2),'LineWidth',1.5)
plot(x071017_2_scan1_C01(:,1), x071017_2_scan1_C01(:,2),'LineWidth',1.5)
plot(x071017_2_scan2_C01(:,1), x071017_2_scan2_C01(:,2),'LineWidth',1.5)
plot(x071017_2_scan3_C01(:,1), x071017_2_scan3_C01(:,2),'LineWidth',1.5)

plot(x071117_1_scan1_C01(:,1), x071117_1_scan1_C01(:,2),'LineWidth',3)
plot(x071117_1_scan2_C01(:,1), x071117_1_scan2_C01(:,2),'LineWidth',1.5)
plot(x071117_1_scan3_C01(:,1), x071117_1_scan3_C01(:,2),'LineWidth',1.5)

plot(x071217_thickPitted_scan1_C01(:,1), x071217_thickPitted_scan1_C01(:,2),'--','LineWidth',3)
plot(x071217_thickPitted_scan2_C01(:,1), x071217_thickPitted_scan2_C01(:,2),'--','LineWidth',1.5)
plot(x071217_thickPitted_scan3_C01(:,1), x071217_thickPitted_scan3_C01(:,2),'--','LineWidth',1.5)
plot(x071217_thickPitted_scan4_C01(:,1), x071217_thickPitted_scan4_C01(:,2),'--','LineWidth',1.5)
% plot(x071017_2_scan4_C01(:,1), x071017_2_scan4_C01(:,2),':','LineWidth',1.5)
% plot(x071017_2_scan5_C01(:,1), x071017_2_scan5_C01(:,2),':','LineWidth',1.5)

legend('0710\_1\_scan1', '0710\_1\_scan2', '0710\_1\_scan3', ...
    '0710\_2\_scan1', '0710\_2\_scan2', '0710\_2\_scan3', ...
    '0711\_1\_scan1', '0711\_1\_scan2', '0711\_1\_scan3', ...
    '0712\_thickPit\_scan1', '0712\_thickPit\_scan2', '0712\_thickPit\_scan3', '0712\_thickPit\_scan4', ...
    'Location','northeastoutside' )
grid on
grid minor
set(gca,'fontsize',15)
xlabel('Voltage, V','FontSize',18);
ylabel('Current, mA','FontSize',18);
xlim([1,2.3])
ylim([0,25])

%% Electrolyte 2

plot(x071117_2thick_scan1_C01(:,1), x071117_2thick_scan1_C01(:,2),'--','LineWidth',1.5)
hold on
plot(x071117_2thick_scan2_C01(:,1), x071117_2thick_scan2_C01(:,2),'--','LineWidth',1.5)

legend('0710\_2thick\_scan1', '0710\_2thick\_scan2')
grid on
grid minor
set(gca,'fontsize',15)
xlabel('Voltage, V','FontSize',18);
ylabel('Current, mA','FontSize',18);
% xlim([1,2])
% ylim([0,25])

%% Polish data, Electrolyte1 thick sample
plot(x071117_2thick_polish_C01(:,1), x071117_2thick_polish_C01(:,2),'--','LineWidth',1.5)
hold on
plot(x071117_2thick_polish2_C01(:,1), x071117_2thick_polish2_C01(:,2),'--','LineWidth',1.5)
plot(x071117_2thick_polish3_C01(:,1), x071117_2thick_polish3_C01(:,2),'--','LineWidth',1.5)
plot(x071117_2thick_polish4_C01(:,1), x071117_2thick_polish4_C01(:,2),'--','LineWidth',1.5)

plot(x071217_thickPitted_polish1_C01(:,1), x071217_thickPitted_polish1_C01(:,2),'--','LineWidth',1.5)
plot(x071217_thickPitted_polish2_C01(:,1), x071217_thickPitted_polish2_C01(:,2),'--','LineWidth',1.5)
plot(x071217_thickPitted_polish3_C01(:,1), x071217_thickPitted_polish3_C01(:,2),'--','LineWidth',1.5)
plot(x071217_thickPitted_polish4_C01(:,1), x071217_thickPitted_polish4_C01(:,2),'--','LineWidth',1.5)

plot(x071217_1thick_Polish_C01(:,1), x071217_1thick_Polish_C01(:,2),'--','LineWidth',1.5)

legend('0711\_2thick\_polish1, 1.84V, GAS', '0711\_2thick\_polish2, 1.85V, GAS', '0711\_2thick\_polish3, 1.68V, GAS', '0711\_2thick\_polish4, 1.68V, GAS', ...
        '0712\_tP\_polish1, 2.0V, GAS', '0712\_tP\_polish2, 1.9V, GAS', '0712\_tP\_polish3, 1.8V, GAS', '0712\_tP\_polish4, 1.9V, GAS', ...
        '0712\_1thick\_polish, 1.5V', 'Location','northeastoutside' )
grid on
grid minor
set(gca,'fontsize',15)
xlabel('Time, sec','FontSize',18);
ylabel('Current, mA','FontSize',18);
% xlim([50,400]) 
%% Polish data

plot(x071917_6round_polish1_02_CA_C01(:,1), x071917_6round_polish1_02_CA_C01(:,2),'LineWidth',1.5)
hold on
plot(x071917_6round_polish2_02_CA_C01(:,1), x071917_6round_polish2_02_CA_C01(:,2),'LineWidth',1.5)

legend('0719\_6round\_polish1, 1.6V', '0719\_6round\_polish2, 1.6V')

grid on
grid minor
set(gca,'fontsize',15)
xlabel('Time, second','FontSize',18);
ylabel('Current, mA','FontSize',18);
% xlim([100,600])
% ylim([5,20])


%% Polish data, 071617
% GP tweezers: --
% thick line just to differ
% Gas sample 
plot(x071617_1round_polish1_C01(:,1), x071617_1round_polish1_C01(:,2),'--','LineWidth',1.5)
hold on
plot(x071617_2roll_polish1_C01(:,1), x071617_2roll_polish1_C01(:,2),'--','LineWidth',1.5)
plot(x071617_1round_polish2_C01(:,1), x071617_1round_polish2_C01(:,2),'LineWidth',1.5)

plot(x071617_3round_polish1_C01(:,1), x071617_3round_polish1_C01(:,2),'LineWidth',1.5)
plot(x071617_3round_polish2_C01(:,1), x071617_3round_polish2_C01(:,2),'LineWidth',1.5)
plot(x071617_3round_polish3_C01(:,1), x071617_3round_polish3_C01(:,2),'--','LineWidth',1.5)

plot(x071617_4round_polish1_C01(:,1), x071617_4round_polish1_C01(:,2),'LineWidth',1.5)
plot(x071617_4round_polish2_C01(:,1), x071617_4round_polish2_C01(:,2),'LineWidth',3)

plot(x071617_5roll_polish1_C01(:,1), x071617_5roll_polish1_C01(:,2),'LineWidth',3)

plot(x071617_1round_polish3_C01(:,1), x071617_1round_polish3_C01(:,2),'LineWidth',3)
plot(x071617_1round_polish4_C01(:,1), x071617_1round_polish4_C01(:,2),'LineWidth',3)
plot(x071617_1round_polish5_C01(:,1), x071617_1round_polish5_C01(:,2),'LineWidth',3)

plot(x071617_6round_polish1_02_CA_C01(:,1), x071617_6round_polish1_02_CA_C01(:,2),'LineWidth',3)

 legend('1round\_polish1, 1.7V', '2roll\_polish1, 1.7V', '1round\_polish2, 1.9V, gas', '3round\_polish1, 1.6V', ...
     '3round\_polish2, 1.5V', '3round\_polish3, 1.8V, gas', '4round\_polish1, 1.6V, newEl', '4round\_polish2, 1.7V', ...
     '5roll\_polish1, 1.7V', '1round\_polish3, 2.3V, gas', '1round\_polish4, 2.2V, gas', '1round\_polish5, 1.8V, gas',...
     '6round\_polish1, 1V\_1.6V');
 text(150,120,'GP tweezers: --; CSS tweezers: solid line; thick line no meaning','FontSize',14)
%  text(100,110,'CSS tweezers: solid line','FontSize',14)
grid on
grid minor
set(gca,'fontsize',15)
xlabel('Voltage, V','FontSize',18);
ylabel('Current, mA','FontSize',18);

% axis([100,500,6,20])
% xlim([1,2])

%% Scan, 071617
% GP tweezers: --
plot(x071617_1round_scan1_C01(:,1), x071617_1round_scan1_C01(:,2),'--','LineWidth',1.5)
hold on
plot(x071617_1round_scan2_C01(:,1), x071617_1round_scan2_C01(:,2),'LineWidth',1.5)
plot(x071617_1round_scan3_C01(:,1), x071617_1round_scan3_C01(:,2),'LineWidth',1.5)

plot(x071617_2roll_scan1_C01(:,1), x071617_2roll_scan1_C01(:,2),'LineWidth',1.5)
plot(x071617_2roll_scan2_C01(:,1), x071617_2roll_scan2_C01(:,2),'LineWidth',1.5)
plot(x071617_3round_scan1_C01(:,1), x071617_3round_scan1_C01(:,2),'--','LineWidth',1.5)


 legend('0716\_1round\_scan1', '0716\_1round\_scan2', '0716\_1round\_scan3', ...
    '0716\_2roll\_scan1', '0716\_2roll\_scan2', '0716\_3round\_scan1','Location','northeastoutside')
 
 text(0.5,32,'GP tweezers: --; CSS tweezers: solid line','FontSize',14)
%  text(100,110,'CSS tweezers: solid line','FontSize',14)
grid on
grid minor
set(gca,'fontsize',15)
xlabel('Voltage, V','FontSize',18);
ylabel('Current, mA','FontSize',18);
axis([0,2.5,0,35])

%% Scan, 071917
plot(x071817_1round_scan1_C01(:,1), x071817_1round_scan1_C01(:,2),'LineWidth',1.5)
hold on
plot(x071817_1round_scan2_C01(:,1), x071817_1round_scan2_C01(:,2),'LineWidth',1.5)

legend('0718\_1round\_scan1', '0718\_1round\_scan2','Location','northwest')
 
 text(1.2,16.5,'Electrolyte3','FontSize',14)
%  text(100,110,'CSS tweezers: solid line','FontSize',14)
grid on
grid minor
set(gca,'fontsize',15)
xlabel('Voltage, V','FontSize',18);
ylabel('Current, mA','FontSize',18);
axis([0,2.5,0,18])


%% 071917 electrolyte1
plot(x071917_3round_polish1(:,1), x071917_3round_polish1(:,2),'LineWidth',3)
hold on
plot(x071917_3round_polish2(:,1), x071917_3round_polish2(:,2),'LineWidth',1.5)
plot(x071917_4round_polish1(:,1), x071917_4round_polish1(:,2),'LineWidth',1.5)
plot(x071917_4round_polish2(:,1), x071917_4round_polish2(:,2),'LineWidth',1.5)
plot(x071917_5roll_polish1(:,1), x071917_5roll_polish1(:,2),'--','LineWidth',3)
plot(x071917_5roll_polish2(:,1), x071917_5roll_polish2(:,2),'--','LineWidth',1.5)
plot(x071917_5roll_polish3(:,1), x071917_5roll_polish3(:,2),'--','LineWidth',1.5)
plot(x071917_5roll_polish4(:,1), x071917_5roll_polish4(:,2),'--','LineWidth',1.5)
plot(x071917_5roll_polish5(:,1), x071917_5roll_polish5(:,2),'--','LineWidth',1.5)
plot(x071917_6round_polish1(:,1), x071917_6round_polish1(:,2),'LineWidth',3)
plot(x071917_6round_polish2(:,1), x071917_6round_polish2(:,2),'-.','LineWidth',1.5)

 legend('0719\_3round\_polish1, 1.6V', '0719\_3round\_polish2, 1.6V', '0719\_4round\_polish1, 1.5V', '0719\_4round\_polish2, 1.5V', ...
     '0719\_5roll\_polish1, 1.5V', '0719\_5roll\_polish2, 1.5V', '0719\_5roll\_polish3, 1.7V', '0719\_5roll\_polish4, 1.7V', '0719\_5roll\_polish5, 1.5V', ...
     '0719\_6round\_polish1, 1.6V', '0719\_6round\_polish2, 1.6V','Location','northeastoutside')
 
%  text(200,15,'Electrolyte3','FontSize',14)
 text(250,4,'Temp: 18-21C; bold line: fresh electrolyte1; -. line: no meaning','FontSize',14)
grid on
grid minor
set(gca,'fontsize',15)
xlabel('Voltage, V','FontSize',18);
ylabel('Current, mA','FontSize',18);
axis([200,700,3,10])


%%
plot(x071917_1round_polish1_C01(:,1), x071917_1round_polish1_C01(:,2),'-.','LineWidth',1.5)

%  legend('0719\_3round\_polish1, 1.6V', '0719\_3round\_polish2, 1.6V', '0719\_4round\_polish1, 1.5V', '0719\_4round\_polish2, 1.5V', ...
%      '0719\_5roll\_polish1, 1.5V', '0719\_5roll\_polish2, 1.5V', '0719\_5roll\_polish3, 1.7V', '0719\_5roll\_polish4, 1.7V', '0719\_5roll\_polish5, 1.5V', ...
%      '0719\_6round\_polish1, 1.6V', '0719\_6round\_polish2, 1.6V','Location','northeastoutside')
%  
grid on
grid minor
set(gca,'fontsize',15)
xlabel('Voltage, V','FontSize',18);
ylabel('Current, mA','FontSize',18);
axis([200,700,3,10])