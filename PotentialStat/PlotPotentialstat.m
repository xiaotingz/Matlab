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
%% Polish data, Electrolyte1 thin sample

plot(x071217_2thin_Polish1_C01(:,1), x071217_2thin_Polish1_C01(:,2),'LineWidth',1.5)
hold on
plot(x071217_2thin_Polish2_C01(:,1), x071217_2thin_Polish2_C01(:,2),'LineWidth',1.5)
plot(x071217_3thin_Polish1_C01(:,1), x071217_3thin_Polish1_C01(:,2),'LineWidth',1.5)

legend('0712\_2thin\_polish1, 1.50V, IC', '0712\_2thin\_polish2, 1.50V, IC','0712\_3thin\_polish1, 1.50V')

grid on
grid minor
set(gca,'fontsize',15)
xlabel('Voltage, V','FontSize',18);
ylabel('Current, mA','FontSize',18);
xlim([100,600])
ylim([5,20])


%% Polish data, Electrolyte2
% Fresh electrolyte in bold
% thick sample in dash
% Gas sample 
plot(x071117_2thick_polish5_C01(:,1), x071117_2thick_polish5_C01(:,2),'--','LineWidth',3)
hold on
plot(x071117_2thick_polish6_C01(:,1), x071117_2thick_polish6_C01(:,2),'--','LineWidth',1.5)
plot(x071117_2thick_polish7_C01(:,1), x071117_2thick_polish7_C01(:,2),'--','LineWidth',1.5)
plot(x071117_2thick_polish8_C01(:,1), x071117_2thick_polish8_C01(:,2),'--','LineWidth',1.5)

plot(x071117_3thin_polish1_C01(:,1), x071117_3thin_polish1_C01(:,2),'LineWidth',1.5)
plot(x071117_3thin_polish2_C01(:,1), x071117_3thin_polish2_C01(:,2),'LineWidth',1.5)
plot(x071117_4thin_polish1_C01(:,1), x071117_4thin_polish1_C01(:,2),'LineWidth',1.5)

legend('0711\_2thick\_polish5, 1.55V', '0711\_2thick\_polish6, 1.82V', '0711\_2thick\_polish7 1.93V, GAS', ...
    '0711\_2thick\_polish8 1.89V, GAS', '0711\_3thin\_polish1, 1.68V', '0711\_3thin\_polish2, 1.5V', '0711\_4thin\_polish1, 1.5V')
grid on
grid minor
set(gca,'fontsize',15)
xlabel('Voltage, V','FontSize',18);
ylabel('Current, mA','FontSize',18);
% xlim([1,2])


