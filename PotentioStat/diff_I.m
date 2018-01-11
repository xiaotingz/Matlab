
x071917_3round_polish1_I = x071917_3round_polish1(:,2);
x071917_3round_polish2_I = x071917_3round_polish2(:,2);
x071917_4round_polish1_I = x071917_4round_polish1(:,2);
x071917_4round_polish2_I = x071917_4round_polish2(:,2);
x071917_5roll_polish1_I = x071917_5roll_polish1(:,2);
x071917_5roll_polish2_I = x071917_5roll_polish2(:,2);
x071917_5roll_polish3_I = x071917_5roll_polish3(:,2);
x071917_5roll_polish4_I = x071917_5roll_polish4(:,2);
x071917_5roll_polish5_I = x071917_5roll_polish5(:,2);
x071917_6round_polish1_I = x071917_6round_polish1(:,2);
x071917_6round_polish2_I = x071917_6round_polish2(:,2);

x071917_3round_polish1_I_diff = x071917_3round_polish1_I(2:end) - x071917_3round_polish1_I(1:end-1);
x071917_3round_polish2_I_diff = x071917_3round_polish2_I(2:end) - x071917_3round_polish2_I(1:end-1);
x071917_4round_polish1_I_diff = x071917_4round_polish1_I(2:end) - x071917_4round_polish1_I(1:end-1);
x071917_4round_polish2_I_diff = x071917_4round_polish2_I(2:end) - x071917_4round_polish2_I(1:end-1);
x071917_5roll_polish1_I_diff = x071917_5roll_polish1_I(2:end) - x071917_5roll_polish1_I(1:end-1);
x071917_5roll_polish2_I_diff = x071917_5roll_polish2_I(2:end) - x071917_5roll_polish2_I(1:end-1);
x071917_5roll_polish3_I_diff = x071917_5roll_polish3_I(2:end) - x071917_5roll_polish3_I(1:end-1);
x071917_5roll_polish4_I_diff = x071917_5roll_polish4_I(2:end) - x071917_5roll_polish4_I(1:end-1);
x071917_5roll_polish5_I_diff = x071917_5roll_polish5_I(2:end) - x071917_5roll_polish5_I(1:end-1);
x071917_6round_polish1_I_diff = x071917_6round_polish1_I(2:end) - x071917_6round_polish1_I(1:end-1);
x071917_6round_polish2_I_diff = x071917_6round_polish2_I(2:end) - x071917_6round_polish2_I(1:end-1);

x071917_3round_polish1_T = x071917_3round_polish1(:,1);
x071917_3round_polish2_T = x071917_3round_polish2(:,1);
x071917_4round_polish1_T = x071917_4round_polish1(:,1);
x071917_4round_polish2_T = x071917_4round_polish2(:,1);
x071917_5roll_polish1_T = x071917_5roll_polish1(:,1);
x071917_5roll_polish2_T = x071917_5roll_polish2(:,1);
x071917_5roll_polish3_T = x071917_5roll_polish3(:,1);
x071917_5roll_polish4_T = x071917_5roll_polish4(:,1);
x071917_5roll_polish5_T = x071917_5roll_polish5(:,1);
x071917_6round_polish1_T = x071917_6round_polish1(:,1);
x071917_6round_polish2_T = x071917_6round_polish2(:,1);

x071917_3round_polish1_T_diff = x071917_3round_polish1_T(2:end);
x071917_3round_polish2_T_diff = x071917_3round_polish2_T(2:end);
x071917_4round_polish1_T_diff = x071917_4round_polish1_T(2:end);
x071917_4round_polish2_T_diff = x071917_4round_polish2_T(2:end);
x071917_5roll_polish1_T_diff = x071917_5roll_polish1_T(2:end);
x071917_5roll_polish2_T_diff = x071917_5roll_polish2_T(2:end);
x071917_5roll_polish3_T_diff = x071917_5roll_polish3_T(2:end);
x071917_5roll_polish4_T_diff = x071917_5roll_polish4_T(2:end);
x071917_5roll_polish5_T_diff = x071917_5roll_polish5_T(2:end);
x071917_6round_polish1_T_diff = x071917_6round_polish1_T(2:end);
x071917_6round_polish2_T_diff = x071917_6round_polish2_T(2:end);


plot(x071917_3round_polish1_T_diff, x071917_3round_polish1_I_diff,'LineWidth',3)
hold on
plot(x071917_3round_polish2_T_diff, x071917_3round_polish2_I_diff,'LineWidth',1.5)
plot(x071917_4round_polish1_T_diff, x071917_4round_polish1_I_diff,'LineWidth',1.5)
plot(x071917_4round_polish2_T_diff, x071917_4round_polish2_I_diff,'LineWidth',1.5)
plot(x071917_5roll_polish1_T_diff, x071917_5roll_polish1_I_diff,'--','LineWidth',3)
plot(x071917_5roll_polish2_T_diff, x071917_5roll_polish2_I_diff,'--','LineWidth',1.5)
plot(x071917_5roll_polish3_T_diff, x071917_5roll_polish3_I_diff,'--','LineWidth',1.5)
plot(x071917_5roll_polish4_T_diff, x071917_5roll_polish4_I_diff,'--','LineWidth',1.5)
plot(x071917_5roll_polish5_T_diff, x071917_5roll_polish5_I_diff,'--','LineWidth',1.5)
plot(x071917_6round_polish1_T_diff, x071917_6round_polish1_I_diff,'LineWidth',3)
plot(x071917_6round_polish2_T_diff, x071917_6round_polish2_I_diff,'-.','LineWidth',1.5)

 legend('0719\_3round\_polish1, 1.6V', '0719\_3round\_polish2, 1.6V', '0719\_4round\_polish1, 1.5V', '0719\_4round\_polish2, 1.5V', ...
     '0719\_5roll\_polish1, 1.5V', '0719\_5roll\_polish2, 1.5V', '0719\_5roll\_polish3, 1.7V', '0719\_5roll\_polish4, 1.7V', '0719\_5roll\_polish5, 1.5V', ...
     '0719\_6round\_polish1, 1.6V', '0719\_6round\_polish2, 1.6V','Location','northeastoutside')
 
%  text(200,15,'Electrolyte3','FontSize',14)
%  text(250,4,'Temp: 18-21C; bold line: fresh electrolyte1; -. line: no meaning','FontSize',14)
grid on
grid minor
set(gca,'fontsize',15)
xlabel('Voltage, V','FontSize',18);
ylabel('Current, mA','FontSize',18);
axis([200,700,-0.01,0.01])