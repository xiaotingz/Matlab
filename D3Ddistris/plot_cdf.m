set(0,'defaultAxesFontSize',16)

% da = diff_tmp(:,1);
da = dist_left;

da_twin = da(istwin_an4);
da_notwin = da(~istwin_an4);

h1 = cdfplot(da_twin);
hold on
h2 = cdfplot(da_notwin);
xlabel('dist\_left');
set( h1, 'LineWidth', 3);
set( h2, 'LineWidth', 3);
xlim([-10, 10])

legend('twin boundaries', 'general boundaries', 'Location', 'SouthEast')


histogram(da_len_w_an4)


%% 
histogram(da_len_w_an4, 'FaceColor', [0.8, 0.8, 0.8], 'Normalization', 'probability')
xlabel('Weighted Dihedral Angle Of Triple Lines, °')
ylabel('Probability')

