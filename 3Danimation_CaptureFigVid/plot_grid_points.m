x = zeros(size(u));

quiver3(x,x,x,u,v,w)
rotate3d on
pbaspect([1 1 0.5])
title('uniform\_points\_1000', 'FontSize', 18)

%%
x_2 = zeros(size(u_2));

quiver3(x_2,x_2,x_2,u_2,v_2,w_2)
rotate3d on
pbaspect([1 1 0.5])
title('uniform\_points\_300', 'FontSize', 18)

%%
figure()
X_1 = zeros(size(U_1));
quiver3(X_1,X_1,X_1,U_1, V_1, W_1)
rotate3d on
pbaspect([1 1 0.5])
title('grid\_points\_all', 'FontSize', 18)

%%
figure()
X_2 = zeros(size(U_2));
quiver3(X_2,X_2,X_2,U_2, V_2, W_2)
rotate3d on
pbaspect([1 0.5 1])
title('grid\_points', 'FontSize', 18)






