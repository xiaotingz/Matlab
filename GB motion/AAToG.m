function delg = AAToG(theta, u)
    u_cross = u.'*u;
    u_LeviCivita = [0,-u(3),u(2); u(3),0,-u(1);-u(2),u(1),0];
    delg = cos(theta)*eye(3) + sin(theta)*u_LeviCivita + (1-cos(theta))*u_cross;
end