% function delg = AAToG(theta, u)
% %     u_cross = u.'*u;
% %     u_LeviCivita = [0,-u(3),u(2); u(3),0,-u(1);-u(2),u(1),0];
% %     delg = cos(theta)*eye(3) + sin(theta)*u_LeviCivita + (1-cos(theta))*u_cross;
% end
function delg = AAToG (phi, x)
    delg = zeros(3);
	delg(1,1)=cos(phi)+(1.0-cos(phi))*(x(1)^2);
	delg(1,2)=x(1)*x(2)*(1.0-cos(phi))-x(3)*sin(phi);
	delg(1,3)=x(1)*x(3)*(1.0-cos(phi))+x(2)*sin(phi);
	delg(2,1)=x(1)*x(2)*(1.0-cos(phi))+x(3)*sin(phi);
	delg(2,2)=cos(phi)+(1.0-cos(phi))*(x(2)^2);
	delg(2,3)=x(3)*x(2)*(1.0-cos(phi))-x(1)*sin(phi);
	delg(3,1)=x(1)*x(3)*(1.0-cos(phi))-x(2)*sin(phi);
	delg(3,2)=x(2)*x(3)*(1.0-cos(phi))+x(1)*sin(phi);
	delg(3,3)=cos(phi)+(1.0-cos(phi))*(x(3)^2);
end

