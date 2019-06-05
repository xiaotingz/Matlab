% ############################################################################
% assume phi is given in degrees
% ############################################################################
function dg = AAToG(phi, x)
    phi = deg2rad(phi);
    x = x/norm(x);
    
    dg = zeros(3);
	dg(1,1)=cos(phi)+(1.0-cos(phi))*(x(1)^2);
	dg(1,2)=x(1)*x(2)*(1.0-cos(phi))-x(3)*sin(phi);
	dg(1,3)=x(1)*x(3)*(1.0-cos(phi))+x(2)*sin(phi);
	dg(2,1)=x(1)*x(2)*(1.0-cos(phi))+x(3)*sin(phi);
	dg(2,2)=cos(phi)+(1.0-cos(phi))*(x(2)^2);
	dg(2,3)=x(3)*x(2)*(1.0-cos(phi))-x(1)*sin(phi);
	dg(3,1)=x(1)*x(3)*(1.0-cos(phi))-x(2)*sin(phi);
	dg(3,2)=x(2)*x(3)*(1.0-cos(phi))+x(1)*sin(phi);
	dg(3,3)=cos(phi)+(1.0-cos(phi))*(x(3)^2);
    
    dg = dg';
end

