function h = heightUsingTruncatedConeModel(area1, area2, v)
% ############################################################################
% Given two surface areas and the volume, calculate heigt as if the
% geometry is truncated cone
% ############################################################################
% ----------------------- load debug data -----------------------
% area1 = 16;
% area2 = 23.309;
% v = 25.9274;
% ---------------------------------------------------------------
r1 = sqrt(area1/pi);
r2 = sqrt(area2/pi);

h = 3*v/(pi*(r1*r2 + r1^2 + r2^2));
end


