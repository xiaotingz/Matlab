function h = heightUsingTruncatedConeModel(area1, area2, v)
% ############################################################################
% Given two surface areas and the volume, calculate heigt as if the
% geometry is truncated cone
% See http://mathworld.wolfram.com/PyramidalFrustum.html
% ############################################################################
% ----------------------- load debug data -----------------------
% area1 = 16;
% area2 = 23.309;
% v = 25.9274;
% ---------------------------------------------------------------
    h = 3*v/(area1 + area2 + sqrt(area1*area2));
end


