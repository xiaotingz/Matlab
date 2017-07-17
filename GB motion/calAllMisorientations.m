numbin = 9;
step_phi1 = (pi/2)/numbin;
step_PHI = 1/numbin;
step_phi2 = (pi/2)/numbin;
EAs = zeros(numbin*numbin*numbin,3);

cnt = 1;
for i = 1:numbin
    for j = 1:numbin
        for k = 1:numbin
            EAs(cnt,1) = (i-0.5) * step_phi1;
            EAs(cnt,2) = asin((j-0.5) * step_PHI);
            EAs(cnt,3) = (k-0.5) * step_phi2;
            cnt = cnt + 1;
        end
    end
end

axis_raw = zeros(numbin*numbin*numbin,3);
angle = zeros(numbin*numbin*numbin,1);
for i = 1:length(EAs)
    phi1 = EAs(i,1);
    PHI = EAs(i,2);
    phi2 = EAs(i,3);
    Z1 = [cos(phi1),sin(phi1),0;-sin(phi1),cos(phi1),0;0,0,1];
    X = [1,0,0;0,cos(PHI),sin(PHI);0,-sin(PHI),cos(PHI)];
    Z2 = [cos(phi2),sin(phi2),0;-sin(phi2),cos(phi2),0;0,0,1];
    g = Z2*X*Z1;
    
    axis_tmp = [(g(2,3)-g(3,2)),(g(3,1)-g(1,3)),(g(1,2)-g(2,1))];
    axis_raw(i,:) = axis_tmp/norm(axis_tmp);
    angle(i) = acosd(0.5*(trace(g)-1));
end

cs = crystalSymmetry('m-3m');
axis = zeros(9*9*9,3);
for i = 1:length(axis_raw)
    tmp = Miller(axis_raw(i,1),axis_raw(i,2),axis_raw(i,3),cs,'hkl');
    tmp2 = round(tmp);
%     because of crystal symmetry, postive and negative is the same. Thus
%     use abs
    axis(i,:) = abs(int8(tmp2.hkl));
end

angle = int8(floor(angle/5) * 5);

axis_angle = [axis,angle];
axis_angle = unique(axis_angle,'rows');
