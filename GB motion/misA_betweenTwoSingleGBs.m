
O = CrysSym; 

% phi1_1 = rad2deg(4.206);
% PHI_1 = rad2deg(0.866);
% phi2_1 = rad2deg(2.252);
% Z1_1 = [cosd(phi1_1),sind(phi1_1),0;-sind(phi1_1),cosd(phi1_1),0;0,0,1];
% X_1 = [1,0,0;0,cosd(PHI_1),sind(PHI_1);0,-sind(PHI_1),cosd(PHI_1)];
% Z2_1 = [cosd(phi2_1),sind(phi2_1),0;-sind(phi2_1),cosd(phi2_1),0;0,0,1];
% g1 = Z2_1*X_1*Z1_1;
%            
% phi1_2 = rad2deg(5.589);
% PHI_2 = rad2deg(0.168);
% phi2_2 = rad2deg(6.201);
% Z1_2 = [cosd(phi1_2),sind(phi1_2),0;-sind(phi1_2),cosd(phi1_2),0;0,0,1];
% X_2 = [1,0,0;0,cosd(PHI_2),sind(PHI_2);0,-sind(PHI_2),cosd(PHI_2)];
% Z2_2 = [cosd(phi2_2),sind(phi2_2),0;-sind(phi2_2),cosd(phi2_2),0;0,0,1];
% g2 = Z2_2*X_2*Z1_2;
	
% GB_see1 = 22.9@[-4,16,-5]
% GB_see2 = 38.3@[-7,6,1]
% GB_see3 = 51.6@[-2,6,9]   
angle1 = deg2rad(43);
axis1 = [1,0,0]; axis1 = axis1/norm(axis1);
angle2 = deg2rad(46);			
axis2 = [1,0,0] ; axis2 = axis2/norm(axis2);
g1 = AAToG(angle1, axis1);
g2 = AAToG(angle2, axis2);

RF1 = axis1 * tan(angle1/2);
RF1inFZ = ( (sum(RF1 >= 0) == 3) && (RF1(1) < (sqrt(2) - 1)) && ...
            (RF1(1) > RF1(2)) && (RF1(2) > RF1(3)) && (sum(RF1) < 1));
if RF1inFZ
    disp('R1 in FZ');
end


disA = pi/2;
for j = 1:24
    for k = 1:24
        gg1 = O(:,:,j) * g1;
        gg2 = O(:,:,k);

        dg = gg1 * gg2.';
        misA = acos(0.5*(trace(dg)-1));
        RFvec = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA)) *tan(misA/2);
%         if ( misA < disA && sum(misAxis > 0) == 3 )        
        axisInFZ = ( (sum(RFvec >= 0) == 3) && (RFvec(1) <= (sqrt(2) - 1)) && ...
                    (RFvec(1) >= RFvec(2)) && (RFvec(2) >= RFvec(3)) && (sum(RFvec) <= 1));
        if ( misA < disA && axisInFZ )
            g1 = dg;
            disA = misA;
            message = ['misA_1=', num2str(rad2deg(misA)), ' disA_1=', num2str(rad2deg(disA)), ' axis_1=',num2str(RFvec)];
            disp(message);
        end
        dg = gg2 * gg1.';
        misA = acos(0.5*(trace(dg)-1));
        RFvec = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA)) *tan(misA/2);
%         if ( misA < disA && sum(misAxis > 0) == 3 )
        axisInFZ = ( (sum(RFvec >= 0) == 3) && (RFvec(1) <= (sqrt(2) - 1)) && ...
                    (RFvec(1) >= RFvec(2)) && (RFvec(2) >= RFvec(3)) && (sum(RFvec) <= 1));
        if ( misA < disA && axisInFZ )
            g1 = dg;
            disA = misA;
            message = ['misA_1=', num2str(rad2deg(misA)), ' disA_1=', num2str(rad2deg(disA)), ' axis_2=',num2str(RFvec)];
            disp(message);
        end       
    end
end
disA = pi/2;
for j = 1:24
    for k = 1:24
        gg1 = O(:,:,j) * g2;
        gg2 = O(:,:,k);

        dg = gg1 * gg2.';
        misA = acos(0.5*(trace(dg)-1));
        RFvec = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA)) *tan(misA/2);
%         if ( misA < disA && sum(misAxis > 0) == 3 )
        axisInFZ = ( (sum(RFvec >= 0) == 3) && (RFvec(1) <= (sqrt(2) - 1)) && ...
                    (RFvec(1) >= RFvec(2)) && (RFvec(2) >= RFvec(3)) && (sum(RFvec) <= 1));
        if ( misA < disA && axisInFZ )
            g2 = dg;
            disA = misA;
            message = ['misA_2=', num2str(rad2deg(misA)), ' disA_2=', num2str(rad2deg(disA)), ' axis_2=',num2str(RFvec)];
            disp(message);
        end
        dg = gg2 * gg1.';
        misA = acos(0.5*(trace(dg)-1));
        RFvec = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA)) *tan(misA/2);
%         if ( misA < disA && sum(misAxis > 0) == 3 )
        axisInFZ = ( (sum(RFvec >= 0) == 3) && (RFvec(1) <= (sqrt(2) - 1)) && ...
                    (RFvec(1) >= RFvec(2)) && (RFvec(2) >= RFvec(3)) && (sum(RFvec) <= 1));
        if ( misA < disA && axisInFZ )
            g2 = dg;
            disA = misA;
            message = ['misA_2=', num2str(rad2deg(misA)), ' disA_2=', num2str(rad2deg(disA)), ' axis_2=',num2str(RFvec)];
            disp(message);
        end       
    end
end

% find the minimum misorientation angle
% minA = 10000;
% for k = 1:24
%     for l = 1:24
%         gg1 = O(:,:,k)*g1;
%         gg2 = O(:,:,l)*g2;
%         delg_sym = gg1*gg2.';
%         angle = acosd(0.5*(trace(delg_sym)-1));
%         if angle < minA
%             minA = angle;
%         end
%         delg_sym = gg2*gg1.';
%         angle = acosd(0.5*(trace(delg_sym)-1));
%         if angle < minA
%             minA = angle;
%         end
%     end
% end
% minA

% get the list of all misorientations
disA = pi/2;
misInfo = zeros(24*24*2,4);
for j = 1:24
    for k = 1:24
        gg1 = O(:,:,j) * g1;
        gg2 = O(:,:,k) * g2;

        dg = gg1 * gg2.';
        misA = acos(0.5*(trace(dg)-1));
        misAxis = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA));
        misInfo(2*((j-1)*24+k-1)+1, 1) = rad2deg(misA);
        misInfo(2*((j-1)*24+k-1)+1, 2:4) = misAxis;
%         (j-1)*24+k
        if ( misA < disA )
            disA = misA;          
        end
        dg = gg2 * gg1.';
        misA = acos(0.5*(trace(dg)-1));
        misAxis = [(dg(2,3)-dg(3,2)),(dg(3,1)-dg(1,3)),(dg(1,2)-dg(2,1))] /(2*sin(misA));
        misInfo(2*((j-1)*24+k), 1) = rad2deg(misA);
        misInfo(2*((j-1)*24+k), 2:4) = misAxis;
        if ( misA < disA )
            disA = misA;
        end      
    end
end
disA = rad2deg(disA)
misInfo_sorted = sortrows(misInfo);
